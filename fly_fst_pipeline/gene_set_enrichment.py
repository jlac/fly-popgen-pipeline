#!/usr/bin/env python3
"""EA-focused PBSnj gene-set enrichment using chromosome-arm block permutations.

Null: regional PBSnj signal is unrelated to fixed gene-set locations, conditional
on arm. Consecutive score blocks are shuffled *within chromosome arms*, retaining
within-block autocorrelation and the fixed gene annotation / gene lengths.

This does NOT fully control for recombination, inversions, population history,
functional annotation bias, or nonstationary genomic landscapes. Treat p-values
as exploratory positional-null p-values, not selection evidence by themselves.
"""
import argparse
import csv
import gzip
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

import numpy as np


def normalized_arm(arm):
    """Remove a leading chr prefix, compatible with Python 3.7+."""
    arm = str(arm)
    return arm[3:] if arm.startswith('chr') else arm


def open_text(path, mode='rt'):
    """Read ordinary or gzip-compressed text files transparently."""
    path = str(path)
    opener = gzip.open if path.lower().endswith('.gz') else open
    return opener(path, mode, encoding='utf-8')


def read_windows(path, focal):
    rows=[]
    with open_text(path) as h:
        reader=csv.DictReader(h,delimiter='\t')
        key='PBSnj_'+focal
        if key not in reader.fieldnames:
            raise ValueError(f'{key} not found in {path}; rerun pbsnj.py --focal {focal}')
        for r in reader:
            value=float(r[key]) if r[key] else float('nan')
            if r.get('passes_filter','1') != '1': value=float('nan')
            rows.append((normalized_arm(r['arm']),int(r['start_0based']),int(r['end_0based_exclusive']),value))
    if not rows: raise ValueError('No window data found')
    rows.sort(key=lambda x:(x[0],x[1],x[2]))
    return rows


def parse_gff_attributes(s):
    d={}
    for segment in s.split(';'):
        if '=' in segment:
            k,v=segment.split('=',1)
            d[k]=unquote(v)
    return d


def read_genes(path, fmt='auto', id_attr='ID'):
    path=Path(path)
    if fmt=='auto':
        bare_name = path.name.lower()[:-3] if path.name.lower().endswith('.gz') else path.name.lower()
        fmt = 'gff3' if bare_name.endswith(('.gff', '.gff3')) else 'bed'
    genes=[]
    if fmt=='gff3':
        with open_text(path) as h:
            for l in h:
                if not l.strip() or l.startswith('#'):continue
                a=l.rstrip('\n').split('\t')
                if len(a)<9 or a[2]!='gene':continue
                attr=parse_gff_attributes(a[8]); gene=attr.get(id_attr) or attr.get('ID')
                if not gene: continue
                genes.append((normalized_arm(a[0]),int(a[3])-1,int(a[4]),gene,attr.get('Name',gene)))
    else:
        with open_text(path) as h:
            for l in h:
                if not l.strip() or l.startswith('#'):continue
                a=l.rstrip('\n').split('\t')
                if len(a)<4:raise ValueError('Gene BED requires arm, start0, end0, gene_id columns')
                try: start,end=int(a[1]),int(a[2])
                except ValueError:
                    if a[1].lower() in ('start','start0','start_0based'):continue
                    raise
                genes.append((normalized_arm(a[0]),start,end,a[3],a[4] if len(a)>4 else a[3]))
    genes=[g for g in genes if g[1]>=0 and g[2]>g[1]]
    if not genes:raise ValueError('No valid gene intervals found')
    # Gene IDs must uniquely identify intervals. Alternative isoforms are excluded: supply gene features.
    ids=[g[3] for g in genes]
    if len(ids)!=len(set(ids)):raise ValueError('Duplicate gene_id; provide one gene-level interval per unique ID')
    return sorted(genes,key=lambda x:(x[0],x[1],x[2]))


def read_gene_sets(path):
    sets=defaultdict(set)
    with open_text(path) as h:
        for line in h:
            if not line.strip() or line.startswith('#'):continue
            a=line.rstrip('\n').split('\t')
            if str(path).lower().endswith(('.gmt', '.gmt.gz')):
                if len(a)<3:continue
                sets[a[0]].update(g.strip() for g in a[2:] if g.strip())
            else:
                if len(a)<2: raise ValueError('gene sets TSV: set_name TAB gene_id')
                if a[0].lower() in ('set','set_name','gene_set') and a[1].lower() in ('gene','gene_id'):continue
                sets[a[0]].add(a[1])
    if not sets:raise ValueError('No gene sets found')
    return dict(sets)


GO_NAMESPACES = {
    'biological_process': 'Biological Process',
    'molecular_function': 'Molecular Function',
    'cellular_component': 'Cellular Component',
}


def read_go_namespaces(obo_path):
    """Map canonical and alternate GO IDs to namespace from an OBO ontology."""
    namespace_by_id = {}
    current = None

    def commit(term):
        if not term or 'id' not in term or term.get('obsolete'):
            return
        namespace = term.get('namespace')
        if namespace not in GO_NAMESPACES:
            return
        namespace_by_id[term['id']] = namespace
        for alternate in term.get('alt_ids', []):
            namespace_by_id[alternate] = namespace

    with open_text(obo_path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if line == '[Term]':
                commit(current)
                current = {'alt_ids': []}
            elif line.startswith('['):
                commit(current)
                current = None
            elif current is not None:
                if line.startswith('id: '):
                    current['id'] = line[4:]
                elif line.startswith('namespace: '):
                    current['namespace'] = line[11:]
                elif line.startswith('alt_id: '):
                    current['alt_ids'].append(line[8:])
                elif line == 'is_obsolete: true':
                    current['obsolete'] = True
        commit(current)
    if not namespace_by_id:
        raise ValueError('No GO namespaces found in OBO: {}'.format(obo_path))
    return namespace_by_id


def plot_enrichment_results(results, output_path, focal, block_bp, title_prefix='GO gene-set enrichment'):
    if not results:
        return
    try:
        import matplotlib
        matplotlib.use('Agg')
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42
        import matplotlib.pyplot as plt
    except ImportError:
        return
    shown = sorted(results, key=lambda r: r['p_empirical_one_sided'])[:20]
    fig, ax = plt.subplots(figsize=(10, max(4, 0.39 * len(shown) + 1.5)))
    ypos = np.arange(len(shown))
    ax.barh(ypos, [-math.log10(max(r['p_empirical_one_sided'], 1e-300)) for r in shown])
    ax.set_yticks(ypos)
    ax.set_yticklabels([r['gene_set'] for r in shown])
    ax.invert_yaxis()
    ax.set_xlabel('-log10(one-sided empirical permutation p)')
    ax.set_title('{} | focal {} | block {:,} bp'.format(title_prefix, focal, block_bp))
    for y, row in enumerate(shown):
        ax.text(-0.01, y, 'q={:.3g}'.format(row['q_BH']), va='center', ha='right',
                transform=ax.get_yaxis_transform(), fontsize=8)
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches='tight')
    plt.close(fig)


def gene_window_pairs(windows,genes):
    """Flatten overlaps efficiently using sorted half-open genomic intervals."""
    byarm=defaultdict(list)
    for i,(arm,start,end,v) in enumerate(windows):byarm[arm].append((i,start,end))
    arrays={}
    for arm,items in byarm.items():
        arrays[arm]=(np.asarray([x[0] for x in items],int),
                     np.asarray([x[1] for x in items],int),
                     np.asarray([x[2] for x in items],int))
    linked=[];starts=[0]
    for arm,start,end,*_ in genes:
        if arm in arrays:
            indices,ss,ee=arrays[arm]
            # Starts are sorted; select candidate windows whose start precedes gene end.
            hi=int(np.searchsorted(ss,end,side='left'))
            lo=int(np.searchsorted(ee,start,side='right')) if np.all(ee[1:]>=ee[:-1]) else 0
            linked.extend(indices[lo:hi][ee[lo:hi]>start].tolist())
        starts.append(len(linked))
    return np.asarray(linked,dtype=np.int64),np.asarray(starts,dtype=np.int64)


def score_genes(window_scores,linked,starts,mode):
    scores=np.full(len(starts)-1,np.nan)
    if not len(linked):return scores
    lengths=np.diff(starts)
    valid=np.flatnonzero(lengths>0)
    if not len(valid):return scores
    x=window_scores[linked]
    if mode=='max':
        z=np.where(np.isfinite(x),x,-np.inf)
        v=np.maximum.reduceat(z,starts[valid])
        scores[valid]=np.where(np.isfinite(v),v,np.nan)
    else:
        z=np.where(np.isfinite(x),x,0.)
        c=np.isfinite(x).astype(int)
        sums=np.add.reduceat(z,starts[valid]); counts=np.add.reduceat(c,starts[valid])
        scores[valid]=np.divide(sums,counts,out=np.full(len(valid),np.nan),where=counts>0)
    return scores


def percentiles(scores):
    """Tie-aware 0..1 midrank percentiles among finite genes."""
    ans=np.full(len(scores),np.nan)
    ii=np.flatnonzero(np.isfinite(scores))
    if not len(ii):return ans
    v=scores[ii]; order=np.argsort(v,kind='stable'); sortedv=v[order]
    _,left,counts=np.unique(sortedv,return_index=True,return_counts=True)
    ranks=np.empty(len(v),dtype=float)
    for st,c in zip(left,counts):ranks[order[st:st+c]]=(st+(c-1)/2)/max(len(v)-1,1)
    ans[ii]=ranks
    return ans


def permute_blocks(values, windows, block_bp, rng):
    """Shuffle adjacent genomic blocks within each chromosome arm.

    Within each arm, only finite-score windows are shuffled in contiguous
    coordinate blocks. Nonfinite windows remain in place, preserving coverage.
    End blocks can differ in number of windows; concatenation preserves overall
    vector size but boundary alignment is approximate when lengths differ.
    """
    result=values.copy()
    byarm=defaultdict(lambda:defaultdict(list))
    for i,(arm,start,end,_) in enumerate(windows):
        if np.isfinite(values[i]):byarm[arm][start//block_bp].append(i)
    for arm,groups in byarm.items():
        blocks=[np.asarray(groups[k],dtype=int) for k in sorted(groups)]
        if len(blocks)<2: continue
        target=np.concatenate(blocks)
        order=rng.permutation(len(blocks))
        shuffled=np.concatenate([values[blocks[k]] for k in order])
        result[target]=shuffled
    return result


def bh_adjust(p):
    p=np.asarray(p,float); q=np.full(len(p),np.nan)
    idx=np.flatnonzero(np.isfinite(p)); o=idx[np.argsort(p[idx])]
    if len(o):q[o]=np.minimum.accumulate((p[o]*len(o)/np.arange(1,len(o)+1))[::-1])[::-1].clip(0,1)
    return q


def write_tsv(path,fields,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=fields,delimiter='\t',extrasaction='ignore')
        w.writeheader();w.writerows(rows)


def run(windows_file,genes_file,sets_file,outdir,focal='EA',permutations=1000,
        block_bp=250000,score='max',min_genes=5,seed=2026,gene_format='auto',id_attr='ID',go_obo=None):
    if permutations<1 or block_bp<1 or min_genes<1:raise ValueError('Invalid numerical option')
    windows=read_windows(windows_file,focal)
    genes=read_genes(genes_file,gene_format,id_attr)
    sets=read_gene_sets(str(sets_file))
    linked,starts=gene_window_pairs(windows,genes)
    vals=np.asarray([v for *_,v in windows],float)
    gene_ids=[g[3] for g in genes]; symbols=[g[4] for g in genes]
    indexes={s:i for i,s in enumerate(gene_ids)}
    gene_scores=score_genes(vals,linked,starts,score)
    ranks=percentiles(gene_scores)
    finite=np.isfinite(ranks)
    selected=[]
    excluded=[]
    for name,ids in sorted(sets.items()):
        hit=np.asarray(sorted({indexes[x] for x in ids if x in indexes and finite[indexes[x]]}),int)
        if len(hit)<min_genes:
            excluded.append({'gene_set':name,'n_input_genes':len(ids),'n_tested_genes':len(hit),'reason':'below_min_genes'})
        elif len(hit)==sum(finite):
            excluded.append({'gene_set':name,'n_input_genes':len(ids),'n_tested_genes':len(hit),'reason':'no_background_genes'})
        else:selected.append((name,ids,hit))
    if not selected:raise ValueError('No gene sets with >= min_genes mapped to scored genes; check annotation IDs and gene sets')
    # Difference in average within-genome percentile ranks, target vs other genes.
    observed=[]
    for name,ids,ix in selected:
        other=finite.copy();other[ix]=False
        observed.append(float(np.mean(ranks[ix])-np.mean(ranks[other])))
    extreme=np.zeros(len(selected),int)
    nullsum=np.zeros(len(selected),float);nullsq=np.zeros(len(selected),float)
    rng=np.random.default_rng(seed)
    for k in range(permutations):
        newvals=permute_blocks(vals,windows,block_bp,rng)
        nullrank=percentiles(score_genes(newvals,linked,starts,score))
        nullfinite=np.isfinite(nullrank)
        for j,(_,_,ix) in enumerate(selected):
            validix=ix[nullfinite[ix]]
            other=nullfinite.copy();other[ix]=False
            if not len(validix) or not other.any():continue
            stat=float(np.mean(nullrank[validix])-np.mean(nullrank[other]))
            if stat>=observed[j]-1e-12:extreme[j]+=1
            nullsum[j]+=stat; nullsq[j]+=stat*stat
    p=(extreme+1)/(permutations+1)
    q=bh_adjust(p)
    out=Path(outdir);out.mkdir(parents=True,exist_ok=True)
    results=[]
    for j,(name,ids,ix) in enumerate(selected):
        mean_null=nullsum[j]/permutations
        sd_null=max(0.0,nullsq[j]/permutations-mean_null**2)**.5
        results.append({'gene_set':name,'n_input_genes':len(ids),'n_tested_genes':len(ix),
                        'mean_set_percentile':float(np.mean(ranks[ix])),
                        'effect_rank_difference':observed[j],
                        'permutation_null_mean':mean_null,'permutation_null_sd':sd_null,
                        'n_perm_ge_observed':int(extreme[j]),'p_empirical_one_sided':float(p[j]),
                        'q_BH':float(q[j])})
    results.sort(key=lambda r:r['p_empirical_one_sided'])
    write_tsv(out/'enrichment.tsv',list(results[0]),results)
    write_tsv(out/'excluded_sets.tsv',['gene_set','n_input_genes','n_tested_genes','reason'],excluded)
    gtable=[]
    membership=defaultdict(list)
    for name,ids,_ in selected:
        for id in ids:membership[id].append(name)
    for g,v,rank in zip(genes,gene_scores,ranks):
        gtable.append({'arm':g[0],'start_0based':g[1],'end_0based_exclusive':g[2],
                       'gene_id':g[3],'gene_name':g[4],
                       'gene_length':g[2]-g[1], 'score':v,'score_percentile':rank,
                       'gene_sets':';'.join(membership[g[3]])})
    write_tsv(out/'gene_scores.tsv',list(gtable[0]),gtable)
    plot_enrichment_results(results, out/'enrichment.pdf', focal, block_bp)

    if go_obo is not None:
        ns_map = read_go_namespaces(go_obo)
        # GMT from flybase_gaf_to_gmt.py uses GO IDs as set names.
        classified = {ns: [] for ns in GO_NAMESPACES}
        excluded_by_ns = {ns: [] for ns in GO_NAMESPACES}
        unclassified = []
        for result in results:
            namespace = ns_map.get(result['gene_set'])
            if namespace is None:
                unclassified.append(result['gene_set'])
            else:
                classified[namespace].append(result)
        for record in excluded:
            namespace = ns_map.get(record['gene_set'])
            if namespace in excluded_by_ns:
                excluded_by_ns[namespace].append(record)
        if not any(classified.values()):
            raise ValueError('No GO identifiers in the GMT matched the supplied ontology. '
                             'Use the GO-ID-based flybase_all_go.gmt and its go-basic.obo.')
        for ns, label in GO_NAMESPACES.items():
            target = out / ns
            target.mkdir(parents=True, exist_ok=True)
            category_results = [dict(r) for r in classified[ns]]
            category_q = bh_adjust([r['p_empirical_one_sided'] for r in category_results])
            for r, adjusted_q in zip(category_results, category_q):
                r['q_BH_global'] = r['q_BH']
                r['q_BH'] = float(adjusted_q)
                r['go_namespace'] = ns
            fields = list(results[0]) + ['q_BH_global', 'go_namespace']
            write_tsv(target/'enrichment.tsv', fields, category_results)
            write_tsv(target/'excluded_sets.tsv',
                      ['gene_set','n_input_genes','n_tested_genes','reason'], excluded_by_ns[ns])
            plot_enrichment_results(category_results, target/'enrichment.pdf', focal, block_bp, label)
            (target/'methods.json').write_text(json.dumps({
                'go_namespace': ns, 'go_namespace_name': label,
                'go_obo': str(go_obo), 'gene_sets_tested': len(category_results),
                'permutations': permutations, 'block_bp': block_bp,
                'multiple_testing': 'Benjamini-Hochberg within this GO namespace',
                'q_BH_global': 'BH adjusted over all eligible gene sets, provided for comparison',
                'permutation_reuse': 'Same empirical p-values as combined enrichment; no second permutation run',
            }, indent=2))
            print('{}: {} sets -> {}'.format(label, len(category_results), target))
        (out/'unclassified_go_sets.txt').write_text('\n'.join(unclassified) + ('\n' if unclassified else ''))
        if unclassified:
            print('[WARN] {} sets were not found in the OBO; see unclassified_go_sets.txt'.format(len(unclassified)))

    (out/'methods.json').write_text(json.dumps({
        'analysis':'One-sided enrichment of elevated focal-population PBSnj scores in predeclared gene sets',
        'null':'Shuffle consecutive finite window-score blocks within each chromosome arm; missing positions and gene coordinates stay fixed',
        'gene_score':score,'gene_score_rule':'max or mean PBSnj among overlapping windows',
        'test_statistic':'mean genome-wide gene-score percentile(set) minus mean percentile(other eligible genes)',
        'p_value':'(1 + number of permutations with statistic >= observed) / (1 + permutations)',
        'multiple_testing':'Benjamini-Hochberg across eligible gene sets',
        'window_scores':str(windows_file),'gene_annotation':str(genes_file),'gene_sets':str(sets_file),
        'focal':focal,'permutations':permutations,'block_bp':block_bp,'seed':seed,
        'n_windows':len(windows),'n_finite_windows':int(np.isfinite(vals).sum()),
        'n_genes':len(genes),'n_scored_genes':int(finite.sum()),
        'n_tested_sets':len(selected),'n_excluded_sets':len(excluded),
        'go_obo':str(go_obo) if go_obo else None,
        'warning':'Exploratory positional null. Does not explicitly account for recombination, inversions, gene density, or local background selection. The finite-score coverage mask remains fixed. Use stratified controls and block-size sensitivity analyses.'
    },indent=2))
    print(f'Analyzed {len(selected)} gene sets, {int(finite.sum())} scored genes, {permutations} permutations')
    print(f'Output: {out}')
    return results


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pbsnj-tsv',required=True,help='e.g. fst_results/pbsnj/w50000_s10000/pbsnj.tsv')
    p.add_argument('--genes',required=True,help='BDGP5/dm3 gene BED (0-based) or GFF3 (1-based)')
    p.add_argument('--gene-format',choices=['auto','bed','gff3'],default='auto')
    p.add_argument('--gene-id-attribute',default='ID',help='GFF3 gene identifier attribute; match gene sets')
    p.add_argument('--gene-sets',required=True,help='Two-column TSV (set_name, gene_id) or GMT')
    p.add_argument('--outdir',required=True)
    p.add_argument('--focal',default='EA')
    p.add_argument('--permutations',type=int,default=1000)
    p.add_argument('--block-bp',type=int,default=250000)
    p.add_argument('--score',choices=['max','mean'],default='max')
    p.add_argument('--min-genes',type=int,default=5)
    p.add_argument('--seed',type=int,default=2026)
    p.add_argument('--go-obo',help='GO OBO for separate BP, MF, CC enrichment outputs and within-namespace BH correction')
    a=p.parse_args()
    run(a.pbsnj_tsv,a.genes,a.gene_sets,a.outdir,a.focal,a.permutations,a.block_bp,
        a.score,a.min_genes,a.seed,a.gene_format,a.gene_id_attribute,a.go_obo)

if __name__=='__main__':main()
