#!/usr/bin/env python3
"""Multipopulation windowed PBSnj from the fly_fst.py pairwise Hudson FST windows.

This is a *window-level* extension of Schmidt et al. (2019) PBSnj: each window's
pairwise FST is a ratio of site numerator/denominator sums, converted to
T=-log(1-FST), and unrooted neighbor joining estimates terminal branches.
Original Schmidt et al. analyses used SNP-level estimates. No ancestral state
or root is inferred. Negative NJ branch lengths are retained for QC.
"""
import argparse
import bisect
import csv
import itertools
import json
import math
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np


def transformed_fst(fst, epsilon=1e-10):
    """Finite genetic distance; retain negative FST in original source tables."""
    return -math.log1p(-min(max(float(fst), 0.0), 1.0-epsilon))


def nj_terminal_lengths(distance, labels):
    """Saitou-Nei neighbor joining, return terminal branch lengths, no clipping.

    Works for n>=3. Deterministically breaks Q ties by input label order.
    A 3-taxon distance matrix yields the original three-population PBS exactly.
    """
    labels=list(labels)
    d=np.asarray(distance, dtype=float)
    n=len(labels)
    if n<3 or d.shape != (n,n) or not np.all(np.isfinite(d)):
        raise ValueError('NJ requires a finite square matrix with >=3 populations')
    if np.max(np.abs(d-d.T)) > 1e-8 or np.max(np.abs(np.diag(d)))>1e-8:
        raise ValueError('NJ distance matrix must be symmetric with zero diagonal')
    if np.min(d)<-1e-10: raise ValueError('NJ distances cannot be negative')
    dist={(i,j):float(d[i,j]) for i in range(n) for j in range(i+1,n)}
    active=list(range(n))
    next_id=n
    terminal={}
    def get(i,j):
        if i==j: return 0.0
        return dist[tuple(sorted((i,j)))]
    while len(active)>3:
        m=len(active)
        sums={i:sum(get(i,j) for j in active if i!=j) for i in active}
        i,j=min(itertools.combinations(active,2),key=lambda ij:
                ((m-2)*get(*ij)-sums[ij[0]]-sums[ij[1]],ij[0],ij[1]))
        delta=(sums[i]-sums[j])/(m-2)
        limb_i=(get(i,j)+delta)/2
        limb_j=(get(i,j)-delta)/2
        if i<n: terminal[labels[i]]=limb_i
        if j<n: terminal[labels[j]]=limb_j
        new_dist={k:(get(i,k)+get(j,k)-get(i,j))/2 for k in active if k not in (i,j)}
        active=[k for k in active if k not in (i,j)]
        for k,v in new_dist.items(): dist[tuple(sorted((k,next_id)))]=v
        next_id+=1
        active.append(next_id-1)
    i,j,k=active
    lengths={
        i:(get(i,j)+get(i,k)-get(j,k))/2,
        j:(get(i,j)+get(j,k)-get(i,k))/2,
        k:(get(i,k)+get(j,k)-get(i,j))/2
    }
    for node,val in lengths.items():
        if node<n: terminal[labels[node]]=val
    if len(terminal)!=n:
        raise RuntimeError('NJ internal error: not every leaf received a branch length')
    return terminal


def load_pair_windows(outdir, tag, arm, pop1, pop2):
    a,b=sorted((pop1,pop2))
    file=outdir/'windows'/tag/f'{a}_vs_{b}.{arm}.tsv'
    if not file.is_file(): raise FileNotFoundError(file)
    with file.open(newline='') as f:
        return { (int(r['start_0based']),int(r['end_0based_exclusive'])):r
                 for r in csv.DictReader(f, delimiter='\t') }


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore')
        w.writeheader();w.writerows(rows)


def calculate_window(pair_rows, pops, min_snps, min_call_fraction, epsilon):
    """Calculate pair distance matrix from all pairs for a shared genomic window."""
    matrix=np.zeros((len(pops),len(pops)),dtype=float)
    max_snps=0
    min_pair_snps=float('inf')
    min_pair_callable=float('inf')
    saturated=0
    fst_values={}
    min_pair_name=''
    for i,a in enumerate(pops):
        for j in range(i+1,len(pops)):
            b=pops[j]
            r=pair_rows[(a,b)]
            snps=int(r['n_snps'])
            callable_sites=int(r['n_callable_sites'])
            span=int(r['end_0based_exclusive'])-int(r['start_0based'])
            if callable_sites < min_pair_callable:
                min_pair_callable=callable_sites
                min_pair_name=f'{a}_vs_{b}'
            min_pair_snps=min(min_pair_snps,snps)
            max_snps=max(max_snps,snps)
            if span<=0 or callable_sites<span*min_call_fraction:
                return {'min_pair_callable_sites':min_pair_callable,
                        'lowest_callability_pair':min_pair_name,
                        'lowest_callability_fraction':callable_sites/span if span>0 else float('nan'),
                        'max_pair_snps':max_snps,'min_pair_snps':0,'n_saturated_pairs':saturated},'low_callability'
            f=float(r['fst'])
            # An invariant well-called population pair has FST undefined (0/0),
            # but contributes zero evolutionary distance in this calculation.
            if snps==0 and callable_sites>0 and not np.isfinite(f):f=0.0
            if not np.isfinite(f):return None,'undefined_fst'
            if f>=1:saturated+=1
            matrix[i,j]=matrix[j,i]=transformed_fst(f,epsilon)
            fst_values[(a,b)]=f
    if max_snps<min_snps:return None,'low_polymorphism'
    branches=nj_terminal_lengths(matrix,pops)
    return {'branches':branches,'lowest_callability_pair':min_pair_name,
            'lowest_callability_fraction':min_pair_callable/span, 'max_pair_snps':max_snps,
            'min_pair_snps':min_pair_snps,'min_pair_callable_sites':min_pair_callable,
            'n_saturated_pairs':saturated,'distance_matrix':matrix},'pass'


def annotate_windows(rows, genes_file, gene_format='auto', gene_id_attr='ID'):
    """Add gene overlaps and nearest-to-midpoint gene, with 0-based half-open coordinates.

    Gene IDs are taken from the GFF3 ID attribute (typically FBgn IDs), and
    descriptions are taken from available FlyBase GFF3 attributes.
    """
    from gene_set_enrichment import read_genes, parse_gff_attributes, open_text

    genes = read_genes(genes_file, gene_format, gene_id_attr)
    notes = {}
    if gene_format == 'gff3' or (gene_format == 'auto' and
        str(genes_file).lower().endswith(('.gff', '.gff3', '.gff.gz', '.gff3.gz'))):
        with open_text(genes_file) as fh:
            for ln in fh:
                if ln.startswith('#') or not ln.strip():
                    continue
                a = ln.rstrip('\n').split('\t')
                if len(a) < 9 or a[2] != 'gene':
                    continue
                attrs = parse_gff_attributes(a[8])
                gid = attrs.get(gene_id_attr) or attrs.get('ID')
                if gid:
                    desc = (attrs.get('description') or attrs.get('Description')
                            or attrs.get('Note') or '')
                    biotype = (attrs.get('gene_biotype') or attrs.get('gene_type')
                               or attrs.get('biotype') or '')
                    notes[gid] = (desc.replace('\t', ' ').replace('\n', ' ')[:240], biotype)
    by_arm = defaultdict(list)
    for arm, start, end, gid, name in genes:
        by_arm[arm].append((start, end, gid, name))
    indices = {}
    for arm, gg in by_arm.items():
        gg.sort()
        starts = [g[0] for g in gg]
        ends_sorted = sorted((g[1], g) for g in gg)
        indices[arm] = (gg, starts, ends_sorted, [v[0] for v in ends_sorted],
                        max((g[1]-g[0] for g in gg), default=0))

    fields = ['n_overlapping_genes', 'overlapping_gene_ids',
              'overlapping_gene_names', 'nearest_gene_id', 'nearest_gene_name',
              'nearest_gene_distance_bp', 'nearest_gene_description',
              'nearest_gene_biotype']
    for row in rows:
        arm = str(row['arm'])
        arm = arm[3:] if arm.startswith('chr') else arm
        start, end = int(row['start_0based']), int(row['end_0based_exclusive'])
        center = (start + end) / 2.0
        row.update({k: '' for k in fields})
        row['n_overlapping_genes'] = 0
        if arm not in indices:
            continue
        gg, starts, ends_sorted, sorted_ends, maxlen = indices[arm]
        hi = bisect.bisect_left(starts, end)
        lo = bisect.bisect_right(starts, start-maxlen)
        overlaps = [g for g in gg[lo:hi] if g[1] > start]
        row['n_overlapping_genes'] = len(overlaps)
        row['overlapping_gene_ids'] = ';'.join(g[2] for g in overlaps)
        row['overlapping_gene_names'] = ';'.join(g[3] for g in overlaps)

        # Find a gene containing or closest to the window midpoint. The nearest
        # gene may differ from a gene overlapping only the edge of the window.
        ci = bisect.bisect_right(starts, center)
        minlo = bisect.bisect_left(starts, center-maxlen)
        containing = [g for g in gg[minlo:ci] if g[1] > center]
        if containing:
            nearest = min(containing, key=lambda g: (g[1]-g[0], g[2]))
            distance = 0
        else:
            candidates = []
            if ci < len(gg):
                candidates.append(gg[ci])
            ei = bisect.bisect_right(sorted_ends, center)-1
            if ei >= 0:
                candidates.append(ends_sorted[ei][1])
            if not candidates:
                continue
            nearest = min(candidates, key=lambda g: (max(g[0]-center, center-g[1], 0), g[2]))
            distance = int(max(nearest[0]-center, center-nearest[1], 0))
        row['nearest_gene_id'], row['nearest_gene_name'] = nearest[2], nearest[3]
        row['nearest_gene_distance_bp'] = distance
        description, biotype = notes.get(nearest[2], ('', ''))
        row['nearest_gene_description'] = description
        row['nearest_gene_biotype'] = biotype
    return fields


def distinct_peak_labels(rows, focal, arm, cutoff, max_labels=5, min_separation_bp=100000):
    """Collapse overlapping top-quantile windows into genomic peak regions."""
    score_key = 'PBSnj_' + focal
    candidates = sorted((r for r in rows if r['arm'] == arm and
                         math.isfinite(float(r[score_key])) and
                         float(r[score_key]) >= cutoff),
                        key=lambda r: (int(r['start_0based']), int(r['end_0based_exclusive'])))
    clusters = []
    for row in candidates:
        start = int(row['start_0based'])
        end = int(row['end_0based_exclusive'])
        if not clusters or start > clusters[-1]['end']:
            clusters.append({'start': start, 'end': end, 'rows': [row]})
        else:
            clusters[-1]['end'] = max(clusters[-1]['end'], end)
            clusters[-1]['rows'].append(row)
    peaks = []
    for c in clusters:
        best = max(c['rows'], key=lambda r: float(r[score_key]))
        tied = [r for r in c['rows'] if float(r[score_key]) == float(best[score_key])]
        # Center the annotation on the maximum-scoring portion of the peak.
        centroid = sum((int(r['start_0based'])+int(r['end_0based_exclusive']))/2 for r in tied)/len(tied)
        name = best.get('nearest_gene_name', '') or best.get('nearest_gene_id', '')
        if not name:
            name = 'intergenic'
        peaks.append({'x': centroid, 'score': float(best[score_key]),
                      'label': name, 'start': c['start'], 'end': c['end']})
    peaks.sort(key=lambda a: (-a['score'], a['x']))
    selected = []
    for peak in peaks:
        if all(abs(peak['x']-other['x']) >= min_separation_bp for other in selected):
            selected.append(peak)
            if len(selected) >= max_labels:
                break
    return selected



def parse_gene_label_requests(gene_names='', gene_file=None):
    """Gene symbols or FBgn identifiers; comma-separated inline or one per line file."""
    requested=[]
    if gene_names:
        requested.extend(x.strip() for x in gene_names.replace('\n', ',').split(',') if x.strip())
    if gene_file:
        with Path(gene_file).open() as fh:
            for line in fh:
                line=line.split('#',1)[0].strip()
                if line:
                    requested.extend(x.strip() for x in line.replace('\t', ',').split(',') if x.strip())
    return list(dict.fromkeys(requested))


def custom_gene_labels(rows, focal, arm, genes, requests):
    """Position each requested gene at its highest finite overlapping PBSnj window."""
    score_key='PBSnj_' + focal
    by_name=defaultdict(list)
    for garm, start, end, gid, name in genes:
        if garm == arm:
            by_name[gid.lower()].append((start,end,gid,name))
            by_name[name.lower()].append((start,end,gid,name))
    eligible=[]
    for row in rows:
        if row['arm'] != arm:
            continue
        try:
            score=float(row[score_key])
        except (KeyError, ValueError, TypeError):
            continue
        if math.isfinite(score):
            eligible.append((int(row['start_0based']), int(row['end_0based_exclusive']), score))
    labels=[]
    unresolved=[]
    seen=set()
    for requested in requests:
        matches=by_name.get(requested.lower(),[])
        if not matches:
            unresolved.append((requested,'not found in annotation on this arm'))
            continue
        valid=[]
        for start,end,gid,name in matches:
            overlapping=[(ws,we,score) for ws,we,score in eligible if ws<end and we>start]
            if overlapping:
                ws,we,score=max(overlapping,key=lambda w:(w[2],-abs((w[0]+w[1])-(start+end))))
                valid.append((score,(ws+we)/2,gid,name))
        if not valid:
            unresolved.append((requested,'no passing PBSnj windows overlap gene'))
            continue
        score,x,gid,name=max(valid)
        if gid in seen:
            continue
        seen.add(gid)
        labels.append({'x':x,'score':score,'label':name or gid,'gene_id':gid,'manual':True})
    return labels, unresolved


def select_plot_labels(rows, focal, arm, cutoff, max_labels, min_separation_bp,
                       genes=None, requests=None, manual_only=False):
    """Custom annotations take priority; auto labels follow quantile/score threshold."""
    requested=requests or []
    custom, missing=custom_gene_labels(rows,focal,arm,genes or [],requested) if requested else ([],[])
    for gene,reason in missing:
        # A gene on another arm is expected; warn only about names absent from all arms outside here.
        if reason!='not found in annotation on this arm':
            print('[WARN] %s on %s: %s' % (gene,arm,reason),file=sys.stderr)
    automatic=[] if manual_only or max_labels==0 or not math.isfinite(cutoff) else distinct_peak_labels(
        rows,focal,arm,cutoff,max_labels,min_separation_bp)
    selected=list(custom)
    for peak in automatic:
        if any(abs(peak['x']-item['x']) < min_separation_bp for item in selected):
            continue
        selected.append(dict(peak,manual=False))
    return selected

def scan(outdir, focal='EA', populations=None, min_snps=10,
         min_call_fraction=0.5, epsilon=1e-10, quantile=0.99, do_plots=True,
         genes_file=None, gene_format='auto', gene_id_attr='ID',
         peak_labels=5, peak_min_separation_bp=100000,
         plot_label_quantile=0.99, plot_label_min_score=None,
         plot_label_genes='', plot_label_genes_file=None, plot_manual_only=False):
    if not (0 < plot_label_quantile < 1):
        raise ValueError('plot_label_quantile must be between 0 and 1')
    requests=parse_gene_label_requests(plot_label_genes,plot_label_genes_file)
    if requests and not genes_file:
        raise ValueError('Explicit plot gene labels require --genes / --pbsnj-genes or --enrich-genes')
    label_genes=[]
    if requests:
        from gene_set_enrichment import read_genes
        label_genes=read_genes(genes_file,gene_format,gene_id_attr)
        known=set(v.lower() for gene in label_genes for v in (gene[3],gene[4]))
        for name in requests:
            if name.lower() not in known:
                print('[WARN] Requested plot gene not found in annotations: '+name,file=sys.stderr)
    outdir=Path(outdir)
    cfg=json.loads((outdir/'run_config.json').read_text())
    available=sorted(cfg['populations'])
    pops=sorted(populations or available)
    if len(set(pops))!=len(pops) or not set(pops)<=set(available) or len(pops)<3:
        raise ValueError('Selected populations must contain >=3 distinct populations in run_config.json')
    focals=pops if focal=='all' else [focal]
    if not set(focals)<=set(pops):raise ValueError('Focal population must be in selected populations')
    resultdir=outdir/'pbsnj'
    allmeta=[]
    for window,step in cfg['window_sizes_and_steps']:
        tag=f'w{window}_s{step}'
        folder=outdir/'windows'/tag
        if not folder.is_dir():raise FileNotFoundError(folder)
        # Detect arms from the first pair, no assumptions about input chromosome set.
        a,b=pops[0],pops[1]
        arms=sorted(f.name.split('.')[-2] for f in folder.glob(f'{a}_vs_{b}.*.tsv'))
        if not arms:raise FileNotFoundError(f'No FST windows for pair {a}, {b} at {tag}')
        allrows=[]
        for arm in arms:
            pairs=list(itertools.combinations(pops,2))
            tables={pair:load_pair_windows(outdir,tag,arm,*pair) for pair in pairs}
            # Use exactly the shared coordinate windows across all pair tables.
            windows=sorted(set.intersection(*(set(t) for t in tables.values())))
            for start,end in windows:
                rowpairs={pair:table[(start,end)] for pair,table in tables.items()}
                result,status=calculate_window(rowpairs,pops,min_snps,min_call_fraction,epsilon)
                row={'arm':arm,'start_0based':start,'end_0based_exclusive':end,
                     'passes_filter':int(status=='pass'),'filter_reason':status,
                     'n_populations':len(pops),'populations':','.join(pops),
                     'n_pairs':len(pairs),
                     'max_pair_snps':result['max_pair_snps'] if result else '',
                     'min_pair_snps':result['min_pair_snps'] if result else '',
                     'min_pair_callable_sites':result['min_pair_callable_sites'] if result else '',
                     'n_saturated_pairs':result['n_saturated_pairs'] if result else '',
                     'lowest_callability_pair':result.get('lowest_callability_pair','') if result else '',
                     'lowest_callability_fraction':result.get('lowest_callability_fraction','') if result else '',
                     }
                for p in focals:
                    row[f'PBSnj_{p}']=result['branches'][p] if status=='pass' else float('nan')
                allrows.append(row)
        gene_fields = annotate_windows(allrows, genes_file, gene_format, gene_id_attr) if genes_file else []
        fields=['arm','start_0based','end_0based_exclusive','passes_filter','filter_reason',
                'n_populations','populations','n_pairs','max_pair_snps','min_pair_snps',
                'min_pair_callable_sites','n_saturated_pairs','lowest_callability_pair','lowest_callability_fraction']+[f'PBSnj_{p}' for p in focals]+gene_fields
        write_tsv(resultdir/tag/'pbsnj.tsv',fields,allrows)
        for p in focals:
            values=np.array([r[f'PBSnj_{p}'] for r in allrows],dtype=float)
            finite=values[np.isfinite(values)]
            cutoff=float(np.quantile(finite,quantile)) if len(finite) else float('nan')
            top=[r for r in allrows if np.isfinite(r[f'PBSnj_{p}']) and r[f'PBSnj_{p}']>=cutoff]
            write_tsv(resultdir/tag/f'{p}.top{(1-quantile)*100:g}pct.tsv',fields,top)
            allmeta.append({'tag':tag,'focal':p,'n_valid_windows':len(finite),'n_total_windows':len(allrows),
                            'quantile':quantile,'cutoff':cutoff,'n_outlier_windows':len(top)})
            if do_plots:
                label_cutoff=(float(plot_label_min_score) if plot_label_min_score is not None else
                              float(np.quantile(finite,plot_label_quantile)) if len(finite) else float('nan'))
                plot_scan(resultdir/tag/f'{p}.pbsnj.pdf',allrows,p,arms,label_cutoff,
                          peak_labels,peak_min_separation_bp,label_genes,requests,plot_manual_only)
    write_tsv(resultdir/'summary.tsv',
              ['tag','focal','n_valid_windows','n_total_windows','quantile','cutoff','n_outlier_windows'],allmeta)
    (resultdir/'methods.json').write_text(json.dumps({
        'method':'Windowed neighbor-joining PBSnj, Schmidt et al. 2019-inspired; not identical to SNP-level PBSnj',
        'citation':'Schmidt et al. PLoS Genetics 2019 DOI:10.1371/journal.pgen.1008485',
        'input':'per-window ratio-of-sums Hudson pairwise FST from fly_fst.py',
        'distance':'T=-ln(1-max(0,min(FST,1-epsilon)))',
        'algorithm':'Saitou-Nei neighbor joining for each window, unrooted',
        'negative_terminal_branches':'preserved',
        'invariant_pairs':'FST=0 if n_snps=0 and callable_sites>0',
        'missing_data':'all pairs must meet callability; undefined FST invalidates window',
        'focal':focal,'populations':pops,'gene_annotation_file':str(genes_file) if genes_file else None,
        'gene_annotation_note':'Overlaps are any gene intersecting a window; nearest gene is closest to window midpoint and not necessarily causal.',
        'peak_labels_per_arm':peak_labels, 'peak_min_separation_bp':peak_min_separation_bp,
        'plot_label_quantile':plot_label_quantile,'plot_label_min_score':plot_label_min_score,
        'plot_label_genes':requests,'plot_manual_only':plot_manual_only,
        'min_snps_max_any_pair':min_snps,
        'min_callable_fraction_per_pair':min_call_fraction,
        'eps':epsilon,'outlier_quantile':quantile,
        'saturated_fst_warning':'FST>=1 is capped at 1-epsilon; high PBSnj can depend on epsilon',
        'notes':'Outliers are descriptive, not selection p-values. Window NJ topology is estimated separately for every window. Differences in chromosome inversions, callability, and population history can cause apparent high PBSnj.'
    },indent=2)+'\n')
    print(f'[DONE] PBSnj ({", ".join(focals)}) in {resultdir}')


def plot_scan(path,rows,pop,arms,cutoff=float('nan'),max_labels=5,min_separation_bp=100000,
              genes=None,requests=None,manual_only=False):
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['pdf.fonttype'] = 42  # embed TrueType fonts for Illustrator editing
    matplotlib.rcParams['ps.fonttype'] = 42
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(len(arms),1,figsize=(13,max(3,2.35*len(arms))),squeeze=False)
    for ax,arm in zip(axes.flat,arms):
        selected=[r for r in rows if r['arm']==arm]
        x=np.array([(r['start_0based']+r['end_0based_exclusive'])/(2e6) for r in selected])
        y=np.array([r[f'PBSnj_{pop}'] for r in selected],dtype=float)
        ax.plot(x,y,linewidth=.85)
        ax.axhline(0,color='gray',lw=.6)
        ax.set_ylabel(f'{arm}\nPBSnj')
        ax.grid(alpha=.25)
        labels=select_plot_labels(rows,pop,arm,cutoff,max_labels,min_separation_bp,
                                  genes,requests,manual_only)
        for j, peak in enumerate(labels):
            ax.annotate(peak['label'],xy=(peak['x']/1e6,peak['score']),
                        xytext=(0,12+14*(j%3)),textcoords='offset points',
                        rotation=30,fontsize=7,ha='left',va='bottom',
                        arrowprops=dict(arrowstyle='-',color='0.45',lw=.65),
                        clip_on=False)
    axes[-1,0].set_xlabel('Position on chromosome arm (Mb)')
    fig.suptitle(f'Windowed PBSnj: {pop} terminal branch')
    fig.tight_layout(rect=(0,0,1,.96))
    path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path,dpi=180)
    plt.close(fig)


def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--outdir',required=True,type=Path,help='Completed fly_fst.py output directory')
    p.add_argument('--focal',default='EA',help='Focal population, or all (default EA)')
    p.add_argument('--populations',help='Comma-separated subset; default all populations')
    p.add_argument('--min-window-snps',type=int,default=10,help='At least this many informative SNPs in one pair')
    p.add_argument('--min-window-callable-fraction',type=float,default=0.5,help='All pairs must meet this fractional callability')
    p.add_argument('--outlier-quantile',type=float,default=.99)
    p.add_argument('--epsilon',type=float,default=1e-10)
    p.add_argument('--no-plots',action='store_true')
    p.add_argument('--genes',type=Path,help='GFF3/GFF3.gz or BED gene annotation for PBSnj tables and plots')
    p.add_argument('--gene-format',choices=['auto','gff3','bed'],default='auto')
    p.add_argument('--gene-id-attr',default='ID',help='GFF3 gene identifier attribute (usually ID)')
    p.add_argument('--peak-labels',type=int,default=5,help='Max distinct top-quantile peak labels per chromosome arm')
    p.add_argument('--peak-min-separation-bp',type=int,default=100000,help='Minimum separation between labelled peaks')
    p.add_argument('--plot-label-quantile',type=float,default=0.99,help='Auto-label peaks above this PBSnj percentile, independent of top1pct.tsv threshold')
    p.add_argument('--plot-label-min-score',type=float,help='Auto-label peaks with PBSnj >= this score; overrides plot-label-quantile')
    p.add_argument('--plot-label-genes',default='',help='Comma-separated gene symbols or FBgn IDs to label regardless of score')
    p.add_argument('--plot-label-genes-file',type=Path,help='Gene names or FBgn IDs, one per line')
    p.add_argument('--plot-manual-only',action='store_true',help='Only label requested genes; disable automatic peak labels')
    args=p.parse_args(argv)
    if not (0<=args.min_window_callable_fraction<=1) or args.min_window_snps<1 or not (0<args.outlier_quantile<1) or not (0<args.epsilon<1):
        p.error('Invalid numeric parameter')
    scan(args.outdir,args.focal,
         [v.strip() for v in args.populations.split(',')] if args.populations else None,
         args.min_window_snps,args.min_window_callable_fraction,args.epsilon,
         args.outlier_quantile,not args.no_plots,
         args.genes,args.gene_format,args.gene_id_attr,args.peak_labels,args.peak_min_separation_bp,
         args.plot_label_quantile,args.plot_label_min_score,args.plot_label_genes,
         args.plot_label_genes_file,args.plot_manual_only)

if __name__=='__main__':main()
