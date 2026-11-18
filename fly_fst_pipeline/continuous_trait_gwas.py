#!/usr/bin/env python3
"""Quantitative-trait association on coordinate-aligned haploid .seq files.

Each trait is analyzed separately. Association tests use generalized least squares
with null-model REML variance components (EMMAX/P3D approximation) and a leave-one-
chromosome-arm-out (LOCO) kinship computed from widely spaced polymorphic sites.

IMPORTANT: .seq files supply one allele per position, not diploid genotypes.
If they are consensus sequences from diploid genomes, heterozygous states are
irretrievably lost. Effect allele is lexicographically second allele, NOT REF/ALT.
Results are exploratory with small and geographically structured cohorts.
"""
import argparse
import bisect
import csv
import gzip
import json
import math
import sys
import heapq
from array import array
from collections import defaultdict
from contextlib import ExitStack
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import t as student_t
from scipy.stats import norm, rankdata

DM3 = {'2L':23011544, '2R':21146708, '3L':24543557, '3R':27905053, 'X':22422827}
BASE_CODES = np.frombuffer(b'ACGT', dtype=np.uint8)
BASE_NAMES = 'ACGT'


def open_tsv(path, mode='rt'):
    path = str(path)
    opener = gzip.open if path.endswith('.gz') else open
    return opener(path, mode, newline='')


def as_name(s):
    return str(s).strip()


def read_metadata(path):
    with open(path, newline='') as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames or not {'sample', 'population'} <= set(reader.fieldnames):
            raise ValueError('Metadata requires sample,population columns')
        rows = [(as_name(r['sample']), as_name(r['population'])) for r in reader]
    if len(rows) != len({s for s, p in rows}):
        raise ValueError('Duplicate metadata sample names')
    return rows


def load_trait_file(path, column, trait_columns, metadata, population_filter=None, match_suffix_n=False):
    """One numeric phenotype per sample; never mistake replicates for distinct genomes."""
    allowed = {s:p for s,p in metadata if population_filter is None or p in population_filter}
    with open(path, newline='') as f:
        r = csv.DictReader(f)
        if not r.fieldnames or column not in r.fieldnames:
            raise ValueError('Trait CSV lacks sample column {!r}'.format(column))
        if not trait_columns:
            raise ValueError('Specify --trait-columns explicitly (exclude replicate SD columns)')
        missing = set(trait_columns) - set(r.fieldnames)
        if missing:
            raise ValueError('Trait columns absent: {}'.format(sorted(missing)))
        seen = {}
        audit = []
        for row in r:
            supplied = as_name(row[column])
            matched = supplied if supplied in allowed else None
            reason = 'exact' if matched else 'unmatched'
            if matched is None and match_suffix_n:
                options = []
                if supplied.endswith('N'):
                    options.append(supplied[:-1])
                else:
                    options.append(supplied + 'N')
                valid = [s for s in options if s in allowed]
                if len(valid) == 1:
                    matched = valid[0]
                    reason = 'confirmed_suffix_N'
            if matched is not None and matched in seen:
                raise ValueError('Multiple rows match genotype {}: {} and {}. '
                                 'Aggregate technical replicates to a single mean first.'
                                 .format(matched, seen[matched], supplied))
            if matched is not None:
                seen[matched] = supplied
            values = {}
            for trait in trait_columns:
                val = as_name(row.get(trait, ''))
                try:
                    num = float(val) if val else float('nan')
                except ValueError:
                    raise ValueError('Non-numeric trait {} for {}: {!r}'.format(trait, supplied, val))
                values[trait] = num if np.isfinite(num) else float('nan')
            audit.append(dict(trait_sample=supplied, genotype_sample=matched or '',
                              match_status=reason, population=allowed.get(matched, ''),
                              traits=values))
    use = [a for a in audit if a['genotype_sample'] and
           any(np.isfinite(a['traits'][t]) for t in trait_columns)]
    if not use:
        raise ValueError('No finite traits match genotype metadata')
    return use, audit


def write_dicts(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_tsv(path, 'wt') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        for row in rows:
            w.writerow(row)


def file_sequence_length(path):
    """Allow only terminal newline/CRLF. Other whitespace is not silently removed."""
    path = Path(path)
    size = path.stat().st_size
    with path.open('rb') as f:
        if size >= 2:
            f.seek(-2, 2)
            last = f.read(2)
        else:
            f.seek(0)
            last = f.read()
    if last.endswith(b'\r\n'):
        return size - 2
    if last.endswith(b'\n'):
        return size - 1
    return size


def arm_files(seq_dir, samples, arm, skip_length=False):
    found = []
    for sample in samples:
        p = Path(seq_dir) / '{}_Chr{}.seq'.format(sample, arm)
        if p.is_file():
            found.append((sample, p))
    if found:
        lengths = {file_sequence_length(p) for _,p in found}
        if len(lengths) != 1:
            raise ValueError('{}: inconsistent nucleotide lengths {}. Check alignment.'.format(arm, lengths))
        length = next(iter(lengths))
        if not skip_length and length != DM3[arm]:
            raise ValueError('{}: {} bases != dm3 {}. Check reference assembly.'.format(arm, length, DM3[arm]))
    else:
        length = 0
    return found, length


def read_arm_chunks(files, length, chunk_size):
    """Yields sequence chunks; file handles stay open only during each pass."""
    if not files:
        return
    with ExitStack() as stack:
        handles = [stack.enter_context(open(p, 'rb')) for _,p in files]
        for offset in range(0, length, chunk_size):
            width = min(chunk_size, length-offset)
            parts = [h.read(width).upper() for h in handles]
            if any(len(d) != width for d in parts):
                raise ValueError('Unexpected EOF; sequences must be one character per site')
            mat = np.frombuffer(b''.join(parts), dtype=np.uint8).reshape(len(files), width)
            invalid = ~((mat==65)|(mat==67)|(mat==71)|(mat==84)|(mat==78))
            if invalid.any():
                i,j = np.argwhere(invalid)[0]
                raise ValueError('Invalid base {} in {} at position {}'.format(
                    chr(int(mat[i,j])), files[i][0], offset+j+1))
            yield offset, mat


def polymorphic_columns(mat):
    """Returns variant indices and coded two alleles from all rows in mat."""
    counts = np.stack([(mat == b).sum(axis=0) for b in BASE_CODES], axis=0)
    n_called = counts.sum(axis=0)
    alleles = counts > 0
    biallelic = (alleles.sum(axis=0)==2)
    idx = np.flatnonzero(biallelic)
    if len(idx)==0:
        return idx, np.empty(0, dtype=np.int8), np.empty(0, dtype=np.int8), np.empty(0, dtype=int)
    first = np.argmax(alleles[:, idx], axis=0)
    second = 3 - np.argmax(alleles[::-1, idx], axis=0)
    return idx, first, second, n_called[idx]


def estimate_kinship(seq_dir, samples, arms, chunk_size, step, skip_length, min_mac):
    """Pruned-by-position standardized haploid markers, mean-imputed missing calls.

    Kinship marker ascertainment is per arm across phenotype-matched individuals.
    No phenotype values are used for SNP selection. Missing arms contribute zero
    for those lines, a limitation described in the methods file.
    """
    n = len(samples)
    pos = {s:i for i,s in enumerate(samples)}
    sums = {}
    counts = {}
    lengths = {}
    for arm in arms:
        files, length = arm_files(seq_dir, samples, arm, skip_length)
        lengths[arm] = length
        accum = np.zeros((n,n), dtype=np.float64)
        used = 0
        if not files:
            sums[arm], counts[arm] = accum, used
            continue
        ix = np.array([pos[s] for s,p in files])
        for offset, mat in read_arm_chunks(files, length, chunk_size):
            w = mat.shape[1]
            sel = np.arange((-offset) % step, w, step, dtype=int)
            if len(sel) == 0:
                continue
            samp = mat[:, sel]
            idx, first, second, called = polymorphic_columns(samp)
            if not len(idx):
                continue
            # >90% callable in the arm's sequenced phenotype lines
            alt = (samp[:,idx] == BASE_CODES[second][None,:])
            alt_count = alt.sum(axis=0)
            mac = np.minimum(alt_count, called-alt_count)
            good = (called >= max(4, math.ceil(0.9*len(files)))) & (mac >= min_mac)
            if not good.any():
                continue
            idx = idx[good]
            alt = alt[:,good].astype(np.float64)
            called = called[good]
            alt_count = alt.sum(axis=0)
            freq = alt_count / called
            valid = np.isin(samp[:,idx], BASE_CODES)
            z = (alt - freq[None,:]) / np.sqrt(freq*(1-freq))[None,:]
            z[~valid] = 0
            full = np.zeros((n, z.shape[1]), dtype=np.float64)
            full[ix,:] = z
            accum += full @ full.T
            used += z.shape[1]
        sums[arm], counts[arm] = accum, used
        print('[GWAS] Kinship arm {}: {:,} spaced SNPs'.format(arm, used), file=sys.stderr, flush=True)
    return sums, counts, lengths


def make_design(populations, population_fixed=False):
    n = len(populations)
    columns = [np.ones(n)]
    covar_names = ['intercept']
    if population_fixed:
        groups = sorted(set(populations))
        for name in groups[1:]:
            columns.append(np.array([float(p==name) for p in populations]))
            covar_names.append('population_' + name)
    x = np.column_stack(columns)
    return x, covar_names


def transform_trait(y, how):
    y = np.asarray(y, dtype=np.float64)
    if how == 'none':
        return y
    if how == 'log':
        if np.any(y<=0):
            raise ValueError('Cannot log-transform zero/negative trait values')
        return np.log(y)
    if how == 'rank-int':
        return norm.ppf((rankdata(y, method='average')-0.5)/len(y))
    raise ValueError(how)


def fit_reml_delta(K, y, X):
    vals, U = np.linalg.eigh(K)
    vals = np.maximum(vals, 0.)
    yt = U.T @ y
    Xt = U.T @ X
    n,k = X.shape

    def reml(log_delta):
        v = vals + np.exp(log_delta)
        w = 1.0 / v
        C = Xt.T @ (w[:,None]*Xt)
        sign, logdetC = np.linalg.slogdet(C)
        if sign<=0:
            return np.inf
        try:
            beta = np.linalg.solve(C, Xt.T @ (w*yt))
        except np.linalg.LinAlgError:
            return np.inf
        res = yt-Xt@beta
        rss = float(np.dot(res*w, res))
        if rss <= 0:
            return np.inf
        return 0.5*((n-k)*np.log(rss/(n-k))+np.log(v).sum()+logdetC)

    opt = minimize_scalar(reml, bounds=(-9, 9), method='bounded',
                          options={'xatol': 0.002})
    delta = float(np.exp(opt.x))
    invsqrt = 1/np.sqrt(vals+delta)
    W = (U * invsqrt[None,:]) @ U.T
    return W, delta


def fit_null_model(y, pops, K, model, population_fixed=False, transform='none'):
    y = transform_trait(y, transform)
    X, covar_names = make_design(pops, population_fixed)
    n,k = X.shape
    if n <= k+3:
        raise ValueError('Not enough residual degrees of freedom: {} samples, {} covariates'.format(n,k))
    if model=='lmm':
        if K is None or not np.all(np.isfinite(K)) or np.trace(K) <= 1e-9:
            raise ValueError('LOCO kinship unavailable; add other arms or use --gwas-model linear')
        K = K * (n / np.trace(K))
        W, delta = fit_reml_delta(K, y, X)
    else:
        W, delta = np.eye(n), None
    yw = W @ y
    Xw = W @ X
    Q, R = np.linalg.qr(Xw, mode='reduced')
    if np.min(np.abs(np.diag(R))) < 1e-9:
        raise ValueError('Fixed-effects design has collinear columns')
    ry = yw - Q @ (Q.T @ yw)
    return dict(W=W, Q=Q, ry=ry, df=n-k-1, n=n, delta=delta,
                n_covariates=k, covariates=covar_names, trait_transform=transform)


def association_batch(model, G, minimum_called, max_missing, min_maf, min_mac):
    """G is 0/1/NaN where 1 codes the lexicographically second nucleotide."""
    n = G.shape[0]
    call = np.isfinite(G)
    called = call.sum(axis=0)
    alt_count = np.nansum(G, axis=0)
    mac = np.minimum(alt_count, called-alt_count)
    maf = np.divide(mac, called, out=np.zeros_like(mac,dtype=float), where=called>0)
    passing = ((called>=minimum_called)&(called >= math.ceil((1-max_missing)*n))&
               (maf>=min_maf)&(mac>=min_mac))
    if not np.any(passing):
        return passing, None
    gv = G[:,passing].copy()
    freq = alt_count[passing]/called[passing]
    miss = ~np.isfinite(gv)
    if miss.any():
        gv[miss] = np.broadcast_to(freq[None,:], gv.shape)[miss]
    W,Q,ry=model['W'],model['Q'],model['ry']
    gw = W@gv
    gp = gw-Q@(Q.T@gw)
    ssx = np.einsum('ij,ij->j', gp, gp)
    numerator = gp.T@ry
    sumy = float(ry@ry)
    valid = ssx>1e-10
    beta = np.divide(numerator, ssx, out=np.full(len(ssx), np.nan), where=valid)
    rss = np.maximum(sumy-beta*numerator, 0)
    se = np.sqrt(np.divide(rss, model['df']*ssx,
                           out=np.full(len(ssx),np.nan),where=valid))
    statistic = np.divide(beta,se,out=np.full(len(ssx),np.nan),where=se>0)
    p = 2*student_t.sf(abs(statistic),model['df'])
    p[~valid]=np.nan
    return passing, dict(beta=beta,se=se,t=statistic,p=p,maf=maf[passing],
                         mac=mac[passing],n_called=called[passing],gp=gp,ssx=ssx)


def make_null_permutations(model, n_perm, rng):
    if n_perm <= 0:
        return None
    ry = model['ry']
    Q = model['Q']
    perms = np.column_stack([ry[rng.permutation(len(ry))] for _ in range(n_perm)])
    return perms-Q@(Q.T@perms)


def update_maxT(permutation_resid, stats, model, maxT):
    if permutation_resid is None or stats is None:
        return
    gp, ssx = stats['gp'], stats['ssx']
    score = gp.T @ permutation_resid
    ry2 = np.einsum('ij,ij->j', permutation_resid, permutation_resid)
    rss = np.maximum(ry2[None,:] - score*score/ssx[:,None], 1e-12)
    denom = np.sqrt(rss*ssx[:,None] / model['df'])
    tests = np.abs(score / denom)
    maxT[:] = np.maximum(maxT, np.nanmax(tests, axis=0))


def bh_adjust(p):
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    ok = np.flatnonzero(np.isfinite(p))
    if not len(ok):
        return out
    order = ok[np.argsort(p[ok], kind='mergesort')]
    ranks = np.arange(1,len(order)+1,dtype=float)
    q = np.minimum.accumulate((p[order]*len(order)/ranks)[::-1])[::-1]
    out[order] = np.clip(q,0,1)
    return out


def annotate_top_variants(rows, gene_file):
    if not gene_file:
        for row in rows:
            row.update(overlapping_gene_ids='',overlapping_gene_names='',nearest_gene_id='',
                       nearest_gene_name='',nearest_gene_distance_bp='')
        return
    from gene_set_enrichment import read_genes
    genes = read_genes(gene_file)
    by_arm=defaultdict(list)
    for arm,start,end,gid,name in genes:
        by_arm[arm].append((start,end,gid,name))
    for r in rows:
        arr = by_arm.get(r['arm'], [])
        pos = int(r['pos_1based'])-1
        overlaps=[g for g in arr if g[0]<=pos<g[1]]
        r['overlapping_gene_ids']=';'.join(g[2] for g in overlaps)
        r['overlapping_gene_names']=';'.join(g[3] for g in overlaps)
        if arr:
            near=min(arr,key=lambda g:(max(g[0]-pos,pos-g[1]+1,0),g[1]-g[0]))
            r['nearest_gene_id']=near[2]
            r['nearest_gene_name']=near[3]
            r['nearest_gene_distance_bp']=max(near[0]-pos,pos-near[1]+1,0)
        else:
            r['nearest_gene_id']=r['nearest_gene_name']=r['nearest_gene_distance_bp']=''


def pbsnj_intervals(pbsnj_tsv, focal):
    """Read passed PBSnj windows; only finite focal branch statistics are used."""
    if not pbsnj_tsv or not Path(pbsnj_tsv).is_file():
        return {}
    by_arm=defaultdict(list)
    with open_tsv(pbsnj_tsv) as f:
        reader=csv.DictReader(f,delimiter='\t')
        field='PBSnj_'+focal
        if not reader.fieldnames or field not in reader.fieldnames:
            raise ValueError('{} lacks {}'.format(pbsnj_tsv,field))
        for r in reader:
            if r.get('passes_filter','1') != '1':
                continue
            try: value=float(r[field])
            except (ValueError,KeyError):continue
            if not np.isfinite(value):continue
            arm=r['arm']
            if arm.startswith('chr'):arm=arm[3:]
            by_arm[arm].append((int(r['start_0based']),int(r['end_0based_exclusive']),value))
    return {arm:sorted(rows) for arm,rows in by_arm.items()}


def annotate_pbsnj_top(rows, pbsnj_tsv=None, pbsnj_top_tsv=None, focal='EA'):
    all_intervals=pbsnj_intervals(pbsnj_tsv,focal)
    top_intervals=pbsnj_intervals(pbsnj_top_tsv,focal)
    meta={}
    for tag,d in [('all',all_intervals),('top',top_intervals)]:
        for arm,items in d.items():
            meta[(tag,arm)]=(items,[x[0] for x in items],max(x[1]-x[0] for x in items))
    for r in rows:
        arm=r['arm'];site=int(r['pos_1based'])-1
        values=[];peaks=[]
        for key,dest in [('all',values),('top',peaks)]:
            info=meta.get((key,arm))
            if not info:continue
            items,starts,maxspan=info
            start=bisect.bisect_left(starts,site-maxspan)
            end=bisect.bisect_right(starts,site)
            dest.extend(value for left,right,value in items[start:end] if left<=site<right)
        r['max_overlapping_PBSnj_'+focal]='{:.8g}'.format(max(values)) if values else ''
        r['overlaps_'+focal+'_top1pct']=int(bool(peaks)) if pbsnj_top_tsv else ''
        r['n_top1pct_PBSnj_windows']=len(peaks) if pbsnj_top_tsv else ''


BASE_FIELDS=['arm','pos_1based','allele0','allele1','n_samples','n_called',
             'minor_allele_count','maf','beta_allele1','se_beta','t','p','delta_lmm']
ADDITIONAL=['q_BH_trait','q_BH_all_traits','p_maxT_arm_FWER',
            'overlapping_gene_ids','overlapping_gene_names','nearest_gene_id',
            'nearest_gene_name','nearest_gene_distance_bp']


def add_qvalues_to_file(path, q_trait, q_global, q_offset, maxT, n_perm):
    path = Path(path)
    tmp = path.with_name(path.stem + '.q.tmp.gz')
    with open_tsv(path) as inf, open_tsv(tmp,'wt') as outf:
        reader = csv.DictReader(inf,delimiter='\t')
        writer = csv.DictWriter(outf, fieldnames=BASE_FIELDS + ['q_BH_trait','q_BH_all_traits','p_maxT_arm_FWER'],delimiter='\t')
        writer.writeheader()
        for i, row in enumerate(reader):
            row['q_BH_trait']='{:.8g}'.format(q_trait[i])
            row['q_BH_all_traits']='{:.8g}'.format(q_global[q_offset+i])
            if maxT is not None:
                t=abs(float(row['t']))
                row['p_maxT_arm_FWER']='{:.8g}'.format((1+np.count_nonzero(maxT>=t))/(n_perm+1))
            else:
                row['p_maxT_arm_FWER']=''
            writer.writerow(row)
    tmp.replace(path)


def generate_plots(outdir, trait, rows, n_tests, max_points=50000):
    try:
        import matplotlib
        matplotlib.use('Agg')
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42
        import matplotlib.pyplot as plt
    except ImportError as exc:
        print('[GWAS] Could not generate plots for {}: {}'.format(trait, exc), file=sys.stderr)
        return
    if not rows:
        return
    # Render at most one high-signal variant per 5-kb interval to keep PDFs editable.
    by_bin={}
    for r in rows:
        arm=r['arm']; pos=int(r['pos_1based']); p=max(float(r['p']),1e-300)
        k=(arm,pos//5000)
        if k not in by_bin or p < by_bin[k][2]:
            by_bin[k]=(arm,pos,p)
    arm_order={a:i for i,a in enumerate(DM3)}
    sorted_bins=sorted(by_bin.values(), key=lambda v:(arm_order.get(v[0],9),v[1]))
    if len(sorted_bins)>max_points:
        top=sorted(sorted_bins,key=lambda v:v[2])[:max_points//4]
        others=sorted_bins[::max(1,len(sorted_bins)//(max_points-len(top)))]
        sorted_bins=list({(a,p): (a,p,v) for a,p,v in top+others}.values())
    offsets={}; start=0
    for arm in DM3:
        offsets[arm]=start; start+=DM3[arm]
    fig,ax=plt.subplots(figsize=(13,4.5))
    colors=['#355c7d','#6c5b7b']
    for arm in DM3:
        arr=[r for r in sorted_bins if r[0]==arm]
        if not arr:continue
        ax.scatter([(offsets[arm]+r[1])/1e6 for r in arr],
                   [-np.log10(max(r[2],1e-300)) for r in arr],
                   s=2.5,color=colors[arm_order[arm]%2],rasterized=False,linewidths=0)
    axes=[(offsets[a]+DM3[a]/2)/1e6 for a in DM3]
    ax.set_xticks(axes);ax.set_xticklabels(list(DM3))
    ax.set_ylabel(r'$-\log_{10}(P)$');ax.set_xlabel('Chromosome arm')
    ax.set_title('{}: quantitative trait association (exploratory)'.format(trait))
    ax.grid(axis='y',alpha=0.2)
    fig.tight_layout()
    fig.savefig(outdir/'manhattan.pdf');plt.close(fig)
    # QQ from full p-values is created separately from input per-arm association tables.


def qq_plot(outdir, trait, all_p, max_points=10000):
    try:
        import matplotlib
        matplotlib.use('Agg')
        matplotlib.rcParams['pdf.fonttype']=42
        import matplotlib.pyplot as plt
    except ImportError:
        return
    pv=np.asarray(all_p,dtype=float)
    pv=np.sort(pv[np.isfinite(pv)])
    if not len(pv):return
    n=len(pv)
    inds=np.unique(np.linspace(0,n-1,min(n,max_points)).astype(int))
    expected=-np.log10((inds+0.5)/n)
    observed=-np.log10(np.maximum(pv[inds],1e-300))
    fig,ax=plt.subplots(figsize=(5.5,5.5))
    ax.scatter(expected,observed,s=3,color='#355c7d',rasterized=False)
    lim=max(np.max(expected),np.max(observed))*1.03
    ax.plot([0,lim],[0,lim],ls='--',lw=1,color='gray')
    ax.set_xlabel('Expected -log10(P)');ax.set_ylabel('Observed -log10(P)')
    ax.set_title('{}: QQ plot (n={:,} tests)'.format(trait,n))
    fig.tight_layout();fig.savefig(outdir/'qq.pdf');plt.close(fig)



def plot_genotype_phenotype_hits(seq_dir, trait, rows, pheno, outdir, max_hits=5):
    """Multipage vector PDF: raw phenotype by observed sequence allele at top hits."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        matplotlib.rcParams['pdf.fonttype']=42
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
    except ImportError:
        return
    if not rows or max_hits<1:
        return
    selected=[]
    for r in rows:
        # avoid near-duplicates from a single LD block
        if any(r['arm']==x['arm'] and abs(int(r['pos_1based'])-int(x['pos_1based']))<50000 for x in selected):
            continue
        selected.append(r)
        if len(selected)>=max_hits:break
    output=Path(outdir)/'genotype_trait_top_hits.pdf'
    with PdfPages(output) as pdf:
        for r in selected:
            allele0,allele1=r['allele0'],r['allele1']
            x0=[];x1=[]
            for sample,values in pheno.items():
                trait_value=values[trait]
                if not np.isfinite(trait_value):continue
                seq=Path(seq_dir)/'{}_Chr{}.seq'.format(sample,r['arm'])
                if not seq.exists():continue
                with seq.open('rb') as h:
                    h.seek(int(r['pos_1based'])-1)
                    base=h.read(1).decode('ascii','replace').upper()
                if base==allele0:x0.append(trait_value)
                if base==allele1:x1.append(trait_value)
            fig,ax=plt.subplots(figsize=(5.5,4.6))
            vals=[x0,x1]
            for k,v in enumerate(vals):
                if not v:continue
                jitter=np.linspace(-0.08,0.08,len(v)) if len(v)>1 else np.array([0.])
                ax.scatter(np.full(len(v),k)+jitter,v,s=26,alpha=.8,color=('#355c7d','#0b9488')[k])
                ax.hlines(np.mean(v),k-.17,k+.17,lw=2.5,color='black')
            ax.set_xticks([0,1]);ax.set_xticklabels(['{} (n={})'.format(allele0,len(x0)),
                                                    '{} (n={})'.format(allele1,len(x1))])
            ax.set_ylabel('{} (original measurement units)'.format(trait))
            ax.set_xlabel('Observed haploid .seq allele (not REF/ALT)')
            ax.set_title('{}:{}  |  P={:.2g}  |  beta={:.3g}'.format(
                r['arm'],r['pos_1based'],float(r['p']),float(r['beta_allele1'])))
            ax.grid(axis='y',alpha=.2)
            fig.tight_layout();pdf.savefig(fig);plt.close(fig)


def genes_from_top_variants(rows, focal='EA'):
    """Exploratory gene summary limited to top-N variant rows, not a gene-based test."""
    grouped={}
    for r in rows:
        ids=r.get('overlapping_gene_ids','').split(';') if r.get('overlapping_gene_ids') else []
        names=r.get('overlapping_gene_names','').split(';') if r.get('overlapping_gene_names') else []
        if not ids and r.get('nearest_gene_id'):
            ids=[r['nearest_gene_id']];names=[r['nearest_gene_name']]
        for i,gid in enumerate(ids):
            key=(gid,r['arm'])
            if key not in grouped:
                grouped[key]={'arm':r['arm'],'gene_id':gid,
                              'gene_name':names[i] if i<len(names) else gid,
                              'n_top_variants':0,'min_p':1.0,'best_position':'',
                              'min_q_BH_trait':1.0,'overlaps_top1pct_PBSnj':0}
            gene=grouped[key];gene['n_top_variants']+=1
            p=float(r['p'])
            if p<gene['min_p']:
                gene['min_p']=p;gene['best_position']=r['pos_1based']
            q=r.get('q_BH_trait','')
            if q!='':gene['min_q_BH_trait']=min(gene['min_q_BH_trait'],float(q))
            if str(r.get('overlaps_'+focal+'_top1pct',''))=='1':
                gene['overlaps_top1pct_PBSnj']=1
    return sorted(grouped.values(),key=lambda r:r['min_p'])

def run(seq_dir, metadata_file, trait_file, trait_columns, outdir, arms=None,
        trait_sample_column='sample', populations=None, model='lmm',
        min_maf=0.05, min_mac=3, min_samples=15, max_missing=0.2,
        kinship_step=5000, min_kinship_snps=20, permutations=0, seed=2026,
        chunk_size=100000, top_n=100, gene_file=None,
        population_covariate=False, trait_transform='none',
        skip_dm3_length_check=False, match_suffix_n=False,
        pbsnj_tsv=None, pbsnj_top_tsv=None, pbsnj_focal='EA',
        do_plots=True, plot_top_hits=5):
    if model not in ('lmm','linear') or trait_transform not in ('none','log','rank-int'):
        raise ValueError('Unknown model or transform')
    if not (0<=min_maf<0.5 and min_mac>=1 and min_samples>=5 and 0<=max_missing<1
            and kinship_step>=1 and chunk_size>=1 and permutations>=0 and top_n>=1):
        raise ValueError('Invalid GWAS filtering settings')
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    arms = arms or list(DM3)
    population_filter=set(populations) if populations else None
    meta=read_metadata(metadata_file)
    matched,audit=load_trait_file(trait_file,trait_sample_column,trait_columns,meta,
                                  population_filter,match_suffix_n)
    samples=[r['genotype_sample'] for r in matched]
    sample_pops={s:r['population'] for s,r in ((r['genotype_sample'],r) for r in matched)}
    pheno={r['genotype_sample']:r['traits'] for r in matched}
    matching=[]
    for r in audit:
        a={k:v for k,v in r.items() if k!='traits'}
        for trait in trait_columns:
            a['has_'+trait]=int(np.isfinite(r['traits'][trait]))
        matching.append(a)
    write_dicts(outdir/'sample_matching.tsv',
                ['trait_sample','genotype_sample','match_status','population']+
                ['has_'+t for t in trait_columns],matching)
    print('[GWAS] {} matched lines; {} separate traits: {}'.format(
        len(samples),len(trait_columns),','.join(trait_columns)),file=sys.stderr,flush=True)
    if model == 'lmm':
        kin_sums,kin_counts,lengths=estimate_kinship(
            seq_dir,samples,arms,chunk_size,kinship_step,skip_dm3_length_check,
            min_mac=min(2,min_mac))
    else:
        kin_sums={a:np.zeros((len(samples),len(samples))) for a in arms}
        kin_counts={a:0 for a in arms}
        lengths={a:arm_files(seq_dir,samples,a,skip_dm3_length_check)[1] for a in arms}
    totals=np.zeros((len(samples),len(samples)))
    total_snps=0
    for arm in arms:
        totals+=kin_sums[arm];total_snps+=kin_counts[arm]
    sample_index={s:i for i,s in enumerate(samples)}
    all_p={t:array('d') for t in trait_columns}
    result_paths={t:[] for t in trait_columns}
    counts={t:0 for t in trait_columns}
    arm_models={}; arm_maxT={}; status=[]
    top_cands={t:[] for t in trait_columns}  # bounded max-heap (-p, counter, row)
    plot_by_bin={t:{} for t in trait_columns}  # best SNP per 5kb for vector PDFs
    row_counter=0
    for arm in arms:
        files,length=arm_files(seq_dir,samples,arm,skip_dm3_length_check)
        if not files:
            print('[GWAS] {}: no matched sequences; skipped'.format(arm),file=sys.stderr)
            continue
        file_names=[s for s,p in files]
        arm_models[arm]={}
        for trait in trait_columns:
            subset=[i for i,s in enumerate(file_names) if np.isfinite(pheno[s][trait])]
            n=len(subset)
            record={'arm':arm,'trait':trait,'n_samples':n,'model':model,
                    'loco_kinship_markers':total_snps-kin_counts[arm],
                    'status':'','lmm_delta':''}
            if n<min_samples:
                record['status']='skipped_too_few_samples'
                status.append(record)
                print('[GWAS] {} {}: skipped, n={} < {}'.format(arm,trait,n,min_samples),file=sys.stderr)
                continue
            snps_loco=total_snps-kin_counts[arm]
            if model=='lmm' and snps_loco<min_kinship_snps:
                record['status']='skipped_insufficient_LOCO_kinship'
                status.append(record)
                print('[GWAS] {} {}: skipped, LOCO markers {} < {}'.format(
                    arm,trait,snps_loco,min_kinship_snps),file=sys.stderr)
                continue
            local_names=[file_names[i] for i in subset]
            y=np.array([pheno[s][trait] for s in local_names])
            p=[sample_pops[s] for s in local_names]
            global_ix=np.array([sample_index[s] for s in local_names])
            K=((totals-kin_sums[arm])/snps_loco)[np.ix_(global_ix,global_ix)] if model=='lmm' else None
            try:
                fm=fit_null_model(y,p,K,model,population_covariate,trait_transform)
            except ValueError as exc:
                record['status']='skipped_null_model: '+str(exc)
                status.append(record)
                print('[GWAS] {} {}: {}'.format(arm,trait,exc),file=sys.stderr)
                continue
            fm['subset']=np.array(subset)
            fm['n_sample_ids']=local_names
            fm['perm_resid']=make_null_permutations(fm,permutations,
                                                     np.random.RandomState(seed+100*list(DM3).index(arm)+list(trait_columns).index(trait)))
            fm['maxT']=np.zeros(permutations,dtype=float) if permutations else None
            fm['results_path']=outdir/trait/'associations.{}.tsv.gz'.format(arm)
            fm['tkey']=trait
            arm_models[arm][trait]=fm
            record['status']='tested'
            record['lmm_delta']=fm['delta'] if fm['delta'] is not None else ''
            status.append(record)
        if not arm_models[arm]:
            continue
        for trait,fm in arm_models[arm].items():
            fm['results_path'].parent.mkdir(parents=True,exist_ok=True)
            fm['fh']=open_tsv(fm['results_path'],'wt')
            fm['writer']=csv.DictWriter(fm['fh'], fieldnames=BASE_FIELDS,delimiter='\t')
            fm['writer'].writeheader()
            result_paths[trait].append((arm,fm['results_path'],len(all_p[trait])))
        print('[GWAS] Scanning {}: {:,} bp, {} traits'.format(
            arm,length,len(arm_models[arm])),file=sys.stderr,flush=True)
        try:
            for offset,mat in read_arm_chunks(files,length,chunk_size):
                idx, first, second, called=polymorphic_columns(mat)
                if not len(idx):continue
                all_alt=(mat[:,idx] == BASE_CODES[second][None,:]).astype(float)
                called_mask=np.isin(mat[:,idx],BASE_CODES)
                all_alt[~called_mask] = np.nan
                for trait,fm in arm_models[arm].items():
                    G=all_alt[fm['subset'],:]
                    passing, stats=association_batch(fm,G,min_samples,max_missing,min_maf,min_mac)
                    if stats is None:continue
                    valid=np.isfinite(stats['p'])
                    if not valid.any():continue
                    good_idx=np.flatnonzero(passing)[valid]
                    for key in ('beta','se','t','p','maf','mac','n_called'):
                        stats[key]=stats[key][valid]
                    if permutations:
                        substats={'gp':stats['gp'][:,valid],'ssx':stats['ssx'][valid]}
                        update_maxT(fm['perm_resid'],substats,fm,fm['maxT'])
                    for j,jidx in enumerate(good_idx):
                        pval=float(stats['p'][j]);global_pos=offset+int(idx[jidx])+1
                        row={'arm':arm,'pos_1based':global_pos,
                             'allele0':BASE_NAMES[int(first[jidx])],
                             'allele1':BASE_NAMES[int(second[jidx])],
                             'n_samples':fm['n'],'n_called':int(stats['n_called'][j]),
                             'minor_allele_count':int(stats['mac'][j]),
                             'maf':'{:.6g}'.format(stats['maf'][j]),
                             'beta_allele1':'{:.8g}'.format(stats['beta'][j]),
                             'se_beta':'{:.8g}'.format(stats['se'][j]),
                             't':'{:.8g}'.format(stats['t'][j]),
                             'p':'{:.10g}'.format(pval),
                             'delta_lmm':'{:.6g}'.format(fm['delta']) if fm['delta'] is not None else ''}
                        fm['writer'].writerow(row)
                        all_p[trait].append(pval)
                        counts[trait]+=1
                        # Bounded top-hit heap; no millions of row dictionaries retained.
                        row_counter += 1
                        key=(-pval, row_counter, row)
                        heap=top_cands[trait]
                        if len(heap)<top_n:
                            heapq.heappush(heap,key)
                        elif pval < -heap[0][0]:
                            heapq.heapreplace(heap,key)
                        bin_key=(arm,global_pos//5000)
                        prev=plot_by_bin[trait].get(bin_key)
                        if prev is None or pval<float(prev['p']):
                            plot_by_bin[trait][bin_key]=row
                if offset==0 or (offset//chunk_size)%100==0:
                    print('  [GWAS] {}: {:,}/{:,} bp'.format(arm,min(offset+mat.shape[1],length),length),
                          file=sys.stderr,flush=True)
        finally:
            for fm in arm_models[arm].values():
                fm['fh'].close()
                del fm['fh'];del fm['writer']
    write_dicts(outdir/'trait_arm_qc.tsv',
                ['arm','trait','n_samples','model','loco_kinship_markers','status','lmm_delta'],status)
    combined=np.concatenate([np.asarray(all_p[t],dtype=float) for t in trait_columns])
    q_global=bh_adjust(combined)
    global_offset=0
    summary=[]
    for trait in trait_columns:
        pvec=np.asarray(all_p[trait],dtype=float)
        qs=bh_adjust(pvec)
        # Output files are in the same arm order as accumulated P-values.
        for i,(arm,path,start) in enumerate(result_paths[trait]):
            end=(result_paths[trait][i+1][2] if i+1<len(result_paths[trait]) else len(pvec))
            fm=arm_models[arm][trait]
            add_qvalues_to_file(path,qs[start:end],q_global[global_offset+start:global_offset+end],
                                0,fm['maxT'],permutations)
        final_top=[entry[2] for entry in sorted(top_cands[trait], key=lambda x:-x[0])][:top_n]
        # Link top variants to their FDR values using a lookup from coordinate to index.
        index={}
        wanted={(r['arm'],int(r['pos_1based'])) for r in final_top}
        for arm,path,start in result_paths[trait]:
            if not any(k[0]==arm for k in wanted):
                continue
            with open_tsv(path) as f:
                for row in csv.DictReader(f,delimiter='\t'):
                    k=(row['arm'],int(row['pos_1based']))
                    if k in wanted:
                        index[k]=(row['q_BH_trait'],row['q_BH_all_traits'],row['p_maxT_arm_FWER'])
        for row in final_top:
            q=index[(row['arm'],int(row['pos_1based']))]
            row['q_BH_trait'],row['q_BH_all_traits'],row['p_maxT_arm_FWER']=q
        annotate_top_variants(final_top,gene_file)
        annotate_pbsnj_top(final_top,pbsnj_tsv,pbsnj_top_tsv,pbsnj_focal)
        folder=outdir/trait
        extra_fields=['max_overlapping_PBSnj_'+pbsnj_focal,
                      'overlaps_'+pbsnj_focal+'_top1pct','n_top1pct_PBSnj_windows']
        fields=BASE_FIELDS+ADDITIONAL+extra_fields
        write_dicts(folder/'top_variants.tsv', fields,final_top)
        if pbsnj_top_tsv:
            write_dicts(folder/'top_PBSnj_overlap.tsv',fields,
                        [r for r in final_top if r.get('overlaps_'+pbsnj_focal+'_top1pct')==1])
        gene_summary=genes_from_top_variants(final_top,pbsnj_focal)
        write_dicts(folder/'candidate_genes.tsv',
                    ['arm','gene_id','gene_name','n_top_variants','min_p','best_position',
                     'min_q_BH_trait','overlaps_top1pct_PBSnj'],gene_summary)
        if do_plots and len(pvec):
            qq_plot(folder,trait,pvec)
            generate_plots(folder,trait,list(plot_by_bin[trait].values()),len(pvec))
            plot_genotype_phenotype_hits(seq_dir,trait,final_top,pheno,folder,plot_top_hits)
        summary.append(dict(trait=trait,n_tests=len(pvec),minimum_p=float(pvec.min()) if len(pvec) else '',
                            n_q_below_0_05=int(np.sum(qs<0.05)),
                            n_global_q_below_0_05=int(np.sum(q_global[global_offset:global_offset+len(pvec)]<0.05)),
                            n_arms_tested=len(result_paths[trait])))
        global_offset+=len(pvec)
        print('[GWAS] {}: {:,} variants; min p {}; q<0.05: {}'.format(
            trait,len(pvec),summary[-1]['minimum_p'],summary[-1]['n_q_below_0_05']),
            file=sys.stderr,flush=True)
    write_dicts(outdir/'summary.tsv',
                ['trait','n_tests','minimum_p','n_q_below_0_05','n_global_q_below_0_05','n_arms_tested'], summary)
    config=dict(model=model,trait_columns=trait_columns,n_matched_lines=len(samples),
                per_trait_tests={s['trait']:s['n_tests'] for s in summary},
                min_mac=min_mac,min_maf=min_maf,min_samples=min_samples,max_missing=max_missing,
                kinship_step=kinship_step,kinship_markers_by_arm=kin_counts,
                min_kinship_snps=min_kinship_snps,permutations=permutations,
                population_covariate=population_covariate,trait_transform=trait_transform,
                pbsnj_focal=pbsnj_focal,
                pbsnj_tsv=str(pbsnj_tsv) if pbsnj_tsv else None,
                pbsnj_top_tsv=str(pbsnj_top_tsv) if pbsnj_top_tsv else None,
                plot_top_hits=plot_top_hits,do_plots=do_plots,
                sample_basis='one genotype and aggregated phenotype mean per genetic line',
                alleles='haploid .seq: allele1 is lexicographically second of the observed pair, NOT reference allele',
                mixed_model='REML null variance parameter, GLS Wald t test, LOCO thinned marker kinship (EMMAX/P3D approximation)',
                fdr='BH within each trait and across all trait-variant tests',
                optional_permutations='Freedman-Lane residual permutations after GLS whitening; arm-wide max-T family-wise correction, approximate',
                caveats=['Exploratory at small n; structure, linked sites and inversions may confound',
                         'Only assayed/phenotyped lines are used',
                         'Mean-imputation at missing SNP calls; effects are allelic contrasts of haploid sequence',
                         'LOCO kinship from thinned sites does not guarantee LD independence',
                         'X chromosome ploidy/sex and chromosome-arm availability require independent checks',
                         'maxT (when requested) is within chromosome arm and trait, not genomewide across arms or traits'])
    (outdir/'methods.json').write_text(json.dumps(config,indent=2))
    return summary


def make_parser():
    ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('--seq-dir',type=Path,required=True)
    ap.add_argument('--metadata',type=Path,required=True)
    ap.add_argument('--outdir',type=Path,required=True)
    ap.add_argument('--trait-file',type=Path,required=True,help='CSV with one phenotype row per genetic line')
    ap.add_argument('--trait-columns',required=True,help='Comma-separated numeric traits, tested separately')
    ap.add_argument('--trait-sample-column',default='sample')
    ap.add_argument('--arms',default='2L,2R,3L,3R,X')
    ap.add_argument('--populations',help='Comma-separated population subset')
    ap.add_argument('--gwas-model',choices=['lmm','linear'],default='lmm')
    ap.add_argument('--gwas-min-maf',type=float,default=0.05)
    ap.add_argument('--gwas-min-mac',type=int,default=3)
    ap.add_argument('--gwas-min-samples',type=int,default=15)
    ap.add_argument('--gwas-max-missing',type=float,default=0.2)
    ap.add_argument('--gwas-kinship-step',type=int,default=5000)
    ap.add_argument('--gwas-min-kinship-snps',type=int,default=20)
    ap.add_argument('--gwas-permutations',type=int,default=0)
    ap.add_argument('--gwas-seed',type=int,default=2026)
    ap.add_argument('--gwas-chunk-size',type=int,default=100000)
    ap.add_argument('--gwas-top-n',type=int,default=100)
    ap.add_argument('--gwas-plot-top-hits',type=int,default=5)
    ap.add_argument('--gwas-no-plots',action='store_true')
    ap.add_argument('--gwas-annotate-genes',type=Path)
    ap.add_argument('--gwas-population-covariate',action='store_true')
    ap.add_argument('--gwas-trait-transform',choices=['none','log','rank-int'],default='none')
    ap.add_argument('--gwas-pbsnj-tsv',type=Path,help='Optional PBSnj window table for annotating top associated variants')
    ap.add_argument('--gwas-pbsnj-top-tsv',type=Path,help='Optional focal top1pct PBSnj windows to flag overlap')
    ap.add_argument('--gwas-pbsnj-focal',default='EA')
    ap.add_argument('--gwas-match-suffix-n',action='store_true',help='Only for user-confirmed N suffix equivalences')
    ap.add_argument('--skip-dm3-length-check',action='store_true',help='Synthetic tests or validated non-dm3 inputs only')
    return ap


def main(argv=None):
    args=make_parser().parse_args(argv)
    traits=[s.strip() for s in args.trait_columns.split(',') if s.strip()]
    arms=[s.strip().replace('Chr','').replace('chr','') for s in args.arms.split(',')]
    if len(traits)!=len(set(traits)) or not traits:
        raise ValueError('Duplicate or empty trait list')
    if len(arms)!=len(set(arms)) or not set(arms)<=set(DM3):
        raise ValueError('Invalid --arms')
    pops=[x.strip() for x in args.populations.split(',')] if args.populations else None
    run(args.seq_dir,args.metadata,args.trait_file,traits,args.outdir,arms,
        args.trait_sample_column,pops,args.gwas_model,args.gwas_min_maf,
        args.gwas_min_mac,args.gwas_min_samples,args.gwas_max_missing,
        args.gwas_kinship_step,args.gwas_min_kinship_snps,args.gwas_permutations,
        args.gwas_seed,args.gwas_chunk_size,args.gwas_top_n,args.gwas_annotate_genes,
        args.gwas_population_covariate,args.gwas_trait_transform,
        args.skip_dm3_length_check,args.gwas_match_suffix_n,
        args.gwas_pbsnj_tsv,args.gwas_pbsnj_top_tsv,args.gwas_pbsnj_focal,
        not args.gwas_no_plots,args.gwas_plot_top_hits)


if __name__=='__main__':main()
