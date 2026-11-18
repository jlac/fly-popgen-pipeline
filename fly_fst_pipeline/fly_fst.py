#!/usr/bin/env python3
"""Population genetics from coordinate-aligned, one-character-per-base haploid .seq files.

Requires numpy; scipy and matplotlib for plots. Uses bounded-memory chunked I/O.
Each sample contributes ONE allele per site. N denotes missingness.
No VCF or reference FASTA required for within-dataset differentiation.
"""
import argparse, csv, gzip, itertools, json, math, re, sys
from collections import defaultdict
from pathlib import Path

import numpy as np

DM3 = {'2L':23011544, '2R':21146708, '3L':24543557,
       '3R':27905053, 'X':22422827}
BASES = 'ACGT'
BASE_CODES = np.array([ord(b) for b in BASES], dtype=np.uint8)


def sequence_base_length(path):
    """Count nucleotide bytes, ignoring only a terminal LF or CRLF (not internal gaps)."""
    size=path.stat().st_size
    with path.open('rb') as h:
        h.seek(max(size-2,0))
        trailing=h.read()
    if trailing.endswith(b'\r\n'): return size-2
    if trailing.endswith(b'\n'): return size-1
    return size


def parse_windows(value):
    out = []
    for elem in value.split(','):
        s, t = elem.strip().split(':')
        a, b = int(s), int(t)
        if a <= 0 or b <= 0 or a % b != 0:
            raise ValueError('Window sizes must be positive multiples of step sizes')
        out.append((a, b))
    if not out:
        raise ValueError('At least one window is required')
    return out


def write_table(path, header, rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with open(path,'w',newline='') as h:
        w=csv.writer(h,delimiter='\t'); w.writerow(header); w.writerows(rows)


def read_metadata(path):
    with open(path,newline='') as h:
        r=csv.DictReader(h)
        if not r.fieldnames or not {'sample','population'}.issubset(set(r.fieldnames)):
            raise ValueError('Metadata CSV must have columns sample,population')
        data=[]
        for row in r:
            sample=row['sample'].strip(); population=row['population'].strip()
            if not sample or not population or not re.fullmatch(r'[A-Za-z0-9_.-]+',sample):
                raise ValueError(f'Invalid metadata row: {row}')
            data.append((sample,population))
    if len(set(x[0] for x in data)) != len(data):
        raise ValueError('Duplicate sample names in metadata')
    if len(set(p for _,p in data)) < 2:
        raise ValueError('At least two populations are required')
    return data


def tajima_d(pi_total, S, n):
    if n < 3 or S < 2: return float('nan')
    a1=sum(1.0/i for i in range(1,n)); a2=sum(1.0/(i*i) for i in range(1,n))
    b1=(n+1)/(3*(n-1)); b2=2*(n*n+n+3)/(9*n*(n-1))
    c1=b1-1/a1; c2=b2-(n+2)/(a1*n)+a2/(a1*a1)
    e1=c1/a1; e2=c2/(a1*a1+a2)
    var=e1*S+e2*S*(S-1)
    return (pi_total-S/a1)/math.sqrt(var) if var>0 else float('nan')


def bincount_into(target, positions, values, resolution):
    if len(positions):
        target += np.bincount(positions//resolution, weights=values,
                              minlength=len(target))[:len(target)]


def window_rows(arrays, length, window, step, resolution):
    # Bin boundaries align with window boundaries. The terminal partial bin
    # contributes only to the shortened final window.
    pref={key:np.concatenate(([0.0],np.cumsum(v))) for key,v in arrays.items()}
    for start in range(0,length,step):
        end=min(length,start+window)
        a=start//resolution; b=(end+resolution-1)//resolution
        val={key:float(p[b]-p[a]) for key,p in pref.items()}
        yield start, end, val


def main(argv=None):
    ap=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument('--seq-dir',required=True,type=Path)
    ap.add_argument('--metadata',required=True,type=Path,help='CSV with sample,population')
    ap.add_argument('--outdir',required=True,type=Path)
    ap.add_argument('--arms',default='2L,2R,3L,3R,X',help='Comma-separated arms')
    ap.add_argument('--populations',help='Comma-separated subset of metadata populations for FST and downstream analyses; default all')
    ap.add_argument('--min-population-lines',type=int,default=2,help='Require this many sequenced lines per selected population/arm')
    ap.add_argument('--pbsnj-populations',help='Comma-separated subset of the populations included in the FST run; default all selected populations')
    ap.add_argument('--windows',default='5000:2500,10000:5000,50000:10000,100000:25000')
    ap.add_argument('--chunk-size',type=int,default=100000)
    ap.add_argument('--min-call-fraction',type=float,default=0.8)
    ap.add_argument('--min-called',type=int,default=2,help='Absolute min called haploid lines per population')
    ap.add_argument('--outlier-quantile',type=float,default=0.99,help='Descriptive upper-tail FST window threshold per pair and size')
    ap.add_argument('--min-window-snps',type=int,default=10,help='Minimum biallelic SNPs for outlier eligibility')
    ap.add_argument('--min-window-callable-fraction',type=float,default=0.5,help='Minimum fraction of window bases callable in both populations')
    ap.add_argument('--no-sites',action='store_true',help='Omit potentially huge per-SNP gzipped tables')
    ap.add_argument('--no-plots',action='store_true')
    ap.add_argument('--pbs-pbe',action='store_true',help='Compute Yi PBS and Yassin PBE after FST windows')
    ap.add_argument('--pbsnj',action='store_true',help='Compute windowed neighbor-joining PBSnj after FST windows')
    ap.add_argument('--pbsnj-focal',default='EA',help='Focal population for PBSnj, or all')
    ap.add_argument('--pbsnj-genes',type=Path,help='GFF3(.gz) or BED genes for PBSnj tables and peak labels; defaults to --enrich-genes')
    ap.add_argument('--pbsnj-peak-labels',type=int,default=5,help='Maximum top-quantile peaks labeled per arm')
    ap.add_argument('--pbsnj-peak-min-separation-bp',type=int,default=100000,help='Minimum distance between labeled peaks')
    ap.add_argument('--pbsnj-plot-label-quantile',type=float,default=0.99,help='Automatic plot label threshold as PBSnj quantile')
    ap.add_argument('--pbsnj-plot-label-min-score',type=float,help='Automatic label threshold as an absolute PBSnj score')
    ap.add_argument('--pbsnj-plot-label-genes',default='',help='Comma-separated gene symbols or FBgn identifiers to label')
    ap.add_argument('--pbsnj-plot-label-genes-file',type=Path,help='Requested gene names or FBgn IDs, one per line')
    ap.add_argument('--pbsnj-plot-manual-only',action='store_true',help='Label requested genes only; suppress automatic peak labels')
    ap.add_argument('--enrich-genes',type=Path,help='dm3 gene intervals BED or GFF3 for PBSnj gene-set enrichment (requires --pbsnj)')
    ap.add_argument('--enrich-gene-sets',type=Path,help='Two-column set/gene TSV or GMT gene sets (requires --enrich-genes)')
    ap.add_argument('--enrich-window',default='50000:10000',help='Window:step resolution for PBSnj enrichment')
    ap.add_argument('--enrich-permutations',type=int,default=1000)
    ap.add_argument('--enrich-block-bp',type=int,default=250000)
    ap.add_argument('--enrich-min-genes',type=int,default=5)
    ap.add_argument('--enrich-score',choices=['max','mean'],default='max')
    ap.add_argument('--enrich-seed',type=int,default=2026)
    ap.add_argument('--enrich-go-obo',type=Path,help='GO go-basic.obo: separate BP, MF, CC enrichment with within-ontology FDR')
    ap.add_argument('--branch-triplets',type=Path,help='CSV focal,nonfocal1,nonfocal2; default every 3-population orientation')
    ap.add_argument('--branch-populations',help='Comma-separated subset for --pbs-pbe, otherwise --populations selection')
    ap.add_argument('--branch-baseline',choices=['arm','genome'],default='arm',help='Median calibration scale for PBE')
    # Independent continuous-trait association module (one test per trait, per SNP).
    ap.add_argument('--continuous-trait-gwas',action='store_true',help='Run quantitative trait association after population statistics')
    ap.add_argument('--trait-file',type=Path,help='CSV with one row per phenotyped genetic line')
    ap.add_argument('--trait-columns',help='Comma-separated continuous traits; each is tested separately')
    ap.add_argument('--trait-sample-column',default='sample')
    ap.add_argument('--gwas-model',choices=['lmm','linear'],default='lmm')
    ap.add_argument('--gwas-min-maf',type=float,default=0.05)
    ap.add_argument('--gwas-min-mac',type=int,default=3)
    ap.add_argument('--gwas-min-samples',type=int,default=15)
    ap.add_argument('--gwas-max-missing',type=float,default=0.2)
    ap.add_argument('--gwas-kinship-step',type=int,default=5000)
    ap.add_argument('--gwas-min-kinship-snps',type=int,default=20)
    ap.add_argument('--gwas-permutations',type=int,default=0,help='Optional, expensive per-arm residual maxT permutations; default off')
    ap.add_argument('--gwas-seed',type=int,default=2026)
    ap.add_argument('--gwas-top-n',type=int,default=100)
    ap.add_argument('--gwas-plot-top-hits',type=int,default=5)
    ap.add_argument('--gwas-annotate-genes',type=Path,help='GFF3(.gz) gene annotations for top SNPs; falls back to --enrich-genes')
    ap.add_argument('--gwas-pbsnj-window',help='Optional window:step used for PBSnj overlap, e.g. 10000:5000')
    ap.add_argument('--gwas-population-covariate',action='store_true',help='Add fixed population indicators to genetic relatedness correction')
    ap.add_argument('--gwas-trait-transform',choices=['none','log','rank-int'],default='none')
    ap.add_argument('--gwas-match-suffix-n',action='store_true',help='Permit known same-line sample aliases with N suffix')
    ap.add_argument('--allow-incomplete-arms',action='store_true',help='Skip samples lacking a given arm (QC notes this)')
    ap.add_argument('--skip-dm3-length-check',action='store_true',help='Only for testing or validated non-dm3 input')
    args=ap.parse_args(argv)
    if (args.chunk_size<1 or args.min_population_lines<2 or not (0<args.min_call_fraction<=1) or args.min_called<2
        or not (0<args.outlier_quantile<1) or not (0<=args.min_window_callable_fraction<=1)
        or args.min_window_snps<1):
        ap.error('Invalid chunk size or minimum calls settings')
    if args.enrich_genes or args.enrich_gene_sets:
        if not (args.pbsnj and args.enrich_genes and args.enrich_gene_sets):
            ap.error('--enrich-genes and --enrich-gene-sets require --pbsnj and both enrichment paths')
        if args.pbsnj_focal=='all':ap.error('Set a single --pbsnj-focal for enrichment')
        if args.enrich_permutations<1 or args.enrich_block_bp<1:ap.error('Invalid enrichment settings')
    if args.continuous_trait_gwas and (not args.trait_file or not args.trait_columns):
        ap.error('--continuous-trait-gwas requires --trait-file and --trait-columns')
    if args.trait_file and not args.continuous_trait_gwas:
        ap.error('--trait-file requires --continuous-trait-gwas')
    metadata=read_metadata(args.metadata)
    available_pops={p for _,p in metadata}
    if args.populations:
        chosen=[x.strip() for x in args.populations.split(',') if x.strip()]
        if len(chosen)!=len(set(chosen)):ap.error('Duplicate names in --populations')
        unknown=set(chosen)-available_pops
        if unknown:ap.error(f'Unknown --populations: {sorted(unknown)}; available: {sorted(available_pops)}')
        metadata=[(s,p) for s,p in metadata if p in chosen]
    pops=sorted({p for _,p in metadata})
    if len(pops)<2:ap.error('At least 2 selected populations required for FST')
    from collections import Counter
    nper=Counter(p for _,p in metadata)
    sparse={p:nper[p] for p in pops if nper[p]<max(args.min_called,args.min_population_lines)}
    if sparse:ap.error(f'Insufficient samples in selected populations: {sparse}. '
                       'Remove them with --populations or add sequences; Hudson FST needs >=2 per population.')
    if args.pbsnj and len(pops)<3:ap.error('PBSnj requires >=3 populations')
    if args.pbsnj_populations and not args.pbsnj:ap.error('--pbsnj-populations requires --pbsnj')
    if args.branch_populations and not args.pbs_pbe:ap.error('--branch-populations requires --pbs-pbe')
    selected_branch=None
    if args.branch_populations:
        selected_branch=[x.strip() for x in args.branch_populations.split(',') if x.strip()]
        if len(selected_branch)!=len(set(selected_branch)) or len(selected_branch)<3 or not set(selected_branch)<=set(pops):
            ap.error('--branch-populations requires >=3 distinct populations in --populations')
    selected_pbsnj=None
    if args.pbsnj_populations:
        selected_pbsnj=[x.strip() for x in args.pbsnj_populations.split(',') if x.strip()]
        if len(selected_pbsnj)!=len(set(selected_pbsnj)) or not set(selected_pbsnj)<=set(pops) or len(selected_pbsnj)<3:
            ap.error('--pbsnj-populations requires >=3 distinct populations included in --populations')
        if args.pbsnj_focal!='all' and args.pbsnj_focal not in selected_pbsnj:
            ap.error('--pbsnj-focal must appear in --pbsnj-populations')
    if args.pbsnj and args.pbsnj_focal!='all' and args.pbsnj_focal not in pops:
        ap.error(f'--pbsnj-focal {args.pbsnj_focal} not in selected populations: {pops}')
    pairs=list(itertools.combinations(pops,2))
    windows=parse_windows(args.windows)
    resolution=math.gcd(*[v for win in windows for v in win])
    arms=[s.strip().replace('chr','').replace('Chr','') for s in args.arms.split(',')]
    if not all(arm in DM3 for arm in arms): ap.error('Unknown arm in --arms')
    args.outdir.mkdir(parents=True,exist_ok=True)
    bypop={p:[] for p in pops}
    for sample,pop in metadata: bypop[pop].append(sample)
    config={'populations':pops,'metadata_n':len(metadata),'pairs':len(pairs),
            'metadata_path':str(args.metadata),'selected_pbsnj_populations':selected_pbsnj,
            'haploid':True,'assembly_hypothesis':'dm3/BDGP5; length check does not verify alignment',
            'window_sizes_and_steps':windows,'resolution':resolution,'min_call_fraction':args.min_call_fraction,
            'min_called':args.min_called,'sitewise_fst':'Hudson 2-pop unbiased biallelic estimator',
            'window_fst':'sum(numerator)/sum(denominator); unweighted by SNP density',
            'tajima_d':'complete-case sites only; distinct from missingness-adjusted pi',
            'snp_outliers':'descriptive ranking only, not significance tests'}
    (args.outdir/'run_config.json').write_text(json.dumps(config,indent=2))
    qcrows=[]; pairtot=defaultdict(lambda:[0.0,0.0,0,0]); win_paths=[]
    all_arm_summaries=[]
    for arm in arms:
        print(f'[INFO] Processing {arm}',file=sys.stderr,flush=True)
        available=[]; absent=[]
        for sample,p in metadata:
            file=args.seq_dir/f'{sample}_Chr{arm}.seq'
            if file.exists(): available.append((sample,p,file))
            else: absent.append(sample)
        if absent and not args.allow_incomplete_arms:
            raise FileNotFoundError(f'{arm}: missing sequence(s): {absent[:15]} (use --allow-incomplete-arms to skip)')
        if not available:
            print(f'[WARN] No sequences for {arm}; skipped',file=sys.stderr); continue
        sizes={sequence_base_length(file) for _,_,file in available}
        if len(sizes)!=1:
            raise ValueError(f'{arm}: different file lengths: {sizes}; verify coordinate alignment')
        length=sizes.pop()
        if not args.skip_dm3_length_check and length!=DM3[arm]:
            raise ValueError(f'{arm}: length {length} != dm3 {DM3[arm]}. Do not force mismatched coordinates.')
        sample_to_pos={s:i for i,(s,_,_) in enumerate(available)}
        groups={p:np.array([sample_to_pos[s] for s in bypop[p] if s in sample_to_pos],dtype=int) for p in pops}
        npop={p:len(groups[p]) for p in pops}
        insufficient={p:npop[p] for p in pops if npop[p]<max(args.min_called,args.min_population_lines)}
        if insufficient:raise ValueError(f'{arm}: fewer than required sequenced samples per population: {insufficient}')
        thresholds={p:max(args.min_called,math.ceil(args.min_call_fraction*npop[p])) for p in pops}
        n_bins=(length+resolution-1)//resolution
        buffers={}
        for p in pops:
            for key in ('callable','pi_sum','complete_pi','complete_S','complete_sites'):
                buffers[(p,key)]=np.zeros(n_bins,dtype=np.float64)
        for a,b in pairs:
            for key in ('num','den','snps','callable','excluded_multiallelic'):
                buffers[((a,b),key)]=np.zeros(n_bins,dtype=np.float64)
        site_handles={}
        if not args.no_sites:
            for a,b in pairs:
                path=args.outdir/'snp_fst'/f'{a}_vs_{b}.{arm}.tsv.gz'
                path.parent.mkdir(parents=True,exist_ok=True)
                h=gzip.open(path,'wt',compresslevel=3)
                h.write('arm\tpos_1based\tallele0\tallele1\tpop1\tpop2\tn_called1\tn_called2\tfreq_allele1_pop1\tfreq_allele1_pop2\thudson_num\thudson_den\thudson_fst\n')
                site_handles[(a,b)]=h
        call_counts=np.zeros(len(available),dtype=np.int64)
        try:
            with __import__('contextlib').ExitStack() as stack:
                fh=[stack.enter_context(open(file,'rb')) for _,_,file in available]
                for offset in range(0,length,args.chunk_size):
                    width=min(args.chunk_size,length-offset)
                    data=[h.read(width).upper() for h in fh]
                    if any(len(x)!=width for x in data):raise ValueError('Unexpected EOF in sequence')
                    mat=np.frombuffer(b''.join(data),dtype=np.uint8).reshape(len(available),width)
                    valid=(mat==65)|(mat==67)|(mat==71)|(mat==84)
                    invalid=(~valid)&(mat!=78)
                    if invalid.any():
                        i,j=np.argwhere(invalid)[0]; raise ValueError(f'Invalid base in {available[i][0]} {arm}:{offset+j+1}: {chr(int(mat[i,j]))!r}')
                    call_counts+=valid.sum(axis=1)
                    locations=np.arange(offset,offset+width,dtype=np.int64)
                    counts={}; ncalled={}
                    for p in pops:
                        inds=groups[p]
                        counts[p]=np.stack([(mat[inds]==base).sum(axis=0) for base in BASE_CODES],axis=0) if len(inds) else np.zeros((4,width),dtype=np.int32)
                        ncalled[p]=counts[p].sum(axis=0)
                        ok=ncalled[p]>=thresholds[p]
                        bincount_into(buffers[(p,'callable')],locations,ok.astype(float),resolution)
                        # Unbiased pairwise nucleotide diversity from the sampled haploid alleles
                        n=ncalled[p].astype(float)
                        n_safe=np.maximum(n,1)
                        freq=counts[p]/n_safe
                        pi=(1.0-(freq*freq).sum(axis=0))*n_safe/np.maximum(n_safe-1,1)
                        pi[n<2]=0
                        pi[~ok]=0
                        bincount_into(buffers[(p,'pi_sum')],locations,pi,resolution)
                        complete=(ncalled[p]==npop[p]) & (npop[p]>=3)
                        seg=complete & ((counts[p]>0).sum(axis=0)>1)
                        bincount_into(buffers[(p,'complete_pi')],locations,np.where(complete,pi,0),resolution)
                        bincount_into(buffers[(p,'complete_S')],locations,seg.astype(float),resolution)
                        bincount_into(buffers[(p,'complete_sites')],locations,complete.astype(float),resolution)
                    for a,b in pairs:
                        if npop[a]<args.min_called or npop[b]<args.min_called:
                            continue
                        ca=counts[a];cb=counts[b];na=ncalled[a];nb=ncalled[b]
                        ok=(na>=thresholds[a])&(nb>=thresholds[b])
                        bincount_into(buffers[((a,b),'callable')],locations,ok.astype(float),resolution)
                        pooled=ca+cb
                        count_alleles=(pooled>0).sum(axis=0)
                        multi=ok&(count_alleles>2)
                        bincount_into(buffers[((a,b),'excluded_multiallelic')],locations,multi.astype(float),resolution)
                        selected=ok&(count_alleles==2)
                        idx=np.flatnonzero(selected)
                        if not len(idx):continue
                        is_present=pooled[:,idx]>0
                        first=np.argmax(is_present,axis=0)
                        second=3-np.argmax(is_present[::-1,:],axis=0)
                        pa=ca[second,idx]/na[idx]
                        pb=cb[second,idx]/nb[idx]
                        num=(pa-pb)**2-pa*(1-pa)/(na[idx]-1)-pb*(1-pb)/(nb[idx]-1)
                        den=pa*(1-pb)+pb*(1-pa)
                        good=den>0
                        idx=idx[good]; first=first[good];second=second[good]
                        pa=pa[good]; pb=pb[good];num=num[good];den=den[good]
                        bincount_into(buffers[((a,b),'num')],locations[idx],num,resolution)
                        bincount_into(buffers[((a,b),'den')],locations[idx],den,resolution)
                        bincount_into(buffers[((a,b),'snps')],locations[idx],np.ones(len(idx)),resolution)
                        if (a,b) in site_handles:
                            h=site_handles[(a,b)]
                            for k,j in enumerate(idx):
                                h.write(f'{arm}\t{offset+j+1}\t{BASES[first[k]]}\t{BASES[second[k]]}\t{a}\t{b}\t{na[j]}\t{nb[j]}\t{pa[k]:.6g}\t{pb[k]:.6g}\t{num[k]:.8g}\t{den[k]:.8g}\t{num[k]/den[k]:.8g}\n')
                    if offset==0 or ((offset//args.chunk_size)%100==0):
                        print(f'  {arm}: {min(offset+width,length):,}/{length:,} bp',file=sys.stderr,flush=True)
        finally:
            for h in site_handles.values():h.close()
        for (sample,p,_),n in zip(available,call_counts):
            qcrows.append((sample,p,arm,length,int(n),length-int(n),round(n/length,6)))
        for a,b in pairs:
            key=(a,b)
            v=[buffers[(key,z)].sum() for z in ('num','den','snps','callable')]
            for i,x in enumerate(v):pairtot[key][i]+=x
            all_arm_summaries.append((arm,a,b,int(npop[a]),int(npop[b]),int(v[2]),int(v[3]),
                                      v[0]/v[1] if v[1]>0 else float('nan')))
            for window,step in windows:
                tag=f'w{window}_s{step}'
                path=args.outdir/'windows'/tag/f'{a}_vs_{b}.{arm}.tsv'
                rows=[]
                for start,end,vals in window_rows({z:buffers[(key,z)] for z in ('num','den','snps','callable','excluded_multiallelic')},length,window,step,resolution):
                    den=vals['den'];fst=vals['num']/den if den>0 else float('nan')
                    rows.append((arm,start,end,a,b,round(fst,8),int(round(vals['snps'])),int(round(vals['callable'])),
                                 int(round(vals['excluded_multiallelic'])),round(vals['num'],8),round(den,8)))
                write_table(path,['arm','start_0based','end_0based_exclusive','pop1','pop2','fst','n_snps','n_callable_sites','n_excluded_multiallelic','sum_numerator','sum_denominator'],rows)
                win_paths.append(path)
        for p in pops:
            for window,step in windows:
                tag=f'w{window}_s{step}'
                path=args.outdir/'diversity'/tag/f'{p}.{arm}.tsv'
                rows=[]
                arr={k:buffers[(p,k)] for k in ('callable','pi_sum','complete_pi','complete_S','complete_sites')}
                for start,end,v in window_rows(arr,length,window,step,resolution):
                    callable_sites=int(round(v['callable']))
                    pi=v['pi_sum']/callable_sites if callable_sites else float('nan')
                    S=int(round(v['complete_S']))
                    d=tajima_d(v['complete_pi'],S,npop[p])
                    rows.append((arm,start,end,p,npop[p],callable_sites,round(pi,10),
                                 int(round(v['complete_sites'])),S,round(d,7)))
                write_table(path,['arm','start_0based','end_0based_exclusive','population','n_lines',
                                  'n_callable_sites','pi_per_callable_bp','n_complete_sites','segregating_complete_sites','tajima_d_complete_sites'],rows)
    write_table(args.outdir/'qc_by_sample.tsv',['sample','population','arm','length','called_bases','missing_bases','called_fraction'],qcrows)
    write_table(args.outdir/'fst_by_arm.tsv',['arm','pop1','pop2','n_lines1','n_lines2','n_biallelic_snps','n_callable_sites','hudson_fst'],all_arm_summaries)
    write_table(args.outdir/'fst_genomewide.tsv',['pop1','pop2','hudson_fst','n_biallelic_snps','n_callable_sites','sum_numerator','sum_denominator'],
                [(a,b,num/den if den>0 else float('nan'),int(snp),int(call),num,den)
                 for (a,b),(num,den,snp,call) in sorted(pairtot.items())])
    summarize_outliers(args.outdir,pairs,arms,windows,args.outlier_quantile,
                       args.min_window_snps,args.min_window_callable_fraction)
    if not args.no_plots:
        try: plot_results(args.outdir,pops,pairs,arms,windows)
        except ImportError as err: print(f'[WARN] Plotting unavailable: {err}',file=sys.stderr)
    if args.pbs_pbe:
        from branch_stats import run as run_branch_stats
        run_branch_stats(args.outdir, args.branch_triplets, args.branch_baseline,
                         args.min_window_snps,args.min_window_callable_fraction,do_plots=not args.no_plots,
                         populations=selected_branch)
    if args.pbsnj:
        from pbsnj import scan as run_pbsnj
        run_pbsnj(args.outdir, focal=args.pbsnj_focal, populations=selected_pbsnj, min_snps=args.min_window_snps,
                  min_call_fraction=args.min_window_callable_fraction,
                  quantile=args.outlier_quantile, do_plots=not args.no_plots,
                  genes_file=args.pbsnj_genes or args.enrich_genes,
                  peak_labels=args.pbsnj_peak_labels,
                  peak_min_separation_bp=args.pbsnj_peak_min_separation_bp,
                  plot_label_quantile=args.pbsnj_plot_label_quantile,
                  plot_label_min_score=args.pbsnj_plot_label_min_score,
                  plot_label_genes=args.pbsnj_plot_label_genes,
                  plot_label_genes_file=args.pbsnj_plot_label_genes_file,
                  plot_manual_only=args.pbsnj_plot_manual_only)
    if args.enrich_genes:
        from gene_set_enrichment import run as run_enrichment
        sizes=parse_windows(args.enrich_window)
        if len(sizes)!=1:ap.error('--enrich-window must be one window:step')
        win,step=sizes[0]
        tag=f'w{win}_s{step}'
        path=args.outdir/'pbsnj'/tag/'pbsnj.tsv'
        if not path.is_file():ap.error(f'Missing requested PBSnj window file: {path}')
        run_enrichment(path,args.enrich_genes,args.enrich_gene_sets,
                       args.outdir/'enrichment'/tag,args.pbsnj_focal,
                       args.enrich_permutations,args.enrich_block_bp,args.enrich_score,
                       args.enrich_min_genes,args.enrich_seed,go_obo=args.enrich_go_obo)
    if args.continuous_trait_gwas:
        pbsnj_win_path=None
        pbsnj_top_path=None
        if args.pbsnj and args.pbsnj_focal != 'all':
            gwas_window=args.gwas_pbsnj_window or (args.enrich_window if args.enrich_genes else
                         '{}:{}'.format(windows[0][0],windows[0][1]))
            dimensions=parse_windows(gwas_window)
            if len(dimensions)!=1:ap.error('--gwas-pbsnj-window needs one window:step')
            w,st=dimensions[0]
            base=args.outdir/'pbsnj'/'w{}_s{}'.format(w,st)
            pbsnj_win_path=base/'pbsnj.tsv'
            pbsnj_top_path=base/'{}.top1pct.tsv'.format(args.pbsnj_focal)
            if not pbsnj_win_path.exists():pbsnj_win_path=None
            if not pbsnj_top_path.exists():pbsnj_top_path=None
        from continuous_trait_gwas import run as run_continuous_trait_gwas
        traits=[v.strip() for v in args.trait_columns.split(',') if v.strip()]
        if not traits or len(traits)!=len(set(traits)):
            ap.error('--trait-columns must contain distinct nonempty names')
        run_continuous_trait_gwas(
            seq_dir=args.seq_dir, metadata_file=args.metadata, trait_file=args.trait_file,
            trait_columns=traits, outdir=args.outdir/'continuous_trait_gwas',
            arms=arms, trait_sample_column=args.trait_sample_column,
            populations=pops, model=args.gwas_model, min_maf=args.gwas_min_maf,
            min_mac=args.gwas_min_mac, min_samples=args.gwas_min_samples,
            max_missing=args.gwas_max_missing, kinship_step=args.gwas_kinship_step,
            min_kinship_snps=args.gwas_min_kinship_snps,
            permutations=args.gwas_permutations, seed=args.gwas_seed,
            chunk_size=args.chunk_size, top_n=args.gwas_top_n,
            gene_file=args.gwas_annotate_genes or args.enrich_genes or args.pbsnj_genes,
            population_covariate=args.gwas_population_covariate,
            trait_transform=args.gwas_trait_transform,
            skip_dm3_length_check=args.skip_dm3_length_check,
            match_suffix_n=args.gwas_match_suffix_n,
            pbsnj_tsv=pbsnj_win_path,pbsnj_top_tsv=pbsnj_top_path,
            pbsnj_focal=args.pbsnj_focal if args.pbsnj_focal!='all' else 'EA',
            do_plots=not args.no_plots,plot_top_hits=args.gwas_plot_top_hits)
    print(f'[DONE] Results in {args.outdir}',file=sys.stderr)



def summarize_outliers(outdir,pairs,arms,windows,quantile,min_snps,min_callable_frac):
    """Rank empirical FST window tails; does not claim selection or statistical significance."""
    for window,step in windows:
        tag=f'w{window}_s{step}'
        results=[]; regions=[]
        for a,b in pairs:
            candidates=[]
            for arm in arms:
                path=outdir/'windows'/tag/f'{a}_vs_{b}.{arm}.tsv'
                if not path.exists():continue
                with open(path,newline='') as f:
                    for r in csv.DictReader(f,delimiter='\t'):
                        try:fst=float(r['fst'])
                        except ValueError:continue
                        n_snps=int(r['n_snps']);call=int(r['n_callable_sites'])
                        span=int(r['end_0based_exclusive'])-int(r['start_0based'])
                        if (np.isfinite(fst) and n_snps>=min_snps
                                and call>=span*min_callable_frac):
                            candidates.append((arm,int(r['start_0based']),int(r['end_0based_exclusive']),fst,n_snps,call))
            if not candidates:continue
            cutoff=float(np.quantile([c[3] for c in candidates],quantile))
            tops=sorted((c for c in candidates if c[3]>=cutoff),key=lambda x:(arms.index(x[0]),x[1]))
            for arm,start,end,fst,ns,call in tops:
                results.append((a,b,arm,start,end,fst,ns,call,cutoff))
            for arm in arms:
                group=[t for t in tops if t[0]==arm]
                if not group:continue
                start,end=group[0][1:3]; vals=[group[0][3]]; count=1
                for r in group[1:]:
                    if r[1] <= end:
                        end=max(end,r[2]);vals.append(r[3]);count+=1
                    else:
                        regions.append((a,b,arm,start,end,count,max(vals),float(np.mean(vals)),cutoff))
                        start,end=r[1:3];vals=[r[3]];count=1
                regions.append((a,b,arm,start,end,count,max(vals),float(np.mean(vals)),cutoff))
        od=outdir/'outliers'/tag
        write_table(od/'outlier_windows.tsv',['pop1','pop2','arm','start_0based','end_0based_exclusive','fst','n_snps','n_callable_sites','empirical_fst_cutoff'],results)
        write_table(od/'merged_outlier_regions.tsv',['pop1','pop2','arm','start_0based','end_0based_exclusive','n_outlier_windows','max_fst','mean_fst','empirical_fst_cutoff'],regions)


def plot_results(outdir,pops,pairs,arms,windows):
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['pdf.fonttype'] = 42  # embed TrueType fonts for Illustrator editing
    matplotlib.rcParams['ps.fonttype'] = 42
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    try:
        from scipy.cluster.hierarchy import linkage, leaves_list
        from scipy.spatial.distance import squareform
    except ImportError:
        linkage=leaves_list=squareform=None
    def read_tsv(p):
        with open(p,newline='') as h:return list(csv.DictReader(h,delimiter='\t'))
    rec=read_tsv(outdir/'fst_genomewide.tsv')
    mat=np.zeros((len(pops),len(pops)));mat[:]=np.nan
    np.fill_diagonal(mat,0)
    for r in rec:
        i=pops.index(r['pop1']);j=pops.index(r['pop2'])
        mat[i,j]=mat[j,i]=float(r['hudson_fst'])
    order=np.arange(len(pops))
    if linkage is not None and np.all(np.isfinite(mat)) and len(pops)>2:
        dist=np.maximum(mat,0)
        order=leaves_list(linkage(squareform(dist,checks=False),method='average'))
    fig,ax=plt.subplots(figsize=(7,6))
    v=mat[np.ix_(order,order)]
    im=ax.imshow(v,cmap='viridis')
    ax.set_xticks(range(len(pops)),[pops[i] for i in order],rotation=45,ha='right')
    ax.set_yticks(range(len(pops)),[pops[i] for i in order])
    for i in range(len(pops)):
        for j in range(len(pops)):
            if np.isfinite(v[i,j]):ax.text(j,i,f'{v[i,j]:.3f}',ha='center',va='center',fontsize=8,
                                            color='white' if v[i,j] < np.nanmax(v)*0.55 else 'black')
    fig.colorbar(im,ax=ax,label='Genome-wide Hudson FST')
    ax.set_title('Pairwise population differentiation\nHierarchically ordered (average linkage)')
    fig.tight_layout()
    (outdir/'plots').mkdir(exist_ok=True)
    fig.savefig(outdir/'plots'/'pairwise_fst_heatmap.pdf',dpi=180);plt.close(fig)
    # choose nearest 50 kb window (if absent use largest specified)
    w,s=min(windows,key=lambda x:abs(x[0]-50000))
    for a,b in pairs:
        fig,axes=plt.subplots(len(arms),1,figsize=(12,2.1*len(arms)),sharey=True)
        axes=np.atleast_1d(axes)
        for ax,arm in zip(axes,arms):
            f=outdir/'windows'/f'w{w}_s{s}'/f'{a}_vs_{b}.{arm}.tsv'
            if not f.exists():continue
            rows=read_tsv(f)
            xx=np.array([(int(r['start_0based'])+int(r['end_0based_exclusive']))/2/1e6 for r in rows])
            yy=np.array([float(r['fst']) for r in rows])
            ax.plot(xx,yy,lw=0.8)
            ax.set_ylabel(f'{arm}\nFST')
            ax.grid(alpha=0.15)
        axes[-1].set_xlabel('Chromosome-arm position (Mb)')
        fig.suptitle(f'{a} vs {b}: Hudson FST, {w//1000} kb windows, {s//1000} kb step')
        fig.tight_layout(rect=(0,0,1,0.96))
        fig.savefig(outdir/'plots'/f'fst_scan_{a}_vs_{b}.pdf',dpi=150)
        plt.close(fig)

if __name__=='__main__':main()
