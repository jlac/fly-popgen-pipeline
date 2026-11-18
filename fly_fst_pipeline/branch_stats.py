#!/usr/bin/env python3
"""Yi et al. PBS and Yassin et al. PBE from paired window FST tables.

Requires three populations per contrast. For a focal population A and nonfocal B,C:
 T_XY = -log(1-FST_XY)
 PBS_A = (T_AB + T_AC - T_BC)/2
 PBE_A = PBS_A - T_BC * median(PBS_A)/median(T_BC)
The PBE medians are calculated from the specified baseline (chromosome arm or genome).
"""
import argparse
import csv
import itertools
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np


def fst_to_t(fst, eps=1e-10):
    """Negative Hudson estimates map to T=0; FST >=1 is saturated at 1-eps."""
    v=np.asarray(fst,dtype=float)
    return -np.log1p(-np.clip(v,0,1-eps))


def branch_values(f_ab,f_ac,f_bc,eps=1e-10):
    ta,tb,tc=(fst_to_t(x,eps) for x in (f_ab,f_ac,f_bc))
    return (ta+tb-tc)/2,tc,ta,tb


def read_tsv(path):
    with open(path,newline='') as f: return list(csv.DictReader(f,delimiter='\t'))


def write_tsv(path,header,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with open(path,'w',newline='') as f:
        w=csv.writer(f,delimiter='\t');w.writerow(header);w.writerows(rows)


def get_pair(a,b):return tuple(sorted((a,b)))


def read_triplets(path,pops):
    if path is None:
        return [(a,b,c) for a in sorted(pops) for b,c in itertools.combinations(sorted(p for p in pops if p!=a),2)]
    with open(path,newline='') as f:
        r=csv.DictReader(f)
        if not r.fieldnames or not {'focal','nonfocal1','nonfocal2'}<=set(r.fieldnames):
            raise ValueError('Triplet CSV requires columns focal,nonfocal1,nonfocal2')
        triples=[(row['focal'].strip(),row['nonfocal1'].strip(),row['nonfocal2'].strip()) for row in r]
    out=[]; seen=set()
    for a,b,c in triples:
        if len({a,b,c})!=3 or not {a,b,c}<=set(pops):raise ValueError('Invalid triplet '+str((a,b,c)))
        b,c=sorted((b,c));key=(a,b,c)
        if key not in seen:out.append(key);seen.add(key)
    return out


def run(outdir, triplets_file=None, baseline='arm', min_snps=10, min_callable_fraction=.5, eps=1e-10, do_plots=True, populations=None):
    outdir=Path(outdir)
    cfg=json.loads((outdir/'run_config.json').read_text())
    pops=sorted(populations or cfg['populations'])
    if len(pops)<3 or len(pops)!=len(set(pops)) or not set(pops)<=set(cfg['populations']):
        raise ValueError('Select at least three distinct populations present in run_config.json')
    arms=[p.name.split('.')[-2] for p in (outdir/'windows'/f'w{cfg["window_sizes_and_steps"][0][0]}_s{cfg["window_sizes_and_steps"][0][1]}').glob('*.tsv')]
    arms=sorted(set(arms))
    triples=read_triplets(triplets_file,pops)
    meta=[]
    plot_widths=[w for w,s in cfg['window_sizes_and_steps'] if w>=50000]
    plot_width=min(plot_widths) if plot_widths else max(w for w,s in cfg['window_sizes_and_steps'])
    for win,step in cfg['window_sizes_and_steps']:
        tag=f'w{win}_s{step}'
        pair_rows={}
        for a,b in itertools.combinations(pops,2):
            pathbase=outdir/'windows'/tag
            for arm in arms:
                f=pathbase/f'{a}_vs_{b}.{arm}.tsv'
                if not f.exists():continue
                pair_rows[(a,b,arm)]={ (int(r['start_0based']),int(r['end_0based_exclusive'])):r for r in read_tsv(f) }
        for a,b,c in triples:
            records=[]
            for arm in arms:
                keys=[(get_pair(x,y)+(arm,)) for x,y in ((a,b),(a,c),(b,c))]
                if any(k not in pair_rows for k in keys): continue
                d1,d2,d3=(pair_rows[k] for k in keys)
                common=sorted(set(d1)&set(d2)&set(d3))
                for start,end in common:
                    r1,r2,r3=d1[start,end],d2[start,end],d3[start,end]
                    pair=[r1,r2,r3]
                    ns=[int(r['n_snps']) for r in pair]
                    calls=[int(r['n_callable_sites']) for r in pair]
                    # A nonfocal pair may be invariant while the focal population differs.
                    # Interpret a fully callable invariant pair as zero differentiation.
                    vals=np.array([0.0 if n==0 and call>0 else float(r['fst'])
                                   for r,n,call in zip(pair,ns,calls)],dtype=float)
                    span=end-start
                    valid=(np.all(np.isfinite(vals)) and max(ns)>=min_snps
                           and all(call>=span*min_callable_fraction for call in calls))
                    if valid:
                        pbs,tbc,tab,tac=branch_values(*vals,eps)
                        pbs=float(pbs);tbc=float(tbc);tab=float(tab);tac=float(tac)
                    else:pbs=tbc=tab=tac=float('nan')
                    records.append(dict(arm=arm,start=start,end=end,valid=valid,pbs=pbs,tbc=tbc,tab=tab,tac=tac,
                                        fst_ab=vals[0],fst_ac=vals[1],fst_bc=vals[2],
                                        snps_min=min(ns),callable_min=min(calls),saturated=sum(v>=1 for v in vals)))
            baseline_groups=defaultdict(list)
            for r in records:
                if r['valid']:
                    baseline_groups[r['arm'] if baseline=='arm' else 'genome'].append(r)
            calibration={}
            for group, rows in baseline_groups.items():
                pbsmed=float(np.median([r['pbs'] for r in rows])); tmed=float(np.median([r['tbc'] for r in rows]))
                calibration[group]=(pbsmed,tmed, pbsmed/tmed if tmed>0 else float('nan'),len(rows))
            out=[]
            for r in records:
                group=r['arm'] if baseline=='arm' else 'genome'
                pm,tm,ratio,n=calibration.get(group,(float('nan'),)*3+(0,))
                expected=r['tbc']*ratio if np.isfinite(ratio) and r['valid'] else float('nan')
                pbe=r['pbs']-expected if np.isfinite(expected) else float('nan')
                out.append((r['arm'],r['start'],r['end'],a,b,c,int(r['valid']),r['snps_min'],r['callable_min'],
                            r['fst_ab'],r['fst_ac'],r['fst_bc'],r['tab'],r['tac'],r['tbc'],r['pbs'],expected,pbe,pm,tm,n,r['saturated']))
            name=f'{a}__{b}__{c}'
            dest=outdir/'branch_stats'/tag/f'{name}.tsv'
            header=['arm','start_0based','end_0based_exclusive','focal','nonfocal1','nonfocal2',
                    'passes_filter','min_pair_snps','min_pair_callable_sites','fst_focal_nonfocal1','fst_focal_nonfocal2',
                    'fst_nonfocal1_nonfocal2','T_focal_nonfocal1','T_focal_nonfocal2','T_nonfocal1_nonfocal2',
                    'PBS','PBS_expected','PBE','baseline_median_PBS','baseline_median_T_nonfocal',
                    'baseline_n_windows','n_pair_fst_at_or_above_1']
            write_tsv(dest,header,out)
            finite=[r for r in out if np.isfinite(r[17])]
            cutoff=float(np.quantile([r[17] for r in finite],.99)) if finite else float('nan')
            write_tsv(outdir/'branch_stats'/tag/f'{name}.top1pct_PBE.tsv', header,
                      [r for r in out if np.isfinite(r[17]) and r[17]>=cutoff])
            for group,(pm,tm,ratio,n) in calibration.items():
                meta.append((tag,a,b,c,group,baseline,pm,tm,ratio,n,cutoff))
            if do_plots and win==plot_width:
                plot_triplet(outdir/'branch_stats'/tag/f'{name}.pdf',out,a,b,c,arms)
    write_tsv(outdir/'branch_stats'/'calibration.tsv',
              ['window_tag','focal','nonfocal1','nonfocal2','baseline_group','baseline_mode',
               'median_PBS','median_T_nonfocal','PBS_per_T_ratio','n_valid_windows','PBE_top_1pct_cutoff'],meta)
    (outdir/'branch_stats'/'methods.json').write_text(json.dumps({
        'PBS':'(T_AB+T_AC-T_BC)/2; T=-log(1-FST)',
        'PBE':'PBS - T_BC * median(PBS)/median(T_BC)',
        'baseline':baseline,'min_snps_in_any_pair_per_window':min_snps,
        'min_callable_fraction_per_pair_per_window':min_callable_fraction,
        'FST_transform':'clip negative Hudson FST to zero; FST >=1 to 1-epsilon, only for transformed branch stats',
        'epsilon':eps,'triplet_count':len(triples),
        'notes':'PBS and PBE are descriptive. All three pairs must pass callability filters, and at least one must meet the SNP count; pair FST=0 is used for invariant fully-callable pairs. PBE undefined if median T_BC=0. Medians computed across qualifying windows, including overlapping windows if requested.'},indent=2))
    print(f'[DONE] PBS/PBE for {len(triples)} focal triplets: {outdir / "branch_stats"}')


def plot_triplet(path,rows,a,b,c,arms):
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams['pdf.fonttype'] = 42  # embed TrueType fonts for Illustrator editing
    matplotlib.rcParams['ps.fonttype'] = 42
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(len(arms),1,figsize=(12,2.2*len(arms)))
    axes=np.atleast_1d(axes)
    for ax,arm in zip(axes,arms):
        r=[v for v in rows if v[0]==arm]
        xs=[(v[1]+v[2])/2/1e6 for v in r]
        ys=[v[17] for v in r]
        ax.plot(xs,ys,lw=.8)
        ax.set_ylabel(f'{arm}\nPBE')
        ax.grid(alpha=.2)
    axes[-1].set_xlabel('Position on chromosome arm (Mb)')
    fig.suptitle(f'PBE — focal {a}, comparators {b} / {c}')
    fig.tight_layout(rect=(0,0,1,.96));path.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(path,dpi=160);plt.close(fig)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--outdir',type=Path,required=True,help='Existing fly_fst.py results directory')
    p.add_argument('--triplets',type=Path,help='Optional CSV focal,nonfocal1,nonfocal2; default all triplets')
    p.add_argument('--populations',help='Comma-separated subset of populations in existing FST run')
    p.add_argument('--baseline',choices=['arm','genome'],default='arm')
    p.add_argument('--min-window-snps',type=int,default=10)
    p.add_argument('--min-window-callable-fraction',type=float,default=.5)
    p.add_argument('--no-plots',action='store_true')
    a=p.parse_args()
    run(a.outdir,a.triplets,a.baseline,a.min_window_snps,a.min_window_callable_fraction,do_plots=not a.no_plots,
        populations=[v.strip() for v in a.populations.split(',') if v.strip()] if a.populations else None)


if __name__=='__main__': main()
