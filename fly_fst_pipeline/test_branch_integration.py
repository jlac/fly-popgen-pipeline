#!/usr/bin/env python3
"""End-to-end test for branch calibration from a miniature FST output."""
import csv, json, math, tempfile
from pathlib import Path
import numpy as np
from branch_stats import run

with tempfile.TemporaryDirectory() as td:
    root=Path(td)
    (root/'run_config.json').write_text(json.dumps({'populations':['A','B','C'],'window_sizes_and_steps':[[20,20]]}))
    base=root/'windows'/'w20_s20';base.mkdir(parents=True)
    # Three windows; all numerators are already represented as window FST.
    branch_lengths=[(.3,.5,.4),(.4,.5,.3),(.5,.7,.2)]
    for a,b,idx in [('A','B',0),('A','C',1),('B','C',2)]:
        with open(base/f'{a}_vs_{b}.2L.tsv','w',newline='') as f:
            w=csv.writer(f,delimiter='\t')
            w.writerow(['start_0based','end_0based_exclusive','fst','n_snps','n_callable_sites'])
            for k,triple in enumerate(branch_lengths):
                w.writerow([k*20,(k+1)*20,1-math.exp(-triple[idx]),10,20])
    run(root,baseline='arm',min_snps=3,min_callable_fraction=.5,do_plots=False)
    fp=root/'branch_stats'/'w20_s20'/'A__B__C.tsv'
    with open(fp) as f:records=list(csv.DictReader(f,delimiter='\t'))
    assert len(records)==3
    assert np.allclose([float(x['PBS']) for x in records],[.2,.3,.5])
    expected=(.2 - .4 * (.3/.3))
    assert np.isclose(float(records[0]['PBE']),expected)
    print('End-to-end PBS/PBE tests passed')
