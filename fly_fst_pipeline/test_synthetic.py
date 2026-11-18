#!/usr/bin/env python3
"""Regression test using bundled miniature sequences."""
import csv, gzip, math, subprocess, sys, tempfile
from pathlib import Path
root=Path(__file__).resolve().parent
with tempfile.TemporaryDirectory() as t:
    out=Path(t)/'out'
    subprocess.run([sys.executable,str(root/'fly_fst.py'),
        '--seq-dir',str(root/'test_input'),'--metadata',str(root/'test_input'/'metadata.csv'),
        '--outdir',str(out),'--arms','2L','--windows','20:10,50:10',
        '--chunk-size','13','--skip-dm3-length-check','--min-window-snps','1'],check=True)
    with gzip.open(out/'snp_fst'/'EA_vs_ED.2L.tsv.gz','rt') as f:
        rec=list(csv.DictReader(f,delimiter='\t'))
    r=next(x for x in rec if x['pos_1based']=='10')
    assert r['allele0']=='A' and r['allele1']=='G'
    assert float(r['freq_allele1_pop1'])==0
    assert float(r['freq_allele1_pop2'])==1
    assert math.isclose(float(r['hudson_fst']),1.0)
    with (out/'fst_genomewide.tsv').open() as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    fst={(r['pop1'],r['pop2']):float(r['hudson_fst']) for r in rows}
    assert math.isclose(fst['EA','ED'],0.75,abs_tol=1e-9)
    assert math.isnan(fst['ED','SP'])
    with (out/'windows'/'w20_s10'/'EA_vs_ED.2L.tsv').open() as f:
        w=list(csv.DictReader(f,delimiter='\t'))
    assert math.isclose(float(w[0]['fst']),0.75,abs_tol=1e-9)
    assert int(w[2]['n_callable_sites'])==19
    assert (out/'plots'/'pairwise_fst_heatmap.pdf').exists()
    assert (out/'outliers'/'w20_s10'/'outlier_windows.tsv').exists()
print('Synthetic regression tests passed.')
