#!/usr/bin/env python3
"""Integration smoke test: existing FST/PBSnj workflow + new trait-GWAS flag."""
import csv
import subprocess
import sys
import tempfile
from pathlib import Path


def test_continuous_trait_cli():
    root=Path(__file__).resolve().parent
    with tempfile.TemporaryDirectory() as t:
        work=Path(t)
        metadata=root/'test_input'/'metadata.csv'
        with metadata.open(newline='') as h:
            names=[r['sample'] for r in csv.DictReader(h)]
        pheno=work/'traits.csv'
        with pheno.open('w',newline='') as h:
            w=csv.writer(h)
            w.writerow(['sample','Size','Color'])
            for i, name in enumerate(names):
                w.writerow([name,i+1,(i*7)%9])
        cmd=[sys.executable,str(root/'fly_fst.py'),
             '--seq-dir',str(root/'test_input'),'--metadata',str(metadata),
             '--outdir',str(work/'result'),'--arms','2L',
             '--windows','20:10','--chunk-size','19','--skip-dm3-length-check',
             '--no-plots','--no-sites','--pbsnj','--pbsnj-focal','EA',
             '--continuous-trait-gwas','--trait-file',str(pheno),
             '--trait-columns','Size,Color','--gwas-model','linear',
             '--gwas-min-samples','8','--gwas-min-mac','2',
             '--gwas-kinship-step','1']
        result=subprocess.run(cmd,cwd=str(root),text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        if result.returncode:
            raise AssertionError('CLI failed:\n'+result.stderr+'\n'+result.stdout)
        assert (work/'result'/'fst_genomewide.tsv').exists()
        assert (work/'result'/'pbsnj'/'w20_s10'/'pbsnj.tsv').exists()
        assert (work/'result'/'continuous_trait_gwas'/'summary.tsv').exists()
        assert (work/'result'/'continuous_trait_gwas'/'Size'/'top_variants.tsv').exists()
        assert (work/'result'/'continuous_trait_gwas'/'Color'/'top_variants.tsv').exists()
        print('PASS: full FST + PBSnj + two independent continuous-trait GWAS analyses')


if __name__=='__main__':test_continuous_trait_cli()
