import unittest
import tempfile
import csv
import json
import math
from pathlib import Path
import numpy as np
from pbsnj import nj_terminal_lengths, transformed_fst, scan

class PBSNJTests(unittest.TestCase):
    def test_three_taxa_equals_pbs(self):
        d=np.array([[0,.4,.6],[.4,0,.5],[.6,.5,0]])
        t=nj_terminal_lengths(d,['EA','ED','ZI'])
        self.assertAlmostEqual(t['EA'],(.4+.6-.5)/2)
        self.assertAlmostEqual(t['ED'],(.4+.5-.6)/2)
        self.assertAlmostEqual(t['ZI'],(.6+.5-.4)/2)
    def test_seven_taxa_additive_tree(self):
        # All seven tips attach directly to a common unresolved center;
        # NJ returns exact pendant branch lengths even when ties occur.
        names=['EA','ED','ZI','SP','Ug','Fr','CO']
        lengths=np.array([.51,.08,.11,.07,.2,.1,.09])
        d=lengths[:,None]+lengths[None,:]
        np.fill_diagonal(d,0)
        t=nj_terminal_lengths(d,names)
        for p,length in zip(names,lengths):self.assertAlmostEqual(t[p],length,places=9)
    def test_relabeling_and_order(self):
        names=['EA','ED','ZI','SP']
        lengths=np.array([.41,.2,.1,.3]);d=lengths[:,None]+lengths[None,:]
        np.fill_diagonal(d,0)
        t=nj_terminal_lengths(d,names)
        order=[2,0,3,1]
        u=nj_terminal_lengths(d[np.ix_(order,order)],[names[i] for i in order])
        self.assertAlmostEqual(t['EA'],u['EA'])
    def test_transform(self):
        self.assertAlmostEqual(transformed_fst(.5),math.log(2))
        self.assertEqual(transformed_fst(-.1),0)
        self.assertTrue(np.isfinite(transformed_fst(1)))
    def test_existing_synthetic_fst_run(self):
        src=Path(__file__).parent/'test_run'
        if not (src/'run_config.json').exists():self.skipTest('No existing synthetic run')
        with tempfile.TemporaryDirectory() as td:
            # copy only the existing FST windows and config
            import shutil
            tmp=Path(td)
            shutil.copy(src/'run_config.json',tmp/'run_config.json')
            shutil.copytree(src/'windows',tmp/'windows')
            scan(tmp,focal='EA',min_snps=1,do_plots=False)
            with (tmp/'pbsnj'/'w20_s10'/'pbsnj.tsv').open() as fh:
                rows=list(csv.DictReader(fh,delimiter='\t'))
            self.assertTrue(rows)
            self.assertTrue(any(r['passes_filter']=='1' for r in rows))
            self.assertTrue((tmp/'pbsnj'/'methods.json').is_file())

if __name__=='__main__':unittest.main()
