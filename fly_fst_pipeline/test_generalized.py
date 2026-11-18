import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parent

class GeneralizedTests(unittest.TestCase):
    def test_prefix_and_mapping(self):
        with tempfile.TemporaryDirectory() as tmp:
            t=Path(tmp);d=t/'seq';d.mkdir()
            for sample in ['CK1','EA49','UG2','Other_A']:
                (d/f'{sample}_Chr2R.seq').write_text('ACGTN')
            mapping=t/'map.csv';mapping.write_text('sample,population\nOther_A,XX\n')
            output=t/'metadata.csv'
            subprocess.run([sys.executable,str(ROOT/'make_metadata.py'),'--seq-dir',str(d),'--output',str(output),'--sample-map',str(mapping)],check=True,stdout=subprocess.DEVNULL)
            with output.open() as f:rows=list(csv.DictReader(f))
            self.assertEqual({r['sample']:r['population'] for r in rows},{'CK1':'CK','EA49':'EA','UG2':'UG','Other_A':'XX'})
    def test_subset_runs(self):
        with tempfile.TemporaryDirectory() as tmp:
            t=Path(tmp);out=t/'out'
            subprocess.run([sys.executable,str(ROOT/'fly_fst.py'),'--seq-dir',str(ROOT/'test_input'),
              '--metadata',str(ROOT/'test_input/metadata.csv'),'--outdir',str(out),'--arms','2L',
              '--populations','EA,ED','--windows','20:10','--skip-dm3-length-check',
              '--chunk-size','13','--no-plots','--no-sites'],check=True,stdout=subprocess.DEVNULL)
            with (out/'fst_genomewide.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual(len(rows),1)
            self.assertEqual((rows[0]['pop1'],rows[0]['pop2']),('EA','ED'))
    def test_integrated_pbsnj_subset(self):
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'out'
            subprocess.run([sys.executable,str(ROOT/'fly_fst.py'),'--seq-dir',str(ROOT/'test_input'),
                '--metadata',str(ROOT/'test_input/metadata.csv'),'--outdir',str(out),'--arms','2L',
                '--pbsnj','--pbsnj-focal','EA','--pbsnj-populations','EA,ED,SP',
                '--windows','20:10','--skip-dm3-length-check','--chunk-size','13',
                '--min-window-snps','1','--no-plots','--no-sites'],check=True,stdout=subprocess.DEVNULL)
            with (out/'pbsnj/w20_s10/pbsnj.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
            self.assertEqual(rows[0]['n_populations'],'3')
            self.assertIn('lowest_callability_pair',rows[0])
if __name__=='__main__':unittest.main()
