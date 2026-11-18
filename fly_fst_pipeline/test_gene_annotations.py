import csv
import gzip
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from pbsnj import annotate_windows, distinct_peak_labels
ROOT=Path(__file__).resolve().parent

class AnnotationsTests(unittest.TestCase):
    def test_annotation_gzip_and_nearest(self):
        with tempfile.TemporaryDirectory() as td:
            fn=Path(td)/'genes.gff.gz'
            with gzip.open(str(fn),'wt') as fh:
                fh.write('##gff-version 3\n')
                fh.write('2L\tx\tgene\t11\t35\t.\t+\t.\tID=FBgn0001;Name=metal1;description=metal%20transport;gene_biotype=protein_coding\n')
                fh.write('chr2L\tx\tgene\t101\t140\t.\t-\t.\tID=FBgn0002;Name=metal2;Note=ion%20response\n')
            rows=[{'arm':'2L','start_0based':10,'end_0based_exclusive':20,'PBSnj_EA':0.4},
                  {'arm':'2L','start_0based':65,'end_0based_exclusive':75,'PBSnj_EA':0.8},
                  {'arm':'2L','start_0based':130,'end_0based_exclusive':150,'PBSnj_EA':0.5}]
            fields=annotate_windows(rows,fn)
            self.assertIn('overlapping_gene_ids',fields)
            self.assertEqual(rows[0]['overlapping_gene_names'],'metal1')
            self.assertEqual(rows[0]['nearest_gene_description'],'metal transport')
            self.assertEqual(rows[1]['n_overlapping_genes'],0)
            self.assertEqual(rows[1]['nearest_gene_name'],'metal2')
            self.assertEqual(rows[2]['overlapping_gene_names'],'metal2')
            self.assertEqual(rows[2]['nearest_gene_distance_bp'],0)
            peaks=distinct_peak_labels(rows,'EA','2L',0.4,10,1)
            self.assertTrue(peaks)

    def test_integrated_run_with_annotation(self):
        with tempfile.TemporaryDirectory() as td:
            t=Path(td)
            gff=t/'genes.gff.gz'
            with gzip.open(str(gff),'wt') as fh:
                fh.write('2L\tx\tgene\t1\t40\t.\t+\t.\tID=FBgn101;Name=GeneOne\n')
                fh.write('2L\tx\tgene\t50\t90\t.\t+\t.\tID=FBgn102;Name=GeneTwo\n')
            out=t/'out'
            subprocess.run([sys.executable,str(ROOT/'fly_fst.py'),
               '--seq-dir',str(ROOT/'test_input'),'--metadata',str(ROOT/'test_input/metadata.csv'),
               '--outdir',str(out),'--arms','2L','--pbsnj','--pbsnj-focal','EA',
               '--pbsnj-populations','EA,ED,SP','--pbsnj-genes',str(gff),
               '--pbsnj-peak-labels','3','--windows','20:10','--skip-dm3-length-check',
               '--chunk-size','13','--min-window-snps','1','--no-sites'],
               check=True,stdout=subprocess.DEVNULL)
            for f in ('pbsnj.tsv','EA.top1pct.tsv'):
                with (out/'pbsnj/w20_s10'/f).open() as fh:
                    rows=list(csv.DictReader(fh,delimiter='\t'))
                self.assertIn('nearest_gene_name',rows[0] if rows else {})
            self.assertTrue((out/'pbsnj/w20_s10/EA.pbsnj.pdf').is_file())

if __name__=='__main__': unittest.main()
