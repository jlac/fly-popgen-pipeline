import unittest
import tempfile
from pathlib import Path
from pbsnj import parse_gene_label_requests, custom_gene_labels, select_plot_labels, plot_scan

class CustomLabelsTests(unittest.TestCase):
    def setUp(self):
        self.genes=[('2L',0,30,'FBgn0001','heph'),('2L',40,70,'FBgn0002','Drak'),('2L',80,100,'FBgn0003','Ir11a')]
        self.rows=[{'arm':'2L','start_0based':a,'end_0based_exclusive':b,
                    'PBSnj_EA':score,'nearest_gene_name':name}
                   for a,b,score,name in [(0,20,.1,'heph'),(20,40,.9,'heph'),(40,60,.3,'Drak'),(60,80,.2,'Drak'),(80,100,.01,'Ir11a')]]
    def test_manual_low_gene_and_id(self):
        labels,errors=custom_gene_labels(self.rows,'EA','2L',self.genes,['Ir11a','FBgn0002'])
        self.assertEqual({p['label'] for p in labels},{'Ir11a','Drak'})
        self.assertFalse(errors)
    def test_manual_only(self):
        labels=select_plot_labels(self.rows,'EA','2L',.8,5,100000,self.genes,['Ir11a'],True)
        self.assertEqual([p['label'] for p in labels],['Ir11a'])
    def test_auto_threshold(self):
        labels=select_plot_labels(self.rows,'EA','2L',.8,5,1,self.genes,[],False)
        self.assertEqual(len(labels),1)
    def test_file(self):
        with tempfile.TemporaryDirectory() as td:
            fn=Path(td)/'names.txt';fn.write_text('heph\n# comment\nDrak\n')
            self.assertEqual(parse_gene_label_requests('Ir11a',fn),['Ir11a','heph','Drak'])
    def test_pdf(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'test.pdf'
            plot_scan(path,self.rows,'EA',['2L'],.8,5,1,self.genes,['Drak'],True)
            self.assertEqual(path.read_bytes()[:4],b'%PDF')

if __name__=='__main__':unittest.main()
