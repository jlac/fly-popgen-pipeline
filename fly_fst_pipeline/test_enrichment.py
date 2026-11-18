import csv
import tempfile
import unittest
from pathlib import Path
import numpy as np
from gene_set_enrichment import gene_window_pairs, score_genes, permute_blocks, run, bh_adjust

class GeneEnrichmentTests(unittest.TestCase):
    def test_overlap_and_gene_scores(self):
        w=[('2L',0,100,1.),('2L',50,150,3.),('2L',100,200,5.)]
        g=[('2L',0,60,'g1','g1'),('2L',100,150,'g2','g2'),('2L',400,500,'g3','g3')]
        links, starts=gene_window_pairs(w,g)
        self.assertEqual(links.tolist(),[0,1,1,2])
        self.assertTrue(np.allclose(score_genes(np.array([1.,3.,5.]),links,starts,'max')[:2],[3.,5.]))
        self.assertTrue(np.allclose(score_genes(np.array([1.,3.,5.]),links,starts,'mean')[:2],[2.,4.]))
        self.assertTrue(np.isnan(score_genes(np.array([1.,3.,5.]),links,starts,'max')[2]))
    def test_block_permutation_preserves_multiset(self):
        w=[('2L',i*100,(i+1)*100,float(i)) for i in range(30)]
        x=np.arange(30.,dtype=float)
        y=permute_blocks(x,w,500,np.random.default_rng(11))
        np.testing.assert_array_equal(np.sort(x),np.sort(y))
        self.assertFalse(np.array_equal(x,y))
        # Consecutive original values within each full block remain consecutive.
        for a in range(0,30,5):self.assertTrue(np.all(np.diff(y[a:a+5])==1))
    def test_bh(self):
        self.assertTrue(np.allclose(bh_adjust([.01,.04,.03]),[.03,.04,.04]))
    def test_end_to_end(self):
        with tempfile.TemporaryDirectory() as temp:
            p=Path(temp)
            rng=np.random.default_rng(7)
            with (p/'pbsnj.tsv').open('w',newline='') as h:
                writer=csv.writer(h,delimiter='\t')
                writer.writerow(['arm','start_0based','end_0based_exclusive','passes_filter','PBSnj_EA'])
                for i in range(80):writer.writerow(['2L',i*1000,(i+1)*1000,1,round(float(rng.normal())+(2 if i<10 else 0),5)])
                for i in range(60):writer.writerow(['3R',i*1000,(i+1)*1000,1,round(float(rng.normal()),5)])
            with (p/'genes.bed').open('w') as h:
                for a,n in [('2L',80),('3R',60)]:
                    for i in range(n):h.write(f'{a}\t{i*1000}\t{(i+1)*1000}\t{a}_{i}\n')
            with (p/'sets.tsv').open('w') as h:
                for i in range(10):h.write(f'candidate\t2L_{i}\n')
                for i in range(10,30):h.write(f'background\t3R_{i-10}\n')
            rows=run(p/'pbsnj.tsv',p/'genes.bed',p/'sets.tsv',p/'out',permutations=99,
                     block_bp=5000,min_genes=5,seed=11)
            self.assertEqual(len(rows),2)
            self.assertTrue((p/'out'/'enrichment.tsv').exists())
            self.assertTrue(all(0<r['p_empirical_one_sided']<=1 for r in rows))
            self.assertTrue(rows[0]['q_BH']>=rows[0]['p_empirical_one_sided'])

if __name__=='__main__':unittest.main()
