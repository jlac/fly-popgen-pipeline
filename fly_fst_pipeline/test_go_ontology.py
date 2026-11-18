import csv
import tempfile
import unittest
from pathlib import Path
from gene_set_enrichment import read_go_namespaces, run, bh_adjust

class OntologyTests(unittest.TestCase):
    def test_ontology_separated(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d)
            (p/'go.obo').write_text('''format-version: 1.2
[Term]
id: GO:1111111
name: test process
namespace: biological_process
alt_id: GO:1111199
[Term]
id: GO:2222222
name: test function
namespace: molecular_function
[Term]
id: GO:3333333
name: test component
namespace: cellular_component
[Term]
id: GO:4444444
name: obsolete
namespace: biological_process
is_obsolete: true
''')
            ns = read_go_namespaces(p/'go.obo')
            self.assertEqual(ns['GO:1111199'], 'biological_process')
            self.assertNotIn('GO:4444444', ns)
            with (p/'pbsnj.tsv').open('w') as h:
                h.write('arm\tstart_0based\tend_0based_exclusive\tpasses_filter\tPBSnj_EA\n')
                for i in range(24):
                    h.write('2L\t{}\t{}\t1\t{}\n'.format(i*1000, (i+1)*1000, i/100.))
            with (p/'genes.bed').open('w') as h:
                for i in range(24):h.write('2L\t{}\t{}\tg{}\n'.format(i*1000,(i+1)*1000,i))
            with (p/'sets.gmt').open('w') as h:
                for term,offset in [('GO:1111111',0),('GO:2222222',8),('GO:3333333',16)]:
                    h.write(term+'\ttest\t'+'\t'.join('g{}'.format(i) for i in range(offset,offset+8))+'\n')
            result=run(p/'pbsnj.tsv',p/'genes.bed',p/'sets.gmt',p/'out',
                       permutations=19,block_bp=4000,min_genes=3,go_obo=p/'go.obo')
            self.assertEqual(len(result),3)
            for folder,expected in [('biological_process','GO:1111111'),('molecular_function','GO:2222222'),('cellular_component','GO:3333333')]:
                t=p/'out'/folder
                self.assertTrue((t/'enrichment.pdf').is_file())
                with (t/'enrichment.tsv').open() as h:rows=list(csv.DictReader(h,delimiter='\t'))
                self.assertEqual(len(rows),1)
                self.assertEqual(rows[0]['gene_set'],expected)
                self.assertAlmostEqual(float(rows[0]['q_BH']),float(rows[0]['p_empirical_one_sided']))
                self.assertIn('q_BH_global',rows[0])
            self.assertTrue((p/'out'/'enrichment.pdf').is_file())
            self.assertTrue((p/'out'/'gene_scores.tsv').is_file())

if __name__ == '__main__':unittest.main()
