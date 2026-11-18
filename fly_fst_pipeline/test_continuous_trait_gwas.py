#!/usr/bin/env python3
"""Self-contained synthetic end-to-end test: no real data required."""
import csv
import gzip
import json
import tempfile
from pathlib import Path

import numpy as np

from continuous_trait_gwas import run, file_sequence_length, bh_adjust


def test_continuous_trait_gwas():
    with tempfile.TemporaryDirectory() as temp:
        root=Path(temp)
        inp=root/'all_arms';inp.mkdir()
        lines=['EA1','EA2','EA3','EA4','EA5','ED1','ED2','ED3','ED4',
               'ZI1','ZI2','ZI3','CO1','CO2','CO3','FR1','FR2','FR3','FR4','FR5',
               'UG1','UG2','UG3','SP1','SP2','SP3']
        rng=np.random.RandomState(42)
        n=len(lines)
        gr=rng.binomial(1,.5,n)
        # Artificial genotype/phenotype association on 2L:777.
        genetic={}
        for arm in ('2L','2R','3L'):
            matrix=np.full((n,1900),ord('A'),dtype=np.uint8)
            for p in range(13,1900,13):
                g=rng.binomial(1,.25,n).astype(bool)
                matrix[g,p]=ord('C')
            if arm=='2L':
                matrix[gr.astype(bool),776]=ord('T')
                matrix[~gr.astype(bool),776]=ord('C')
                matrix[0,812]=ord('N')
            genetic[arm]=matrix
            for i,s in enumerate(lines):
                if arm=='3L' and s in ('SP1','CO2'):
                    continue
                with open(inp/'{}_Chr{}.seq'.format(s,arm),'wb') as f:
                    f.write(matrix[i].tobytes())
                    if arm=='2R' and s=='ZI3':
                        f.write(b'\n')  # trailing newline must not count as base
        assert file_sequence_length(inp/'ZI3_Chr2R.seq')==1900
        with (root/'metadata.csv').open('w',newline='') as f:
            w=csv.writer(f);w.writerow(['sample','population'])
            for s in lines:
                w.writerow([s, ''.join(x for x in s if x.isalpha())])
        with (root/'traits.csv').open('w',newline='') as f:
            w=csv.writer(f);w.writerow(['sample','Mn','Zn'])
            for i,s in enumerate(lines):
                # Planted gene signal, independent random second trait
                w.writerow([s,10+5*gr[i]+rng.normal(0,.6),rng.normal(2,1)])
        (root/'genes.gff').write_text('2L\tFlyBase\tgene\t690\t900\t.\t+\t.\tID=FBgn00001;Name=PlantedGene\n')
        (root/'pbsnj.tsv').write_text('arm\tstart_0based\tend_0based_exclusive\tpasses_filter\tPBSnj_EA\n2L\t700\t850\t1\t0.58\n')
        (root/'EA.top1pct.tsv').write_text('arm\tstart_0based\tend_0based_exclusive\tpasses_filter\tPBSnj_EA\n2L\t700\t850\t1\t0.58\n')
        output=root/'gwas'
        result=run(inp,root/'metadata.csv',root/'traits.csv',['Mn','Zn'],output,
                   arms=['2L','2R','3L'],min_mac=3,min_samples=15,
                   kinship_step=13,min_kinship_snps=20,permutations=9,
                   chunk_size=200,top_n=50,gene_file=root/'genes.gff',
                   skip_dm3_length_check=True,
                   pbsnj_tsv=root/'pbsnj.tsv',pbsnj_top_tsv=root/'EA.top1pct.tsv')
        assert len(result)==2
        assert result[0]['trait']=='Mn' and result[1]['trait']=='Zn'
        assert all(x['n_tests']>0 for x in result)
        with (output/'Mn'/'top_variants.tsv').open(newline='') as f:
            rows=list(csv.DictReader(f,delimiter='\t'))
        assert any(r['arm']=='2L' and r['pos_1based']=='777' for r in rows),rows[:3]
        signal=next(r for r in rows if r['arm']=='2L' and r['pos_1based']=='777')
        assert float(signal['p'])<0.01, signal
        assert signal['overlapping_gene_names']=='PlantedGene'
        assert signal['overlaps_EA_top1pct']=='1'
        assert float(signal['max_overlapping_PBSnj_EA'])==0.58
        assert (output/'Mn'/'top_PBSnj_overlap.tsv').exists()
        assert signal['q_BH_trait']!='' and signal['q_BH_all_traits']!=''
        assert signal['p_maxT_arm_FWER']!=''
        with gzip.open(output/'Mn'/'associations.2L.tsv.gz','rt') as f:
            g=list(csv.DictReader(f,delimiter='\t'))
        assert all('q_BH_trait' in r and 'p_maxT_arm_FWER' in r for r in g)
        assert (output/'Mn'/'manhattan.pdf').stat().st_size>500
        assert (output/'Zn'/'qq.pdf').stat().st_size>500
        assert (output/'Mn'/'genotype_trait_top_hits.pdf').stat().st_size>500
        assert (output/'sample_matching.tsv').exists()
        methods=json.loads((output/'methods.json').read_text())
        assert methods['n_matched_lines']==26
        assert methods['kinship_markers_by_arm']['2L']>0
        assert np.allclose(bh_adjust([0.01,0.04,0.03]),[0.03,0.04,0.04])
        print('PASS: two separate traits, LOCO LMM, 3 arms, one missing arm, tail newline, 9 permutations, FDR, annotation, PDFs')


if __name__=='__main__':test_continuous_trait_gwas()
