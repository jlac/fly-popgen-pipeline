# FlyBase GAF → GMT converter

`flybase_gaf_to_gmt.py` converts FlyBase `gene_association.fb.gz` (GAF 2.x) into a GMT file, matching FlyBase `FBgn` gene identifiers used by the Drosophila PBSnj enrichment pipeline.

## Download the required files

```bash
wget -O gene_association.fb.gz \
  https://s3ftp.flybase.org/releases/current/precomputed_files/go/gene_association.fb.gz
wget -O go-basic.obo \
  https://current.geneontology.org/ontology/go-basic.obo
```

## Run

```bash
python flybase_gaf_to_gmt.py \
  --gaf gene_association.fb.gz \
  --obo go-basic.obo \
  --out-gmt flybase_all_go.gmt \
  --ion-gmt flybase_ion_candidate_sets.gmt \
  --min-genes 5 \
  --report flybase_go_gmt_report.json
```

The script has no third-party Python dependencies. The full GMT contains one gene set per GO term (including descendant annotations). The ion GMT unions genes from GO terms whose **names** match element-specific keywords. It is an **exploratory candidate set, not a curated or statistically validated ion pathway database**; inspect `flybase_go_gmt_report.json` for every contributing GO term before analyzing these sets.

### Recommended evidence sensitivity analysis

```bash
python flybase_gaf_to_gmt.py \
  --gaf gene_association.fb.gz --obo go-basic.obo \
  --out-gmt flybase_go_noIEA.gmt \
  --ion-gmt flybase_ion_noIEA.gmt \
  --exclude-iea --min-genes 5
```

By default, the script includes IEA, excludes annotations with `NOT` qualifier, restricts GAF database column to `FB`, and propagates annotations via `is_a` and `part_of` edges. It deliberately does not propagate regulatory edges (`regulates`). The GAF is not tied to a genome coordinate system. The **separate gene annotation GFF3/BED** must match your genomic data's coordinate system (probably dm3/BDGP5); the gene IDs in the GMT must match the GFF3 gene feature `ID` attribute or the BED fourth column.

## Verify that FlyBase gene IDs match your gene annotation

For an uncompressed dm3 GFF3 with `ID=FBgn...`, run:

```bash
python - <<'PY'
import re
from pathlib import Path

gff='dmel-all-no-analysis-r5.57.gff'
gmt='flybase_ion_candidate_sets.gmt'

genes=set()
for line in open(gff):
    if line.startswith('#'): continue
    cols=line.rstrip('\n').split('\t')
    if len(cols)<9 or cols[2]!='gene': continue
    m=re.search(r'(?:^|;)ID=([^;]+)',cols[8])
    if m: genes.add(m.group(1))

members=set()
for line in open(gmt):
    cols=line.rstrip('\n').split('\t')
    members.update(cols[2:])

print('GFF gene IDs:',len(genes))
print('GMT unique genes:',len(members))
print('GMT genes matching GFF IDs:',len(genes & members))
print('Example GFF IDs:', sorted(genes)[:5])
print('Example GMT IDs:', sorted(members)[:5])
PY
```

If overlap is near zero, inspect GFF3 identifiers: they may require an `ID` normalization or mapping step. Do not use symbol-based IDs in one file and FBgn identifiers in the other.

## Run PBSnj enrichment

```bash
python fly_fst_pipeline/gene_set_enrichment.py \
  --pbsnj-tsv chr2L_PBSnj_out/pbsnj/w10000_s5000/pbsnj.tsv \
  --genes dmel-all-no-analysis-r5.57.gff \
  --gene-format gff3 \
  --gene-sets flybase_ion_candidate_sets.gmt \
  --outdir chr2L_PBSnj_out/enrichment \
  --focal EA --permutations 10000 --block-bp 250000
```

For publication, define independent, preregistered gene sets, check GO term overlap/redundancy, and evaluate recombination/inversion and genome coverage biases.
