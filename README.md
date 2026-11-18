# Fly PopGen Pipeline

A chunked, Python-based workflow for exploratory population genomics of *Drosophila melanogaster* using **reference-coordinate-aligned, haploid `.seq` files**. It includes pairwise Hudson FST, sliding-window branch statistics (including PBSnj), GO gene-set enrichment, and per-trait quantitative-phenotype association scans.

**Important:** This is research software, not a validated production GWAS package. Check reference coordinates, per-arm sample availability, population structure, variant encoding, and test assumptions before biological interpretation. PBSnj outliers are not proof of selection; trait associations are not proof of causality.

## Features

- Chunked processing of chromosome-arm `.seq` files, including lines with missing chromosome arms.
- Population differentiation, windowed PBS/PBE/PBSnj, gene-annotated candidate regions.
- GO gene-set enrichment with block permutations and separate BP/MF/CC summaries.
- Editable vector PDF figures and configurable automatic/manual PBSnj gene labels.
- `--continuous-trait-gwas`: **independent** association testing for each numerical phenotype, with a relatedness-aware LMM option.
- A `tools/flybase_gaf_to_gmt.py` converter for FlyBase GO annotation to GMT.

## Requirements

- Python and the packages in `fly_fst_pipeline/requirements.txt` (Python 3.7 compatible constraints are included; Python 3.11 recommended for new environments).
- Input sequence files of consistent length and **one nucleotide call per reference coordinate** on each chromosome arm. These inputs are *not* raw reads and do not retain diploid heterozygosity if generated as haploid consensus sequences.
- Population metadata (`sample,population`), plus an optional numeric trait CSV (`sample` + trait columns).
- Optional matching-assembly GFF3 and GO GMT/OBO files for annotation and enrichment.

```bash
python -m pip install -r fly_fst_pipeline/requirements.txt
python -m pytest -q fly_fst_pipeline
```

## Example: population differentiation and PBSnj

```bash
python fly_fst_pipeline/make_metadata.py --seq-dir all_arms --output all_arms/metadata.csv

python fly_fst_pipeline/fly_fst.py \
  --seq-dir all_arms \
  --metadata all_arms/metadata.csv \
  --outdir PBSnj_out \
  --arms 2L,2R,3L,3R,X \
  --allow-incomplete-arms \
  --populations EA,ED,SP,ZI,UG,FR,CO \
  --windows 10000:5000 \
  --pbsnj --pbsnj-focal EA
```

See [the full pipeline guide](fly_fst_pipeline/README.md) and the [custom gene labeling guide](fly_fst_pipeline/CUSTOM_GENE_LABELS_README.md).

## Example: continuous-trait association scans

Provide a `traits.csv` with one row per independently genotyped line (not one row per technical measurement):

```csv
sample,Ca,Fe,Mn,Zn
EA1,1.3,2.1,0.6,1.8
EA2,1.1,1.8,0.9,1.6
ED1,1.5,2.2,0.7,1.9
```

```bash
python fly_fst_pipeline/continuous_trait_gwas.py \
  --seq-dir all_arms \
  --metadata all_arms/metadata.csv \
  --trait-file traits.csv \
  --trait-columns Ca,Fe,Mn,Zn \
  --trait-sample-column sample \
  --outdir GWAS_out \
  --arms 2L,2R,3L,3R,X \
  --gwas-model lmm
```

Each trait is tested separately. See [GWAS options and caveats](fly_fst_pipeline/CONTINUOUS_TRAIT_GWAS_README.md).

## Example: GO enrichment and plotting

```bash
python fly_fst_pipeline/gene_set_enrichment.py --help
python fly_fst_pipeline/pbsnj.py --help
python tools/flybase_gaf_to_gmt.py --help
```

Read [GO ontology output guidance](fly_fst_pipeline/GO_ONTOLOGY_README.md) and [GMT conversion guidance](docs/GO_TO_GMT.md).

## Data policy

**No real genotype/phenotype data are included.** Only tiny synthetic testing examples are shipped under `fly_fst_pipeline/test_input/`. The `.gitignore` excludes common genomic input formats and generated outputs. Do not upload identifiable or unpublished research datasets to this public repository.

## License

No software license has been selected. Public GitHub visibility does not, by itself, grant reuse or redistribution rights. The owner can add an open-source license later.
