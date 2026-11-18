# Continuous-trait GWAS module for Fly FST / PBSnj

**Command:** `--continuous-trait-gwas` in `fly_fst.py`, or standalone `continuous_trait_gwas.py` to avoid repeating FST / PBSnj. Each numeric trait is analyzed **independently**, never combined into an omnibus response. A separate set of association tests, effects, P-values, BH q-values, tables, and vector PDF figures is produced for each trait.

## 1. Requirements

- `numpy`, `scipy`, `matplotlib`. Supports Python 3.7 with appropriately compatible package versions.
- dm3/BDGP Release 5 coordinate-aligned `.seq` files named `EA49_Chr2L.seq`, etc., containing one `A/C/G/T/N` character per genomic coordinate. Only **one allele per position** is available. These may be consensus sequences from diploid lines, in which case heterozygous genotypes cannot be recovered. Verify genotype calling before interpreting associations.
- Metadata CSV with **`sample,population`**, one row per sequenced line.
- Trait CSV with **one row per genetically distinct line** and numeric trait columns. Multiple technical measurements of a line must be averaged or modeled separately upstream; do NOT enter technical replicates as independent genotype samples.
- Assayed lines must match genotyped sample IDs. Use your own trait table; no real phenotypes or sequencing data are distributed with this public repository. Verify all genotype/phenotype sample aliases independently.

## 2. Recommended: standalone GWAS on already-analyzed genomes

Run this from the directory holding `all_arms/`, `fly_fst_pipeline/`, the gene GFF3, and your phenotype CSV:

```bash
python fly_fst_pipeline/continuous_trait_gwas.py \
  --seq-dir all_arms \
  --metadata all_arms/metadata.csv \
  --trait-file traits.csv \
  --trait-columns Ca,Fe,K,Mg,Mn,Na,P,S,Sr,Zn \
  --trait-sample-column sample \
  --outdir PBSnj_out/continuous_trait_gwas \
  --arms 2L,2R,3L,3R,X \
  --gwas-model lmm \
  --gwas-min-maf 0.05 \
  --gwas-min-mac 3 \
  --gwas-min-samples 15 \
  --gwas-max-missing 0.2 \
  --gwas-kinship-step 5000 \
  --gwas-annotate-genes dmel-all-no-analysis-r5.57.gff \
  --gwas-pbsnj-tsv PBSnj_out/pbsnj/w10000_s5000/pbsnj.tsv \
  --gwas-pbsnj-top-tsv PBSnj_out/pbsnj/w10000_s5000/EA.top1pct.tsv \
  --gwas-pbsnj-focal EA
```

If PBSnj tables do not exist or have no passing finite windows, omit the `--gwas-pbsnj-*` options; the association analysis runs independently. Missing `.seq` files for an arm automatically exclude **that line for that arm only**. An arm/trait is skipped and reported if fewer than `--gwas-min-samples` lines are available.

A fast exploratory linear model is available via `--gwas-model linear` (consider `--gwas-population-covariate`), but the LMM is preferable for related lines. To emphasize within-population evidence, you can add `--gwas-population-covariate`, but with 26 lines and 7 populations this is likely low-powered.

## 3. Full existing pipeline + continuous-trait GWAS

Append these switches to the existing `fly_fst.py` command (while retaining `--seq-dir`, `--metadata`, `--outdir`, `--arms`, `--allow-incomplete-arms`, `--pbsnj`, etc.):

```bash
--continuous-trait-gwas \
--trait-file traits.csv \
--trait-columns Ca,Fe,K,Mg,Mn,Na,P,S,Sr,Zn \
--trait-sample-column sample \
--gwas-model lmm \
--gwas-min-maf 0.05 \
--gwas-min-mac 3 \
--gwas-min-samples 15 \
--gwas-annotate-genes dmel-all-no-analysis-r5.57.gff \
--gwas-pbsnj-window 10000:5000
```

GWAS is invoked after population-statistics, PBSnj, and (if requested) GO enrichment. `--gwas-pbsnj-window` selects which PBSnj window-size results to annotate the top associated variants against. `--no-plots` suppresses all PDF plots. Standalone GWAS has its own `--gwas-no-plots` option.

**Note:** Even if population-level PBSnj fails because a population has too few sequences for one arm, standalone GWAS can still analyze that arm if at least 15 genotype/phenotype-matched lines exist. You do not need to rerun PBSnj to use the quantitative-trait module.

## 4. Association models and filters

- LMM: construct a kinship/relatedness matrix using widely spaced biallelic sites (`--gwas-kinship-step`, default 5000 bp). For each arm, exclude markers from that arm when calculating kinship (**LOCO**, leave-one-arm-out). Fit one REML null variance ratio per trait/arm; run generalized least-squares additive allelic effects with t tests. This is an EMMAX/P3D-style approximation, not an exact refit of the mixed model for every marker.
- Linear: ordinary regression, optionally with fixed population indicators. It does not account for polygenic relatedness unless relevant covariates are specified.
- Default minor allele frequency >= 0.05 and **minor allele count >= 3** among called lines, sample call fraction >= 80%; minimum usable line count = 15. Missing marker genotypes are allele-mean imputed for GLS, and the called count is reported.
- Each SNP is biallelic across the phenotyped lines available for the arm. Multiallelic SNPs are skipped. `allele0` / `allele1` are the two alleles in lexicographic order; **`allele1` is NOT necessarily the alternate/reference-derived allele**. `beta_allele1` is the trait change for having allele1 versus allele0 as represented by the haploid `.seq` file.
- Trait transformations: `--gwas-trait-transform none` (default), `log`, or `rank-int`; effects under transformations are in transformed units.
- Multiple testing: BH q-values among all tested variants **within each trait** (`q_BH_trait`), and **across all ten traits** (`q_BH_all_traits`). Neither is a guarantee of validity if population structure or model calibration is wrong.
- Optional `--gwas-permutations N`: costly per-arm, per-trait Freedman–Lane-style residual permutations after GLS whitening. Produces an **approximate per-arm family-wise maxT P-value**, not genome-wide across arms or all traits. Default **0** because scanning millions of SNPs with 10,000 permutations can be prohibitive; test 50–100 if desired, then scale as resources permit.
- Annotation uses FlyBase gene intervals matching the `.seq` genome assembly. Overlapping genes are reported along with the nearest gene, **not** an assertion of causal mediation.
- GWAS/PBSnj integration labels top associations with the maximum valid focal PBSnj in overlapping windows and whether they overlap a top-1%-PBSnj window. The top-1% threshold is descriptive, not a selection significance test.

## 5. Output layout

```
PBSnj_out/continuous_trait_gwas/
  sample_matching.tsv
  trait_arm_qc.tsv
  summary.tsv
  methods.json
  Mn/
    associations.2L.tsv.gz    # all tested biallelic SNPs on 2L
    associations.2R.tsv.gz    # other arms similarly
    top_variants.tsv           # top N across all tested arms; gene annotations
    candidate_genes.tsv        # genes represented among top N associations
    top_PBSnj_overlap.tsv      # top hits also in top-1% PBSnj windows (if supplied)
    manhattan.pdf              # vector plot, per 5-kb bin best SNP for editability
    qq.pdf                     # vector QQ plot based on all tested SNPs
    genotype_trait_top_hits.pdf # vector multipage raw-trait plots by allele
  Ca/ ... Fe/ ... K/ ... etc.
```

`associations.ARM.tsv.gz` contains `arm`, `pos_1based`, `allele0`, `allele1`, `n_samples`, `n_called`, `minor_allele_count`, `maf`, `beta_allele1`, `se_beta`, `t`, `p`, `delta_lmm`, `q_BH_trait`, `q_BH_all_traits` and `p_maxT_arm_FWER` (blank if permutations disabled).

The Manhattan plot uses one lowest-P SNP in every 5-kb bin (or further downsampling if necessary). This retains peak locations but **does not draw every tested variant**; the full association tables are complete. The QQ plot is computed from all P-values and uses sampled rank points for vector PDF size. Multi-page genotype–trait PDFs include at most 5 top hits at least 50 kb apart.

## 6. Notes on inference

The matched dataset has only **26 phenotyped genetic lines** distributed across seven populations. Millions of SNP tests will have very low power for a stringent genome-wide GWAS, and population-specific trait differences can be confounded with population history. The small-sample LMM approximation does not eliminate this problem. Interpret results as **exploratory candidate associations**, ideally confirmed through targeted variants, cross-validation, experimental follow-up, and sensitivity to population covariates.

Thinning kinship markers by physical distance is not full LD-pruning, and natural inversions/linked selection violate some simplifying assumptions. Mean imputation at missing calls can bias small-sample tests. An arm with few SNPs outside it will be skipped rather than quietly running an uncorrected model. For genome-wide scans, check effective line counts and model convergence in `trait_arm_qc.tsv` and `methods.json`.

For a wider and more definitive GWAS, phenotype additional sequenced lines; the 400 sequenced lines without trait measurements do not add power to association estimates.

## 7. Tests

```
cd fly_fst_pipeline
python test_continuous_trait_gwas.py
python test_continuous_trait_cli.py
```

The synthetic test checks two independently analyzed traits, 3 arms, missing chromosome-arm files, a trailing newline in a `.seq` input, LOCO mixed-model fitting, BH correction, 9 small-sample permutations, gene annotations, and vector PDF files. The second smoke test checks `fly_fst.py` integration with existing population genetics and PBSnj stages.
