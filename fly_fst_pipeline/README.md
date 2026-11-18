# Drosophila population genetic FST pipeline

A runnable, chunked, low-memory analysis for **aligned, haploid, unwrapped sequence files** named `EA49_Chr2L.seq`, `CO2_Chr3R.seq`, etc. This is not a read mapper: the sequence must have one character **per reference coordinate** and have the correct chromosome-arm length. The examples match *D. melanogaster* BDGP Release 5 (UCSC **dm3**), but **length matching does not prove reference identity**. Independently verify this assembly and coordinate system before interpreting loci.

## Install

```bash
python -m venv flyfst-env
source flyfst-env/bin/activate
pip install numpy scipy matplotlib
```

Python 3.9+ recommended. `numpy` is required; `scipy` and `matplotlib` are needed for clustering and plots. No full-genome VCF creation is required.

## Create and verify metadata

Place the `.seq` files in one directory. The sample identifier before `_Chr` must uniquely identify a *line* (one independent haploid genome).

```bash
python make_metadata.py --seq-dir /data/flyseq --output metadata.csv
```

**Inspect and correct** the `sample,population` assignments. The helper guesses populations solely from name prefixes (`EA`, `ED`, `SP`, `ZI`, `Ug`, `Fr`, `CO`) and cannot verify biological metadata. Example:

```csv
sample,population
EA49,EA
CO2,CO
```

All chromosome files for each sample must share the same chromosome coordinates. Default analysis requires all listed samples for each arm. If a sample is missing an entire arm, use `--allow-incomplete-arms`, which drops only that missing sample on the affected arm. Such missing samples are recorded in the QC/metadata counts; inspect the effects carefully.

## Run

```bash
python fly_fst.py \
  --seq-dir /data/flyseq \
  --metadata metadata.csv \
  --outdir fly_fst_output \
  --arms 2L,2R,3L,3R,X \
  --windows 5000:2500,10000:5000,50000:10000,100000:25000 \
  --chunk-size 100000 \
  --min-call-fraction 0.8
```

Each filename must use the exact case pattern `<sample>_Chr<arm>.seq`. Inputs should contain exactly `A,C,G,T,N`, without FASTA headers or whitespace, one base per coordinate. `N` is missing; no genotype is imputed. Lowercase ACGTN is permitted. This script is designed for **haploid, inbred line consensus sequences**, not diploid heterozygous-genotype calls.

If storage is constrained, add `--no-sites` to omit the large SNP-level compressed TSVs; window and diversity outputs are still calculated. Run `--no-plots` to omit plotting dependencies or defer plotting. Use `--help` for additional options.

**Do not use `--skip-dm3-length-check` with actual study data** unless you have independently validated a different sequence assembly; it exists for the bundled tiny synthetic tests. The five expected dm3 arm lengths are: 2L=23,011,544; 2R=21,146,708; 3L=24,543,557; 3R=27,905,053; X=22,422,827.

## Outputs

- `qc_by_sample.tsv`: per-line, per-arm call rate and missing bases.
- `fst_genomewide.tsv`: 21 possible pairwise population-wide Hudson FST estimates (ratio of sums).
- `fst_by_arm.tsv`: pairwise FST by chromosome arm.
- `snp_fst/<pop1>_vs_<pop2>.<arm>.tsv.gz`: **site-wise Hudson FST** for pairwise biallelic polymorphic sites; numerator, denominator, called haploid sample sizes and allele frequency. These files can be large. Multi-allelic sites are excluded from FST, with QC counts in window output.
- `windows/w<size>_s<step>/*.tsv`: site numerator/denominator sums, pairwise callable-site counts, SNP counts, and Hudson FST in 0-based half-open windows.
- `diversity/w<size>_s<step>/*.tsv`: per-population nucleotide diversity (mean pairwise diversity per callable base) and **complete-case Tajima's D** for sites where every sampled line has a valid base. Tajima's D is `NaN` if fewer than three lines or fewer than two segregating complete sites.
- `outliers/w<size>_s<step>/`: upper 1% FST windows after filtering (defaults: at least 10 polymorphic SNPs and 50% pairwise callable bases), and merged overlapping windows. This is **empirical ranking, not a test of selection or calibrated statistical significance**.
- `plots/pairwise_fst_heatmap.pdf`: pairwise population differentiation matrix, ordered by average-linkage clustering of nonnegative genome-wide estimates (negative estimates preserved in tables).
- `plots/fst_scan_*.pdf`: per-pair arm-wise 50-kb window FST scans (or closest requested window width).

### Statistical methods

For biallelic site s and two population samples of independent **haploid** lines, with reference-coded allele frequencies `p1, p2` and callable counts `n1,n2`, the unbiased Hudson-style numerator and denominator are:

```
N_s = (p1-p2)^2 - p1*(1-p1)/(n1-1) - p2*(1-p2)/(n2-1)
D_s = p1*(1-p2) + p2*(1-p1)
FST_site = N_s/D_s  (only where D_s>0)
FST_window = sum_s N_s / sum_s D_s
```

Unbiased nucleotide diversity within each population is `n/(n-1) * [1 - sum_a p_a^2]`, averaged over sites with the required population call rate. Includes A,C,G,T multiallelic sites. Classical Tajima's D is computed **only over fully called sites in each population** (complete-case), to maintain a constant haploid sample size. These two diversity estimates therefore have *different callable-site denominators*. Do not compare Tajima's D with the missingness-adjusted pi as if derived from identical masks.

Missing bases can produce nonrandom spatial patterns: use `n_callable_sites` to filter low-confidence windows. Site-level Hudson estimates can be negative; do not clip their numerators before taking window ratios. Windows with no informative SNPs have NaN FST. **Windows should not be treated as independent hypothesis tests**, especially with overlapping windows and chromosome inversions.

### Caveats

- The files contain no read-depth, genotype-quality, or heterozygote data. A/C/G/T calls are assumed correct; validate the source sequence generation process.
- All per-sample genomes must be the **same reference-coordinate alignment**, not just equal-length strings. Length check is only an initial sanity check. Also confirm orientation and absence of indel-induced coordinate shifts.
- Population assignment is determined *only* from the metadata CSV (not inferred in the statistical script).
- `min-called=2` enforces at least two nonmissing haplotypes in each population; with `min-call-fraction=0.8`, groups of three lines need all three called at a site.
- Chromosome X is analyzed as a haploid allele per line, consistent with the supplied encoding. Compare X to autosomes with due regard to effective population size, sex and demography.
- Strong FST peaks may reflect demography, low recombination, inversions, or ascertainment rather than positive selection. Follow up with annotation/LD, inversions, and appropriate null models.
- **No outgroup or ancestral polarization** is used, and `pi`/Tajima's D are not calibrated for inaccessible/masked genomic bases beyond ACGTN representation.
- Windows are half-open 0-based in tables; SNP coordinates are 1-based.

## Test on bundled synthetic data

```bash
python fly_fst.py \
  --seq-dir test_input --metadata test_input/metadata.csv \
  --outdir test_output --arms 2L --windows 20:10,50:10 \
  --chunk-size 13 --skip-dm3-length-check --min-window-snps 1
```

In the synthetic data, SNP chr2L:10 is fixed A in population EA and G in ED/SP, so its EA-vs-ED SNP FST must be 1.0. Chromosome-arm genome-wide FST over two variable sites is 0.75 for EA-vs-ED, and undefined for ED-vs-SP (no informative sites within those two groups).

## Example Slurm submission

Edit paths in `run_slurm.sh`, then submit with `sbatch run_slurm.sh`. The 32-GB memory setting is a cautious initial request, not a benchmark. I/O and output volume may dominate processing time. The pipeline runs sequentially and uses one CPU for core computations; the Slurm CPU request is a configurable allocation, not a promise of parallel scaling.

## PBS (Yi et al. 2010) and PBE (Yassin et al. 2016)

**New:** `fly_fst.py` can calculate the Population Branch Statistic and Population Branch Excess for **each focal population and a selected pair of nonfocal populations**. These are three-population statistics; choosing biologically appropriate comparators matters. By default all 105 distinct focal orientations of the 35 three-population sets are calculated (7 x choose(6,2)). It is often better to supply an explicit subset of triplets.

Run the full analysis and PBS/PBE together:

```bash
python fly_fst.py \
  --seq-dir /data/flyseq --metadata metadata.csv --outdir fst_results \
  --pbs-pbe --branch-triplets triplets.csv --branch-baseline arm
```

Or **re-use an existing completed FST run without re-reading any `.seq` files**:

```bash
python branch_stats.py --outdir fst_results \
  --triplets triplets.csv --baseline arm \
  --min-window-snps 10 --min-window-callable-fraction 0.5
```

Example `triplets.csv`:

```csv
focal,nonfocal1,nonfocal2
EA,ED,ZI
SP,CO,ZI
```

### Methods

For a focal population A and two nonfocal populations B and C, transform each window's ratio-of-sums Hudson FST:

```
T_AB = -log(1 - FST_AB)
T_AC = -log(1 - FST_AC)
T_BC = -log(1 - FST_BC)
PBS_A = (T_AB + T_AC - T_BC) / 2
PBE_A = PBS_A - T_BC * (median(PBS_A) / median(T_BC))
```

The **median of PBS** and **median of T_BC** are taken over all passing windows in the chosen calibration domain, **separately for each focal triplet and window configuration**. `--baseline arm` uses medians from that chromosome arm (default, similar to Yassin et al.); `--baseline genome` uses all analyzed arms. For overlapping windows, both medians are medians of the resulting overlapping windows. The PBE ratio is *not* the median of per-window PBS/T_BC ratios.

Window FST estimates are computed first; PBS/PBE operate on these window-based estimates, **not on the mean of per-SNP PBS/PBE**. Negative Hudson FST values are set to 0 **only when transforming to T**, without altering original FST outputs; estimates at or above 1 are capped at 1-1e-10 to avoid infinite branch lengths. A pair with no segregating sites and sufficient callable bases is treated as FST=0 for branch statistics. If the median T_BC is zero, PBE is undefined (NaN) rather than dividing by zero. Negative PBS is retained; no arbitrary clipping is applied to PBS or PBE.

All three pairs must have at least the configured fraction of callable genomic positions, and **at least one pair** must have the configured minimum number of informative SNPs. This allows the two nonfocal populations to be identical without discarding a focal-specific change. Sample sizes and masking vary, so consider sensitivity analyses with different thresholds.

### Additional files

- `branch_stats/w<size>_s<step>/<focal>__<nonfocal1>__<nonfocal2>.tsv` — per-window PBS, PBE, expected PBS, pairwise FST, transformed distances, validity indicator and calibration counts.
- `branch_stats/w<size>_s<step>/*.top1pct_PBE.tsv` — descriptive top-1% PBE windows for each triplet, if PBE is definable. **Not calibrated p-values**.
- `branch_stats/w<size>_s<step>/*.pdf` — per-arm PBE profiles for one selected window width (closest requested width >=50 kb; otherwise largest).
- `branch_stats/calibration.tsv` — observed baseline medians, ratio, number of eligible windows and descriptive 1% cutoff.
- `branch_stats/methods.json` — options and numerical conventions.

**Limitations:** All triplets are *not* interchangeable. Population B should typically be a biologically appropriate close comparator and C a further nonfocal reference; verify inferred relationships using genome-wide FST / population history first. Population structure, inversions, linked selection, poor recombination, and missing data can generate outliers independently of adaptive evolution. Sliding windows, especially overlapping windows, are highly dependent. The supplied `.seq` inputs are haploid consensus sequences without genotype quality; the same alignment and missing-data cautions noted above apply.

Papers:

- Yi et al. (2010), *Sequencing of 50 Human Exomes Reveals Adaptation to High Altitude*, Science, DOI: 10.1126/science.1190371.
- Yassin et al. (2016), *Recurrent specialization on a toxic fruit in an island Drosophila population*, PNAS, DOI: 10.1073/pnas.1522559113.

## NEW: multipopulation PBSnj (EA-focused)

`pbsnj.py` implements an **unrooted neighbor-joining terminal branch** statistic for a focal population using all available populations, based on Schmidt et al. (2019), *The impact of genetic adaptation on chimpanzee subspecies differentiation*, PLoS Genetics, DOI: 10.1371/journal.pgen.1008485. **This implementation works on sliding-window ratio-of-sums Hudson FST estimates, whereas the original paper calculated PBSnj on SNPs. It is a window-based adaptation, not an exact replication of the published SNP-level pipeline.**

Run from scratch and calculate PBSnj after the normal FST analysis:

```bash
python fly_fst.py \
  --seq-dir /data/flyseq --metadata metadata.csv \
  --outdir fly_fst_output --arms 2L,2R,3L,3R,X \
  --windows 50000:10000,100000:25000 \
  --pbsnj --pbsnj-focal EA
```

**If you have already run `fly_fst.py`, you do not need to read your chromosome sequence files again.** Run:

```bash
python pbsnj.py --outdir fly_fst_output --focal EA
```

For a multipopulation comparison with all seven populations, `--populations` may be omitted (the default uses all populations in `run_config.json`). To restrict to an evolutionary subset, use e.g. `--populations EA,ED,SP,ZI,Ug`; at least three populations are required. To obtain terminal branches for every population in one run use `--focal all` (or `--pbsnj-focal all` in the full pipeline). You do not need a focal/nonfocal triplets file for PBSnj.

Options available for the standalone script:

```text
--min-window-snps 10                Minimum polymorphic SNP count in at least one pair
--min-window-callable-fraction 0.5  Required in *each* population pair
--outlier-quantile 0.99             Descriptive top 1% threshold
--epsilon 1e-10                     Saturation cap for FST >= 1 in T conversion
--no-plots                          Produce tables without plots
```

### PBSnj computations and statistical cautions

For each genomic window, the script assembles the full symmetric pairwise FST matrix among the selected populations and calculates pairwise divergence `T_ij = -ln(1 - clip(FST_ij, 0, 1-epsilon))`. Standard **Saitou–Nei neighbor joining** estimates an unrooted tree from this matrix. The focal terminal-branch length is recorded as `PBSnj_EA`. All NJ branches are estimated independently for each window, and the unrooted NJ topology may change between windows. For three populations the terminal branch equals the ordinary three-population PBS formula.

For every population pair, a window must have sufficient pairwise callable bases (`--min-window-callable-fraction`). Windows with undefined pairwise FST are omitted. A fully callable *invariant pair* (0 informative SNPs) receives pairwise distance 0. To preserve focal-private polymorphism, the `--min-window-snps` filter requires the count in **at least one** pair, not in every pair. FST=1 is capped at `1-epsilon` before logarithm, so highly divergent windows can be very sensitive to `epsilon`. Negative raw FST is mapped to 0 **only for PBSnj**, and negative *NJ terminal branches* are **retained**, not clipped (they indicate non-additivity/noisy estimates and must not be interpreted as negative time). All original pairwise FST files are unchanged.

The calculation requires the complete distance matrix (all pairs) for each window; do not interpret missing pairs as zero. Pairwise FST across many populations can be noisy for low sample sizes. Recombination-rate variation, inversions, low diversity, and the different sample sizes may inflate rankings. Outlier windows are **descriptive upper-tail rankings, not selection p-values**. Compare different window sizes, account for missing/callability biases and inversions, and assess genomic enrichment with appropriate linkage-aware null models.

### PBSnj output

- `pbsnj/w50000_s10000/pbsnj.tsv`: per-window EA terminal length (`PBSnj_EA`), QC, n pairs, number of saturated pairwise FST estimates, and other requested focal branches.
- `pbsnj/w50000_s10000/EA.top1pct.tsv`: top 1% eligible EA windows (empirical threshold).
- `pbsnj/w50000_s10000/EA.pbsnj.pdf`: EA branch scan across chromosome arms.
- `pbsnj/summary.tsv`: eligible window counts and thresholds per window size / focal population.
- `pbsnj/methods.json`: filters and numerical conventions.

Test implementation:

```bash
python -m unittest test_pbsnj.py -v
```

Tests include: exact three-population PBS equivalence, recovery of terminal lengths from known seven-population additive distances, invariant relabeling, distance transforms, and compatibility with previously generated FST windows.

---

## EA ion-homeostasis gene-set enrichment with genomic block permutations

**New module: `gene_set_enrichment.py`.** Runs on an existing `pbsnj.tsv`; **no genome reprocessing needed**. It asks whether prespecified ion-related gene sets overlap unusually elevated EA PBSnj signals compared with a chromosome-arm block-permutation positional null. This is an *exploratory enrichment test*, not direct evidence of selection or an ion GWAS.

### Required files

1. **PBSnj windows** from a previous run, e.g. `fst_results/pbsnj/w50000_s10000/pbsnj.tsv` with `PBSnj_EA` and `passes_filter` columns.
2. **Gene coordinates** aligned to your `.seq` files: chromosome-arm, 0-based inclusive start, 0-based exclusive end, stable gene ID, optional display name. Example BED row: `2L<TAB>1000<TAB>2000<TAB>FBgn0000001<TAB>GENE_SYMBOL`. Alternatively, supply **GFF3** with feature type `gene` (coordinates are converted from 1-based inclusive), matching the same dm3/BDGP5 reference as the `.seq` files. **Do not use dm6/r6 GFF3 without liftover or validated coordinate conversion.**
3. **Predeclared gene sets** in either GMT (`set_name<TAB>description<TAB>gene1<TAB>gene2...`) or two-column TSV (`set_name<TAB>gene_id`, one association per line). Gene identifiers must match `gene_id` in the annotation (default GFF3 `ID` attribute). Gene symbols in gene sets will *not* automatically match FlyBase `FBgn...` IDs—convert the set identifiers or use a matching annotation ID attribute.

The supplied `example_ion_gene_sets.tsv` is only a tiny **illustrative schema** with gene-symbol identifiers, not a curated, complete or ready-to-test ion-homeostasis ontology. The supplied `example_genes.bed` has **invented coordinates for formatting examples only** and must **never** be used for biological analysis. Assemble biological gene sets from FlyBase/GO (e.g., manganese ion transport, calcium ion transport, iron ion homeostasis, zinc ion transport, sodium/potassium transport, phosphate transport, metal ion response, oxidative stress). Decide the gene sets **before** inspecting PBSnj peaks.

### Standalone run on existing PBSnj outputs

```bash
python gene_set_enrichment.py \
  --pbsnj-tsv fst_results/pbsnj/w50000_s10000/pbsnj.tsv \
  --genes dm3_genes.gff3 --gene-format gff3 \
  --gene-sets ion_gene_sets.gmt \
  --outdir fst_results/enrichment/w50000_s10000 \
  --focal EA --score max \
  --permutations 10000 --block-bp 250000 --seed 2026 \
  --min-genes 5
```

For sensitivity analysis, repeat with `--block-bp 100000`, `500000`, and `1000000`, and with `--score mean`. Nonoverlapping 50-kb windows are preferable as an additional robustness check because the primary scan's overlapping windows are strongly correlated.

### Integrated run (new FST calculation)

```bash
python fly_fst.py \
  --seq-dir /data/flyseq --metadata metadata.csv --outdir fst_results \
  --arms 2L,2R,3L,3R,X --windows 50000:10000,100000:25000 \
  --pbsnj --pbsnj-focal EA \
  --enrich-genes dm3_genes.gff3 \
  --enrich-gene-sets ion_gene_sets.gmt \
  --enrich-window 50000:10000 \
  --enrich-permutations 10000 --enrich-block-bp 250000
```

### Statistical procedure and interpretation

- **Gene score:** highest (`--score max`, default) or mean (`--score mean`) finite PBSnj value among windows *overlapping the gene interval*. Long genes can overlap more windows and therefore have higher maxima; a fixed gene interval is rescored in each permutation to partially control that bias. Genes without a scored window are excluded.
- **Enrichment score:** difference between the mean rank percentile of the set's genes and that of the remaining eligible genes. Larger positive values correspond to elevated focal branch differentiation at genes in the set. The original genome-wide ranked gene universe is reranked independently for each permutation.
- **Null:** divide chromosome arms into nonoverlapping physical blocks (`--block-bp`, default 250 kb); shuffle the *contiguous sequences of finite window scores* across blocks **within each arm**, leaving the actual gene intervals and low-callability/missing-score positions fixed. This approximately retains autocorrelation within blocks and maintains gene length/position and arm distribution, but does not preserve correlations across block boundaries. Truncated/missing-score blocks may have different numbers of windows; score sequences are concatenated and remapped in physical order, so block boundaries may shift slightly.
- **One-sided empirical p:** `(1 + number of permuted scores >= observed) / (1 + number of permutations)`. `q_BH` is the Benjamini–Hochberg correction across eligible gene sets. This assumes a useful positional-exchangeability null **within chromosome arms**; block shuffling is **not** a universally calibrated model of neutral genome evolution.
- **Critical limitations:** The procedure does not match local recombination rate, inversion state, centromeric/telomeric proximity, annotation density or background selection. Different functional categories are often clustered on chromosomes; significant enrichment can reflect this structure. Mask known inversions/low-recombination regions or add recombination-stratified analyses before strong adaptive claims. Use gene-set and block-size sensitivity analyses; avoid post hoc gene-set selection. A high PBSnj value can also be due to drift and demographic history.

### Outputs

- `enrichment.tsv`: set size, tested genes, mean percentile, effect (difference in percentile ranks), permutation-null moments, one-sided empirical p, and BH q.
- `gene_scores.tsv`: arm, interval, stable ID/name, length, gene PBSnj score and rank, tested memberships.
- `excluded_sets.tsv`: gene sets failing mapping/minimum-gene criteria.
- `methods.json`: parameters, scored-gene counts and statistical conventions.

### Validation

```bash
python -m unittest test_enrichment.py -v
```

Tests cover boundary overlaps, gene scoring, arm-specific block score preservation, multiple-testing correction, and an end-to-end synthetic enrichment run. Biological gene coordinates and real ion gene sets are **not provided** by the script; you must supply assembly-matched annotations.
- `enrichment.pdf`: top gene sets ranked by empirical permutation p-value with BH q-value labels (if Matplotlib is installed).

---

## Population-agnostic workflow (updated)

**This release has no hard-coded EA/ED/SP/ZI/UG/FR/CO prefix list.**
`make_metadata.py` infers the letters immediately preceding the first digit in each sample name:
`CK1 -> CK`, `EA49 -> EA`, `ZI90 -> ZI`, `Ug5 -> UG` and `Fr217 -> FR`.
Inferred labels default to UPPERCASE, but are **case-sensitive in downstream commands**. Explicit sample mappings are preserved as written. Inspect the output CSV before running.

```bash
python make_metadata.py --seq-dir chr2R --output chr2R/metadata.csv
```

For names with other conventions, supply mapping files:

```csv
# sample_map.csv (remove this comment line)
sample,population
Sample_A,CK
Sample_B,EA
```

```bash
python make_metadata.py --seq-dir chr2R --output chr2R/metadata.csv \
  --sample-map sample_map.csv
```

Or a `prefix,population` mapping file with `--prefix-map prefix_map.csv`.
You can also use `--prefix-regex '^([A-Za-z]+)[_-]'` for names like `EA_49`.
Unknown naming schemes generate an error rather than silently assigning the wrong population (or `--unknown skip` to omit them). Sample IDs cannot contain unknown population prefixes without explicitly mapping them.

### FST: select any two or more populations

```bash
python fly_fst.py --seq-dir chr2R --metadata chr2R/metadata.csv \
  --outdir chr2R_EA_subset --arms 2R --windows 10000:5000 \
  --populations EA,ED,SP,ZI,UG,FR,CO --no-sites
```

`--populations` controls which samples enter the initial FST and diversity calculations. It avoids unnecessary pairwise calculations for unused populations. Default: all populations present in metadata. Populations must have **at least two** independent haploid sequences per arm; otherwise the program gives a clear error. `--allow-incomplete-arms` allows missing files, but it will not override the minimum sample-size requirement.

### Integrated PBSnj on an arbitrary set of three or more populations

```bash
python fly_fst.py --seq-dir chr2R --metadata chr2R/metadata.csv \
  --outdir chr2R_EA_subset --arms 2R --windows 10000:5000 \
  --populations EA,ED,SP,ZI,UG,FR,CO \
  --pbsnj --pbsnj-focal EA \
  --pbsnj-populations EA,ED,SP,ZI,UG,FR,CO
```

For an existing FST run with a larger population set:

```bash
python pbsnj.py --outdir chr2R_all_population_results \
  --focal EA --populations EA,ED,SP,ZI,UG,FR,CO
```

PBSnj does **not** silently remove populations with poor coverage. All selected pairwise distances must be finite and callable. The output now includes `lowest_callability_pair`, `lowest_callability_fraction`, and `min_pair_callable_sites` for `low_callability` windows, helping diagnose which population pair fails.

### PBS/PBE triplets from any subset

```bash
python branch_stats.py --outdir chr2R_all_population_results \
  --populations EA,ED,SP,ZI,UG,FR,CO
```

Or during the full run use `--pbs-pbe --branch-populations EA,ED,SP,ZI,UG,FR,CO`. For a few evolutionary contrasts, use `--triplets triplets.csv` (standalone) or `--branch-triplets triplets.csv` (integrated). Note: standalone PBSnj/PBS/PBE reruns **overwrite their corresponding output files**; use separate FST output directories to preserve multiple comparison sets.

Population subsets must use the exact case present in your metadata; the pipeline does not use substring/partial matching. The choice of populations determines both power and phylogenetic interpretation of PBSnj, not just computational load.

### Python 3.7 and gzip enrichment compatibility (patch)

The updated `gene_set_enrichment.py` is compatible with Python 3.7: it does not
use `str.removeprefix`, accepts `.gff.gz`/`.gff3.gz` and `.bed.gz` annotations,
can read `.gmt.gz`, and supports older Matplotlib tick-label APIs.

You can replace just `gene_set_enrichment.py` in an existing installation and
rerun the enrichment step from a previously generated PBSnj output without
rerunning WGS processing.

## PBSnj gene annotations and peak labels (updated)

For standalone PBSnj reruns from **existing FST window tables**:

```bash
python fly_fst_pipeline/pbsnj.py \
  --outdir chr2R_PBSnj_out --focal EA \
  --populations EA,ED,SP,ZI,UG,FR,CO \
  --genes dmel-all-no-analysis-r5.57.gff.gz \
  --peak-labels 5 --peak-min-separation-bp 100000
```

This regenerates `pbsnj/w10000_s5000/pbsnj.tsv`, `EA.top1pct.tsv`,
and `EA.pbsnj.pdf` for each configured window resolution; it does **not**
reprocess sequence files or repeat the FST calculation.

In an integrated `fly_fst.py` run, use `--pbsnj-genes` with the same GFF3.
If `--enrich-genes` is supplied, `fly_fst.py` automatically uses it to
annotate the PBSnj tables/plots as well, even without `--pbsnj-genes`.

Added TSV fields:
- `n_overlapping_genes`, `overlapping_gene_ids`, `overlapping_gene_names`:
  all gene intervals intersecting each half-open analysis window.
- `nearest_gene_id`, `nearest_gene_name`, `nearest_gene_distance_bp`:
  nearest gene to the window **midpoint**, with 0 bp when the midpoint
  falls within the gene. Overlapping genes at the window edge may differ
  from the nearest-to-midpoint gene.
- `nearest_gene_description`, `nearest_gene_biotype`: optional GFF3
  attributes when supplied; availability depends on the FlyBase release.

Plot labels are the nearest-to-midpoint gene for the highest-scoring
window in each merged top-quantile peak. Overlapping top-percentile windows
are collapsed into peak regions; up to five highest peaks per arm are
labeled by default, requiring >=100 kb between label positions. A labeled
gene is **positional context, not proof of causality**.

The original PBSnj branch statistic and empirical top-percentile definitions
are unchanged. To reproduce custom PBSnj thresholds during standalone
reprocessing, pass the same `--min-window-snps`,
`--min-window-callable-fraction`, and `--outlier-quantile` settings used in
the original analysis. GFF3 chromosome arms must match the FST window
naming (e.g. `2R` or `chr2R`) and the same dm3 reference assembly.

## New: continuous-trait GWAS

The pipeline supports **`--continuous-trait-gwas`** with **one independent association analysis per numeric trait**. Full implementation, tested example commands, interpretation and outputs are documented in [CONTINUOUS_TRAIT_GWAS_README.md](CONTINUOUS_TRAIT_GWAS_README.md). For existing PBSnj results, run `python fly_fst_pipeline/continuous_trait_gwas.py` separately to avoid recomputing population differentiation.
