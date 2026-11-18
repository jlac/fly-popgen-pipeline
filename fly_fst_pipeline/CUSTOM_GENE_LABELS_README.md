# Controlling PBSnj labels (vector PDF output)

These new settings change **plot annotations only**. The PBSnj values, `pbsnj.tsv`, top-percentile TSV, and enrichment statistics remain unchanged. Gene symbols and FlyBase `FBgn` IDs come from a GFF3 (compressed GFF3 supported) or BED file matching your reference genome. All plot outputs remain vector PDFs.

## Regenerate only PBSnj output (no genome sequence reprocessing)

```bash
python fly_fst_pipeline/pbsnj.py \
  --outdir PBSnj_out --focal EA \
  --populations EA,ED,SP,ZI,UG,FR,CO \
  --genes dmel-all-no-analysis-r5.57.gff.gz \
  --plot-label-quantile 0.95 \
  --peak-labels 15
```

`--plot-label-quantile` controls *automatic plot labels only*; e.g., `0.95` labels from the top 5% of valid PBSnj windows (subject to `--peak-labels` and separation). The original `--outlier-quantile` still controls top-percentile TSV creation.

Instead of a percentile threshold, set an absolute one:

```bash
python fly_fst_pipeline/pbsnj.py \
  --outdir PBSnj_out --focal EA --genes dmel-all-no-analysis-r5.57.gff.gz \
  --plot-label-min-score 0.25 --peak-labels 25
```

Or explicitly request genes, regardless of their PBSnj ranks:

```bash
python fly_fst_pipeline/pbsnj.py \
  --outdir PBSnj_out --focal EA --genes dmel-all-no-analysis-r5.57.gff.gz \
  --plot-label-genes 'heph,Drak,Ir11a,FBgn0000008' \
  --plot-manual-only
```

Alternatively, put gene symbols or FBgn IDs in a text file, one per line (comments beginning with `#` are ignored):

```bash
python fly_fst_pipeline/pbsnj.py \
  --outdir PBSnj_out --focal EA --genes dmel-all-no-analysis-r5.57.gff.gz \
  --plot-label-genes-file ion_gene_labels.txt \
  --plot-manual-only
```

The gene's arrow points to the **highest-scoring finite PBSnj window overlapping its annotated interval**, not its transcription start site. Genes with no overlapping passing window cannot be plotted; nonexistent gene names are warned about. Requested genes bypass the automatic peak count and distance filters; automatic labels are suppressed near manually labeled genes to avoid duplicates. A gene's name might be duplicated in the GFF; use a FBgn ID to target a specific gene.

### Equivalent options in `fly_fst.py`

When running the full analysis, the options are prefixed with `--pbsnj-`:

| PBSnj-only (`pbsnj.py`) | Full pipeline (`fly_fst.py`) |
|---|---|
| `--plot-label-quantile 0.95` | `--pbsnj-plot-label-quantile 0.95` |
| `--plot-label-min-score 0.25` | `--pbsnj-plot-label-min-score 0.25` |
| `--plot-label-genes 'heph,Drak'` | `--pbsnj-plot-label-genes 'heph,Drak'` |
| `--plot-label-genes-file ion_gene_labels.txt` | `--pbsnj-plot-label-genes-file ion_gene_labels.txt` |
| `--plot-manual-only` | `--pbsnj-plot-manual-only` |
| `--peak-labels 15` | `--pbsnj-peak-labels 15` |
| `--peak-min-separation-bp 50000` | `--pbsnj-peak-min-separation-bp 50000` |

If both `--plot-label-min-score` and `--plot-label-quantile` are specified, the absolute threshold takes precedence. The default top-percentile TSV threshold remains 0.99 unless you explicitly set `--outlier-quantile` on `pbsnj.py` or `fly_fst.py`.

**Important:** `pbsnj.py` recomputes window PBSnj from existing pairwise FST TSVs; it does not read previously calculated `pbsnj.tsv` directly. It is far faster than repeating the sequence-based `fly_fst.py` analysis, but will overwrite the PBSnj figures and tables under the existing `outdir`.
