# Separate GO ontology enrichment (BP / MF / CC)

The official GO namespaces are **Biological Process**, **Molecular Function**,
and **Cellular Component** (not Cellular Process).

Provide the *same* `go-basic.obo` used to create `flybase_all_go.gmt`:

```bash
python fly_fst_pipeline/gene_set_enrichment.py \
  --pbsnj-tsv PBSnj_out/pbsnj/w10000_s5000/pbsnj.tsv \
  --genes dmel-all-no-analysis-r5.57.gff \
  --gene-sets flybase_all_go.gmt \
  --go-obo go-basic.obo \
  --outdir PBSnj_out/enrichment/w10000_s5000 \
  --focal EA --permutations 10000 --block-bp 100000 \
  --min-genes 5 --score max
```

Or add `--enrich-go-obo go-basic.obo` to the full `fly_fst.py` command.

In the enrichment output directory:
- `enrichment.tsv`, `enrichment.pdf`: original combined GO results
- `gene_scores.tsv`: unchanged shared gene scores
- `biological_process/enrichment.tsv`, `enrichment.pdf`, `excluded_sets.tsv`, `methods.json`
- `molecular_function/` with the same files
- `cellular_component/` with the same files
- `unclassified_go_sets.txt`: set names not resolvable using the supplied OBO

Each namespace table has `q_BH` adjusted over only that namespace and
`q_BH_global` adjusted over all eligible gene sets. Empirical permutation
P values and gene scores are shared: permutations are calculated once.
GO set names must be GO IDs (GO:nnnnnnn), as produced by the FlyBase converter.

Python 3.7 compatible. Vector PDFs embed TrueType fonts for Illustrator.
