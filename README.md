# Reproducing an RNA-seq differential expression analysis: APOE4 in human iPSC-derived cerebral organoids

Work in progress. Reproducing the organoid RNA-seq analysis from Meyer et al. (2019), *Cell Reports* 26, 1112–1127 (GEO: GSE117588), which compared isogenic APOE4 and APOE3 organoids.

## Roadmap

- [ ] Reproduce the published edgeR result from the authors' count table (target: 16,581 expressed genes, 5,141 at FDR < 0.05)
- [ ] Independent check with DESeq2 (R and PyDESeq2); Benjamini–Hochberg implemented and verified
- [ ] Raw reads to counts with Salmon; comparison with the authors' and NCBI's counts
- [ ] Reproducible Snakemake workflow, tests and fresh-clone check
- [ ] Robustness: how much the result depends on any single sample