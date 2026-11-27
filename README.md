# Reproducing an RNA-seq differential expression analysis: APOE4 in human iPSC-derived cerebral organoids

Work in progress. This project independently reproduces the organoid RNA-seq analysis from Meyer et al. (2019), *Cell Reports* 26, 1112–1127 (GEO: [GSE117588](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117588)), first from the authors' count table and then from the raw sequencing reads, and tests how robust the published result is.

## Roadmap

- [x] Reproduce the published edgeR result from the authors' count table (target: 16,581 expressed genes, 5,141 at FDR < 0.05)
- [ ] Independent check with DESeq2 (R and PyDESeq2); Benjamini–Hochberg implemented and verified
- [ ] Raw reads to counts with Salmon; comparison with the authors' and NCBI's counts
- [ ] Reproducible Snakemake workflow, tests and fresh-clone check
- [ ] Robustness: how much the result depends on any single sample

## The study

APOE4 is the strongest common genetic risk factor for late-onset Alzheimer's disease. Meyer et al. took an iPSC line carrying APOE4 and gene-edited it to APOE3, the common neutral variant. Because both lines come from the same person, they are isogenic: differences between them should mostly reflect APOE itself rather than different genetic backgrounds. Both lines were grown into cerebral organoids and profiled by RNA sequencing, three samples per genotype.

The authors aligned reads with STAR to GRCh38 (Ensembl release 92), kept genes with at least 1 count per million in at least 3 samples, normalised with TMM, and tested for differential expression with edgeR's likelihood ratio test. They report 5,141 genes at FDR < 0.05 (Table S6).

## Week 1: reproducing the published edgeR result

**Notebooks:** `notebooks/01_inspect_counts.ipynb` (Python) and `notebooks/02_edger_reproduction.ipynb` (R)

### Checking the data first

Before any analysis, I checked the authors' count table:

- 58,395 genes × 6 samples, matching the number of annotated genes stated in the paper
- Ensembl gene IDs with no version suffixes, so they match Table S6 directly
- whole-number counts only, as edgeR requires
- library sizes of 9.1–13.0 million counted fragments per sample
- 23,791 genes (about 41%) with no reads in any sample, as expected for a single tissue type

The count columns use the authors' sample names (ADR_5 to ADR_14), not GEO accessions, and they aren't ordered by genotype. I matched each column to its genotype through the GEO sample titles, recorded in `metadata/samples.csv`, and checked the matching in code rather than relying on column order. An early version of my code matched samples by position, which silently assigned three samples to the wrong group; matching by name fixed it.

### Results

| Checkpoint | Paper | This reproduction |
|---|---|---|
| Genes after filtering | 16,581 | **16,581** |
| Significant genes (FDR < 0.05) | 5,141 | 5,127 |
| Higher in APOE4 | 2,457 | 2,454 |
| Lower in APOE4 | 2,684 | 2,673 |

**Filtering matched exactly,** confirming the data, sample labels and filtering rule were all reproduced correctly.

**Fold changes were identical.** Across all 16,581 genes, my log fold changes and the authors' have a Pearson correlation of 1.000.

![Log fold changes: this reproduction against the published values](images/logfc_comparison.png)

**Significant genes overlapped almost completely.** 5,102 genes were significant in both analyses: 99.2% of the paper's significant genes, and 99.5% of mine. The 64 genes that differ (39 significant only in the paper, 25 only in my analysis) all have FDRs between 0.045 and 0.056 in both analyses. They sit on the 0.05 threshold and tip to one side or the other.

**The likely cause is the software version:** the authors used edgeR 3.20.9 (2017), and I used edgeR 4.10.7. Since fold changes are identical, the small difference must come from p-values, which depend on dispersion estimation, an area that has been updated across edgeR versions. Running the original version would confirm this definitively.

### What the data look like

- **Replicates are consistent:** the biological coefficient of variation is 0.165, meaning replicates typically differ by about 16.5% in expression.
- **The genotypes separate clearly.** In the MDS plot, the first dimension explains 75% of the variation and completely separates APOE4 from APOE3. A second, smaller dimension (9%) separates one APOE3 sample, ADR_6, from the other two, which I'll revisit in the robustness analysis.
- **TMM normalisation factors** range from 0.92 to 1.05, with all APOE4 samples near 1.05 and all APOE3 samples below 1, consistent with the large number of genes that differ between the groups.

![MDS plot of the six samples, coloured by genotype](images/mds_plot.png)

## Week 2: an independent check with DESeq2 (in progress)

**Notebook:** `notebooks/03_deseq2.ipynb` (R)

Week 1 asked whether I could get the same answer as the authors. Week 2 asks whether the answer holds with a different statistical method. If a result appears with only one tool, it could be a quirk of that tool; if two independent methods agree, it's more likely to be real.

I ran DESeq2 (1.50.2) on exactly the same 16,581 genes as edgeR, so that only the method differs. To keep the gene sets identical, I recalculated the authors' filter (at least 1 count per million in at least 3 samples) and confirmed it gives the same 16,581 genes. I turned off DESeq2's outlier detection and independent filtering, so that every gene receives a result, as in edgeR.

### The two methods mostly agree

| | edgeR (likelihood ratio test) | DESeq2 (Wald test) |
|---|---|---|
| Significant genes (FDR < 0.05) | 5,127 | 4,548 |
| Higher in APOE4 | 2,454 | 2,124 |
| Lower in APOE4 | 2,673 | 2,424 |

- **4,531 genes are significant in both methods:** 99.6% of DESeq2's significant genes, and 88.4% of edgeR's.
- **The two never disagree on direction:** every gene that goes up in one goes up in the other.
- **The disagreement is one-sided:** 596 genes are significant only in edgeR, and just 17 only in DESeq2.

### Finding the cause by elimination

A differential expression analysis has several parts, and a difference in any of them could explain the gap. I checked each in turn:

| Possible cause | Check | Result |
|---|---|---|
| Different genes | Same filtered gene set in both | Identical: **ruled out** |
| Normalisation | DESeq2 size factors against edgeR's TMM-adjusted library sizes, on the same scale | Within 1% for every sample: **ruled out** |
| Fold changes | Correlation of estimated changes, gene by gene | Pearson r = 0.9997: **ruled out** |
| The statistical test | DESeq2 rerun with a likelihood ratio test, like edgeR's | 4,545 genes, against 4,548 with the Wald test: **ruled out** |
| Dispersion (how variable each method thinks each gene is) | Compared directly, below | **The cause** |

![edgeR and DESeq2 fold changes, gene by gene](images/edgeR_DESeq2_comparison.png)

### The genes the methods disagree on

| Group | Genes | Median expression (baseMean) | Median size of change (\|log2FC\|) |
|---|---|---|---|
| Significant in both | 4,531 | 243 | 0.91 |
| edgeR only | 596 | **34** | 0.81 |
| DESeq2 only | 17 | **1,349** | **0.38** |

The genes only edgeR calls significant are **lowly expressed,** about seven times lower than the genes both methods agree on, while their changes are almost as large. The few genes only DESeq2 calls significant are the opposite: highly expressed genes with small changes, which many reads allow to be measured precisely.

### Confirming the cause: dispersion

Both methods model counts with the negative binomial distribution, in which **dispersion** measures how much a gene varies between replicates beyond counting noise. With only three replicates per group, each gene's own estimate is unreliable, so both methods borrow information across genes, pulling each gene towards a trend for genes of similar expression. They do this differently: edgeR pulls genes firmly towards the trend, while DESeq2 leaves genes that look unusually variable at their own, higher estimate.

I compared the two methods' dispersion estimates for every gene:

| Group | Median ratio of dispersions (DESeq2 ÷ edgeR) |
|---|---|
| All genes | 1.08 |
| Significant in both | 1.06 |
| **edgeR only** | **1.97** |
| DESeq2 only | 0.93 |

![Gene-by-gene dispersion estimates from edgeR and DESeq2](images/logscale_comparison_by_group.png)

*Dispersion estimates from edgeR (x) and DESeq2 (y), log scales; the dashed line marks equal estimates. The 596 genes significant only in edgeR (orange) sit mostly above the line: DESeq2 estimates them as about twice as variable, so it is less confident in their changes. The 17 genes significant only in DESeq2 (green) are stable, highly expressed genes where DESeq2 estimates slightly lower variability.*

### Conclusion

The published result is robust to the choice of method: edgeR and DESeq2 agree on the direction of every change, and on more than 99% of DESeq2's significant genes. Where they differ, it is in lowly expressed genes, where the evidence is thinnest and the two methods' different ways of estimating variability lead to different levels of confidence. Neither is simply right: edgeR is more sensitive for low-count genes, and DESeq2 more cautious.

### Still to do

- PyDESeq2, a Python re-implementation of DESeq2, compared with the R version
- Benjamini–Hochberg correction implemented from scratch, and checked against edgeR's FDR values

## Data

The data are not included in this repository. To reproduce the analysis, download into `data/`:

- `GSE117588_organoids-gene-counts.tsv.gz`, the authors' count table, from the [GEO page](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE117588)
- Table S6 from the paper's supplementary materials, saved as `meyer2019_tableS6_organoid_DEGs.xlsx`

Sample information is in `metadata/samples.csv`.

## Software

R with edgeR 4.10.7 and readxl; Python with pandas. The environment is managed with conda.

## Reference

Meyer K, et al. (2019). REST and neural gene network dysregulation in iPSC models of Alzheimer's disease. *Cell Reports* 26, 1112–1127. https://doi.org/10.1016/j.celrep.2019.01.023

AI tools were used for code review and debugging.