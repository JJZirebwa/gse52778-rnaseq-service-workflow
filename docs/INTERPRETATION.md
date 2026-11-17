# Reading the results

- The input is the processed counts from the Bioconductor `airway` package, so there is no FASTQ QC, alignment or re-quantification step.
- Eight samples from four cell-line pairs is a small study. The paired model removes the baseline differences between cell lines, but the estimates for weaker effects will move with more data.
- The low-count filter (`rowSums(counts >= 10) >= 4`) is a choice. A different filter changes which genes are tested and how many tests the FDR correction covers.
- DESeq2 estimates depend on its negative-binomial model and the `~ cell + dex` design; any covariate not in the design is not accounted for.
- An adjusted p-value controls the false discovery rate. It says nothing about effect size or biological importance, so read it alongside the log2 fold-change.
- The CRISPLD2 check only confirms the contrast runs in the expected direction. It catches a flipped contrast, not a wrong model.
