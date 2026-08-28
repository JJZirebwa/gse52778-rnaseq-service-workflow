# Limitations and interpretation boundaries

- This workflow starts from the processed counts distributed in the Bioconductor `airway` package. It does not perform or claim FASTQ quality control, read alignment, transcript quantification or re-quantification against a current reference.
- The dataset contains eight samples from four cell-line/donor pairs. The paired model controls those observed baseline differences, but the sample size remains small and does not establish broad population generalisability.
- The low-count filter is a declared pragmatic rule for this exercise. Different pre-specified filters can change the tested gene set and multiple-testing burden.
- DESeq2 estimates are conditional on its negative-binomial model and the supplied design. Unmodelled technical or biological covariates may remain.
- Adjusted p-values control the false discovery rate under the analysis assumptions; they do not measure biological importance or guarantee replication.
- The positive CRISPLD2 direction is a narrow sanity check for the stated contrast, not circular evidence that the whole analysis is correct.
- This is a compact service-workflow exercise using a public bulk RNA-seq dataset. It is not a production deployment, a single-cell analysis, a clinical assay or clinically validated evidence.
- Handover outputs require scientific review before any downstream interpretation. Independent data and, where relevant, orthogonal experimental validation are appropriate next steps.

