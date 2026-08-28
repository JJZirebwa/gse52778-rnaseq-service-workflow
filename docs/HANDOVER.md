# Service handover

## Question

This exercise asks which genes differ between dexamethasone-treated and untreated airway smooth-muscle samples in public dataset GSE52778, while controlling for the four donor/cell-line pairs.

## Intake and design

- Input: processed gene-count data from the Bioconductor `airway` package, not FASTQ files.
- Samples: eight libraries, comprising four complete untreated/treated pairs.
- Model: `~ cell + dex`.
- Reference: `dex = untrt`.
- Reported contrast: `dex`, treated (`trt`) versus untreated (`untrt`), so a positive log2 fold-change indicates higher abundance after treatment.
- Pre-model filter: retain genes with a count of at least 10 in at least four samples: `rowSums(counts >= 10) >= 4`.
- Multiple testing: Benjamini-Hochberg adjusted p-values, with the filtered handoff table using `padj <= 0.05`.

## Delivered outputs

- `validated_samples.csv` and `intake_validation.json`: accepted intake and machine-readable validation.
- `library_sizes.png` and `count_distributions.png`: processed-count QC summaries.
- `pca.png` and `sample_distance_heatmap.png`: variance-stabilised sample relationships.
- `full_results.csv`: all tested genes and DESeq2 result fields.
- `fdr_0.05_results.csv`: rows with adjusted p-value at or below 0.05.
- `analysis_metadata.json` and `result_validation.json`: explicit model/filter contract and automated checks.
- `sessionInfo.txt`: complete R session and loaded-package versions.
- `provenance_manifest.json`: accession, commit, container references, software versions and SHA-256 hashes.
- `LIMITATIONS.md`: interpretation boundaries.
- `nextflow/`: runtime report, timeline, trace and DAG when the launcher is used.

## Interpretation and sanity check

The full table estimates the dexamethasone effect after accounting for cell-line differences. CRISPLD2 is required to have a positive log2 fold-change for the stated treated-versus-untreated contrast. That direction check can expose a reversed contrast or gross analysis error; it does not validate the entire model and is not evidence that one gene explains the biological response.

## Suggested next steps

Review the QC plots and sample annotations before interpreting genes or pathways. For a real service request, agree the biological question and decision threshold with the requester, review influential samples and model diagnostics, then consider pathway-level summaries and an independent dataset. Any raw-read reprocessing would be a separately scoped analysis.

