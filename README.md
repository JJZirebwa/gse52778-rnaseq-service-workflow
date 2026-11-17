# GSE52778 paired RNA-seq workflow

[![CI](https://github.com/JJZirebwa/gse52778-rnaseq-service-workflow/actions/workflows/ci.yml/badge.svg)](https://github.com/JJZirebwa/gse52778-rnaseq-service-workflow/actions/workflows/ci.yml)

I wanted to see what a small differential-expression request looks like when it runs as a repeatable workflow, from intake checks to a bundle someone else can pick up. This one takes the eight processed-count samples in GSE52778 (airway smooth muscle cells with and without dexamethasone), runs a paired DESeq2 analysis in Nextflow DSL2, checks the results against the contract below and packages them with a provenance manifest.

## Analysis contract

- Dataset: GSE52778 through Bioconductor `airway` 1.26.0.
- Samples: exactly eight, representing four cell-line pairs with one untreated and one dexamethasone-treated sample each.
- Model: `~ cell + dex`.
- Reference and contrast: `untrt` is the reference; results are `trt` versus `untrt`.
- Low-count rule: retain genes satisfying `rowSums(counts >= 10) >= 4`.
- FDR table: Benjamini-Hochberg `padj <= 0.05`.
- Sanity check: CRISPLD2 must have positive log2 fold-change for treated versus untreated.

The workflow starts from processed counts, so there is no FASTQ QC or alignment step.

## Prerequisites and one-command run

Install Docker and Nextflow 25.04.6 or a compatible newer 25.x release. The runtime is built from the official Bioconductor 3.20 image pinned by digest, with DESeq2 1.46.0 and airway 1.26.0 checked during the build.

From the repository root, run:

```bash
./scripts/run_workflow.sh
```

The launcher refuses to mix with an existing `results/` directory. To choose clean output and work locations:

```bash
./scripts/run_workflow.sh results-run-1 work-run-1
```

The four logical service stages are:

1. `VALIDATE_SERVICE_INTAKE`
2. `RUN_DESEQ2_PAIRED`
3. `VERIFY_RESULTS`
4. `BUILD_HANDOFF_BUNDLE`

Final deliverables are written to `results/handoff/`. Nextflow runtime evidence is written to `results/nextflow/`.

## Validator tests

The validators use only the Python standard library:

```bash
python3 -m unittest discover -s tests -v
```

## Determinism check

Run twice into separate clean locations and compare the primary table:

```bash
./scripts/run_workflow.sh results-run-1 work-run-1
./scripts/run_workflow.sh results-run-2 work-run-2
sha256sum results-run-1/handoff/full_results.csv results-run-2/handoff/full_results.csv
```

Both SHA-256 values must be identical. On macOS, use `shasum -a 256` in place of `sha256sum`.

## Outputs

[docs/HANDOVER.md](docs/HANDOVER.md) lists every output in the bundle, and [docs/INTERPRETATION.md](docs/INTERPRETATION.md) covers what to keep in mind when reading them. The provenance manifest hashes every key handoff file and records the dataset accession, Git commit, container references and software versions.

## Primary sources

- [Bioconductor airway package](https://bioconductor.org/packages/3.20/data/experiment/html/airway.html)
- [Bioconductor DESeq2 package](https://bioconductor.org/packages/3.20/bioc/html/DESeq2.html)
- [Official Bioconductor Docker documentation](https://bioconductor.org/help/docker/)
- [NCBI GEO GSE52778 record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778)

## Licence

Code in this repository is released under the MIT Licence. The upstream dataset and R packages retain their own terms and citations.
