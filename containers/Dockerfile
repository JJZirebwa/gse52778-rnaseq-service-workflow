FROM bioconductor/bioconductor_docker:RELEASE_3_20@sha256:5053f59de1c1a79b5a05d820ba722bae62ded0c1042e547189aa0226f111ea7e

LABEL org.opencontainers.image.title="GSE52778 paired RNA-seq service workflow"
LABEL org.opencontainers.image.description="Pinned DESeq2 1.46.0 and airway 1.26.0 runtime"
LABEL org.opencontainers.image.source="https://github.com/JJZirebwa/gse52778-rnaseq-service-workflow"

RUN R --vanilla -e 'BiocManager::install(c("DESeq2", "airway"), version = "3.20", ask = FALSE, update = FALSE, Ncpus = 2); stopifnot(as.character(packageVersion("DESeq2")) == "1.46.0"); stopifnot(as.character(packageVersion("airway")) == "1.26.0")'
