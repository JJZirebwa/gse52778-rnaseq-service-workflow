#!/usr/bin/env bash
set -euo pipefail

repository_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
output_dir="${1:-${repository_root}/results}"
work_dir="${2:-${repository_root}/work}"
image_name="gse52778-rnaseq-service-workflow:1.0.0"
git_commit="$(git -C "${repository_root}" rev-parse HEAD)"

if [[ -e "${output_dir}" ]]; then
    echo "Refusing to mix outputs: ${output_dir} already exists." >&2
    exit 2
fi

mkdir -p "${output_dir}/nextflow"
docker build --tag "${image_name}" --file "${repository_root}/containers/Dockerfile" "${repository_root}"

nextflow run "${repository_root}/main.nf" \
    -profile docker \
    -work-dir "${work_dir}" \
    --samplesheet "${repository_root}/data/samples.csv" \
    --outdir "${output_dir}" \
    --container_image "${image_name}" \
    --git_commit "${git_commit}" \
    -with-report "${output_dir}/nextflow/report.html" \
    -with-timeline "${output_dir}/nextflow/timeline.html" \
    -with-trace "${output_dir}/nextflow/trace.txt" \
    -with-dag "${output_dir}/nextflow/dag.html"
