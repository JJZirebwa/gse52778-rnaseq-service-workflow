#!/usr/bin/env python3
"""Validate differential-expression outputs and the fixed analysis contract."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Callable


REQUIRED_COLUMNS = {
    "gene_id",
    "gene_symbol",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
}
EXPECTED_CONTRAST = {
    "factor": "dex",
    "numerator": "trt",
    "denominator": "untrt",
}
EXPECTED_FILTER = "rowSums(counts >= 10) >= 4"


def _read_csv(path: Path, label: str, errors: list[str]) -> list[dict[str, str]]:
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            columns = set(reader.fieldnames or [])
            rows = list(reader)
    except (OSError, csv.Error) as exc:
        errors.append(f"could not read {label}: {exc}")
        return []

    missing = sorted(REQUIRED_COLUMNS - columns)
    if missing:
        errors.append(f"{label} missing required columns: {', '.join(missing)}")
    return rows


def _number(value: str | None) -> float | None:
    if value is None or value.strip() in {"", "NA", "NaN"}:
        return None
    try:
        parsed = float(value)
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) else None


def _all_numeric(rows: list[dict[str, str]], column: str, predicate: Callable[[float], bool]) -> bool:
    values = [_number(row.get(column)) for row in rows]
    return bool(values) and all(value is not None and predicate(value) for value in values)


def validate_results(full_path: Path, fdr_path: Path, metadata_path: Path) -> dict:
    errors: list[str] = []
    full_rows = _read_csv(full_path, "full table", errors)
    fdr_rows = _read_csv(fdr_path, "FDR table", errors)

    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        errors.append(f"could not read analysis metadata: {exc}")
        metadata = {}

    if full_rows and not _all_numeric(full_rows, "baseMean", lambda value: value >= 0):
        errors.append("full table baseMean values must be non-negative numeric values")
    if fdr_rows and not _all_numeric(fdr_rows, "padj", lambda value: 0 <= value <= 0.05):
        errors.append("FDR table padj values must be numeric and <= 0.05")

    full_gene_ids = {row.get("gene_id") for row in full_rows}
    if any(row.get("gene_id") not in full_gene_ids for row in fdr_rows):
        errors.append("every FDR table row must also occur in the full table")

    crispld2_values = [
        _number(row.get("log2FoldChange"))
        for row in full_rows
        if row.get("gene_symbol") == "CRISPLD2"
    ]
    crispld2_lfc = next((value for value in crispld2_values if value is not None), None)
    if crispld2_lfc is None or crispld2_lfc <= 0:
        errors.append("CRISPLD2 must have positive log2FoldChange for trt versus untrt")

    if metadata.get("dataset_accession") != "GSE52778":
        errors.append("dataset_accession must be GSE52778")
    if metadata.get("design") != "~ cell + dex":
        errors.append("design must be ~ cell + dex")
    if metadata.get("reference_level") != "untrt":
        errors.append("reference_level must be untrt")
    if metadata.get("contrast") != EXPECTED_CONTRAST:
        errors.append("contrast must be dex: trt versus untrt")

    low_count_filter = metadata.get("low_count_filter", {})
    if (
        low_count_filter.get("rule") != EXPECTED_FILTER
        or low_count_filter.get("applied") is not True
        or not isinstance(low_count_filter.get("genes_before"), int)
        or not isinstance(low_count_filter.get("genes_after"), int)
        or low_count_filter.get("genes_after", 0) <= 0
        or low_count_filter.get("genes_after", 0) >= low_count_filter.get("genes_before", 0)
    ):
        errors.append(f"low-count filter must be applied as {EXPECTED_FILTER}")

    if metadata.get("alpha") != 0.05:
        errors.append("alpha must be 0.05")

    return {
        "status": "valid" if not errors else "invalid",
        "dataset_accession": metadata.get("dataset_accession"),
        "design": metadata.get("design"),
        "reference_level": metadata.get("reference_level"),
        "contrast": metadata.get("contrast"),
        "low_count_filter": low_count_filter,
        "full_result_rows": len(full_rows),
        "fdr_result_rows": len(fdr_rows),
        "crispld2_log2_fold_change": crispld2_lfc,
        "errors": errors,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--full-table", required=True, type=Path)
    parser.add_argument("--fdr-table", required=True, type=Path)
    parser.add_argument("--metadata", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    payload = validate_results(args.full_table, args.fdr_table, args.metadata)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return 0 if payload["status"] == "valid" else 1


if __name__ == "__main__":
    raise SystemExit(main())

