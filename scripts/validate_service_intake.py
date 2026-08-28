#!/usr/bin/env python3
"""Validate the fixed GSE52778 processed-count service intake contract."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path


REQUIRED_COLUMNS = {
    "sample_id",
    "cell",
    "dex",
    "dataset_accession",
    "input_type",
}


def validate_samplesheet(path: Path) -> dict:
    errors: list[str] = []
    try:
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            columns = set(reader.fieldnames or [])
            rows = list(reader)
    except (OSError, csv.Error) as exc:
        return _payload("invalid", 0, 0, [f"could not read sample sheet: {exc}"])

    missing = sorted(REQUIRED_COLUMNS - columns)
    if missing:
        errors.append(f"missing required columns: {', '.join(missing)}")

    if len(rows) != 8:
        errors.append(f"expected exactly 8 samples; found {len(rows)}")

    if REQUIRED_COLUMNS <= columns:
        sample_ids = [row["sample_id"].strip() for row in rows]
        if any(not sample_id for sample_id in sample_ids):
            errors.append("sample_id values must be non-empty")
        if len(sample_ids) != len(set(sample_ids)):
            errors.append("sample_id values must be unique")

        accessions = {row["dataset_accession"].strip() for row in rows}
        if accessions != {"GSE52778"}:
            errors.append("dataset_accession must be GSE52778 for every row")

        input_types = {row["input_type"].strip() for row in rows}
        if input_types != {"processed_counts"}:
            errors.append("input_type must be processed_counts for every row")

        pairs: dict[str, list[str]] = defaultdict(list)
        for row in rows:
            pairs[row["cell"].strip()].append(row["dex"].strip())
        if len(pairs) != 4:
            errors.append(f"expected exactly 4 cell pairs; found {len(pairs)}")
        for cell, treatments in sorted(pairs.items()):
            if Counter(treatments) != Counter({"untrt": 1, "trt": 1}):
                errors.append(f"cell {cell} must contain one untrt and one trt sample")
        pair_count = sum(
            Counter(treatments) == Counter({"untrt": 1, "trt": 1})
            for treatments in pairs.values()
        )
    else:
        pair_count = 0

    return _payload("valid" if not errors else "invalid", len(rows), pair_count, errors)


def _payload(status: str, sample_count: int, pair_count: int, errors: list[str]) -> dict:
    return {
        "status": status,
        "dataset_accession": "GSE52778",
        "input_type": "processed_counts",
        "sample_count": sample_count,
        "pair_count": pair_count,
        "design": "~ cell + dex",
        "reference_level": "untrt",
        "contrast": {
            "factor": "dex",
            "numerator": "trt",
            "denominator": "untrt",
        },
        "errors": errors,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samplesheet", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    payload = validate_samplesheet(args.samplesheet)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return 0 if payload["status"] == "valid" else 1


if __name__ == "__main__":
    raise SystemExit(main())

