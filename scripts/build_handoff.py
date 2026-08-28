#!/usr/bin/env python3
"""Build a self-contained handoff directory and provenance manifest."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import re
import shutil
import subprocess


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def command_output(command: list[str]) -> str | None:
    try:
        return subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def package_version(session_info: str, package: str) -> str | None:
    match = re.search(rf"\b{re.escape(package)}_([0-9.]+)\b", session_info)
    return match.group(1) if match else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--git-commit", required=True)
    parser.add_argument("--container-image", required=True)
    parser.add_argument("--base-container", required=True)
    parser.add_argument("--files", nargs="+", required=True, type=Path)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=False)
    copied: list[Path] = []
    for source in args.files:
        destination = args.output_dir / source.name
        shutil.copy2(source, destination)
        copied.append(destination)

    session_path = args.output_dir / "sessionInfo.txt"
    session_info = session_path.read_text(encoding="utf-8")
    nextflow_version = command_output(["nextflow", "-version"])
    container_image_id = command_output(
        ["docker", "image", "inspect", "--format", "{{.Id}}", args.container_image]
    )

    manifest = {
        "dataset_accession": "GSE52778",
        "input_type": "processed_counts",
        "git_commit": args.git_commit,
        "analysis_container": {
            "reference": args.container_image,
            "local_image_id": container_image_id or "unavailable",
            "pinned_official_base": args.base_container,
        },
        "software": {
            "python": platform.python_version(),
            "nextflow": nextflow_version or "unavailable",
            "DESeq2": package_version(session_info, "DESeq2") or "unavailable",
            "airway": package_version(session_info, "airway") or "unavailable",
        },
        "analysis_contract": {
            "design": "~ cell + dex",
            "reference_level": "untrt",
            "contrast": "dex: trt versus untrt",
            "low_count_filter": "rowSums(counts >= 10) >= 4",
        },
        "sha256": {path.name: sha256(path) for path in sorted(copied)},
    }
    (args.output_dir / "provenance_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
