import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "verify_results.py"
FIXTURES = ROOT / "tests" / "fixtures"
VALID_FULL = FIXTURES / "valid_full_results.csv"
VALID_FDR = FIXTURES / "valid_significant_results.csv"
VALID_METADATA = FIXTURES / "valid_analysis_metadata.json"


class ResultContractValidatorTests(unittest.TestCase):
    def run_validator(self, full=VALID_FULL, fdr=VALID_FDR, metadata=VALID_METADATA):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "result_validation.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--full-table",
                    str(full),
                    "--fdr-table",
                    str(fdr),
                    "--metadata",
                    str(metadata),
                    "--output",
                    str(output),
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            payload = json.loads(output.read_text()) if output.exists() else None
        return completed, payload

    def rewrite_csv(self, source, transform):
        tmpdir = tempfile.TemporaryDirectory()
        path = Path(tmpdir.name) / source.name
        with source.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
            fields = list(rows[0].keys())
        fields, rows = transform(fields, rows)
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        return tmpdir, path

    def rewrite_metadata(self, transform):
        tmpdir = tempfile.TemporaryDirectory()
        path = Path(tmpdir.name) / "analysis_metadata.json"
        payload = json.loads(VALID_METADATA.read_text())
        transform(payload)
        path.write_text(json.dumps(payload))
        return tmpdir, path

    def test_accepts_tables_and_explicit_analysis_contract(self):
        """Catches validators that reject compliant outputs or omit the sanity check."""
        completed, payload = self.run_validator()

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(payload["status"], "valid")
        self.assertEqual(payload["full_result_rows"], 3)
        self.assertEqual(payload["fdr_result_rows"], 1)
        self.assertEqual(payload["crispld2_log2_fold_change"], 2.25)
        self.assertEqual(payload["errors"], [])

    def test_rejects_missing_required_result_column(self):
        """Catches incomplete result exports that omit inferential fields."""
        def remove_padj(fields, rows):
            fields.remove("padj")
            for row in rows:
                row.pop("padj")
            return fields, rows

        tmpdir, full = self.rewrite_csv(VALID_FULL, remove_padj)
        self.addCleanup(tmpdir.cleanup)

        completed, payload = self.run_validator(full=full)

        self.assertEqual(completed.returncode, 1)
        self.assertIn("full table missing required columns: padj", payload["errors"])

    def test_rejects_fdr_value_above_threshold(self):
        """Catches rows incorrectly included in the FDR-filtered table."""
        def raise_padj(fields, rows):
            rows[0]["padj"] = "0.051"
            return fields, rows

        tmpdir, fdr = self.rewrite_csv(VALID_FDR, raise_padj)
        self.addCleanup(tmpdir.cleanup)

        completed, payload = self.run_validator(fdr=fdr)

        self.assertEqual(completed.returncode, 1)
        self.assertIn("FDR table padj values must be numeric and <= 0.05", payload["errors"])

    def test_rejects_nonpositive_crispld2_effect(self):
        """Catches a reversed or scientifically implausible stated contrast."""
        def reverse_effect(fields, rows):
            rows[0]["log2FoldChange"] = "-2.25"
            return fields, rows

        tmpdir, full = self.rewrite_csv(VALID_FULL, reverse_effect)
        self.addCleanup(tmpdir.cleanup)

        completed, payload = self.run_validator(full=full)

        self.assertEqual(completed.returncode, 1)
        self.assertIn(
            "CRISPLD2 must have positive log2FoldChange for trt versus untrt",
            payload["errors"],
        )

    def test_rejects_wrong_design_or_contrast(self):
        """Catches analysis metadata that drops pairing or reverses the contrast."""
        def mutate(payload):
            payload["design"] = "~ dex"
            payload["contrast"]["numerator"] = "untrt"
            payload["contrast"]["denominator"] = "trt"

        tmpdir, metadata = self.rewrite_metadata(mutate)
        self.addCleanup(tmpdir.cleanup)

        completed, payload = self.run_validator(metadata=metadata)

        self.assertEqual(completed.returncode, 1)
        self.assertIn("design must be ~ cell + dex", payload["errors"])
        self.assertIn("contrast must be dex: trt versus untrt", payload["errors"])

    def test_rejects_unapplied_or_wrong_low_count_filter(self):
        """Catches an implicit, skipped, or substituted filtering rule."""
        def mutate(payload):
            payload["low_count_filter"] = {
                "rule": "rowMeans(counts) > 1",
                "applied": False,
                "genes_before": 63677,
                "genes_after": 22724,
            }

        tmpdir, metadata = self.rewrite_metadata(mutate)
        self.addCleanup(tmpdir.cleanup)

        completed, payload = self.run_validator(metadata=metadata)

        self.assertEqual(completed.returncode, 1)
        self.assertIn(
            "low-count filter must be applied as rowSums(counts >= 10) >= 4",
            payload["errors"],
        )


if __name__ == "__main__":
    unittest.main()

