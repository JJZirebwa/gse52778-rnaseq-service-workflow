import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "validate_service_intake.py"
VALID_SAMPLES = ROOT / "tests" / "fixtures" / "valid_samples.csv"


class ServiceIntakeValidatorTests(unittest.TestCase):
    def run_validator(self, samplesheet: Path):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "validation.json"
            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--samplesheet",
                    str(samplesheet),
                    "--output",
                    str(output),
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            payload = json.loads(output.read_text()) if output.exists() else None
        return completed, payload

    def write_rows(self, rows):
        tmpdir = tempfile.TemporaryDirectory()
        path = Path(tmpdir.name) / "samples.csv"
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
            writer.writeheader()
            writer.writerows(rows)
        return tmpdir, path

    def load_valid_rows(self):
        with VALID_SAMPLES.open(newline="") as handle:
            return list(csv.DictReader(handle))

    def test_accepts_exactly_eight_samples_in_four_complete_pairs(self):
        """Catches validators that reject the approved intake or omit contract metadata."""
        completed, payload = self.run_validator(VALID_SAMPLES)

        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertEqual(
            payload,
            {
                "status": "valid",
                "dataset_accession": "GSE52778",
                "input_type": "processed_counts",
                "sample_count": 8,
                "pair_count": 4,
                "design": "~ cell + dex",
                "reference_level": "untrt",
                "contrast": {
                    "factor": "dex",
                    "numerator": "trt",
                    "denominator": "untrt",
                },
                "errors": [],
            },
        )

    def test_rejects_missing_sample(self):
        """Catches a missing exact-eight-sample gate."""
        rows = self.load_valid_rows()[:-1]
        tmpdir, path = self.write_rows(rows)
        self.addCleanup(tmpdir.cleanup)

        completed, payload = self.run_validator(path)

        self.assertEqual(completed.returncode, 1)
        self.assertIn("expected exactly 8 samples; found 7", payload["errors"])

    def test_rejects_duplicate_sample_identifier(self):
        """Catches intake that could silently map two rows to one library."""
        rows = self.load_valid_rows()
        rows[-1]["sample_id"] = rows[0]["sample_id"]
        tmpdir, path = self.write_rows(rows)
        self.addCleanup(tmpdir.cleanup)

        completed, payload = self.run_validator(path)

        self.assertEqual(completed.returncode, 1)
        self.assertIn("sample_id values must be unique", payload["errors"])

    def test_rejects_incomplete_treatment_pair(self):
        """Catches a pair with two treated samples and no untreated control."""
        rows = self.load_valid_rows()
        rows[0]["dex"] = "trt"
        tmpdir, path = self.write_rows(rows)
        self.addCleanup(tmpdir.cleanup)

        completed, payload = self.run_validator(path)

        self.assertEqual(completed.returncode, 1)
        self.assertIn(
            "cell N61311 must contain one untrt and one trt sample",
            payload["errors"],
        )

    def test_rejects_wrong_accession_or_input_type(self):
        """Catches accidental substitution of another dataset or raw-read claims."""
        rows = self.load_valid_rows()
        rows[0]["dataset_accession"] = "GSE00000"
        rows[1]["input_type"] = "fastq"
        tmpdir, path = self.write_rows(rows)
        self.addCleanup(tmpdir.cleanup)

        completed, payload = self.run_validator(path)

        self.assertEqual(completed.returncode, 1)
        self.assertIn("dataset_accession must be GSE52778 for every row", payload["errors"])
        self.assertIn("input_type must be processed_counts for every row", payload["errors"])


if __name__ == "__main__":
    unittest.main()

