import importlib.util
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "build_handoff.py"
SPEC = importlib.util.spec_from_file_location("build_handoff", SCRIPT)
BUILD_HANDOFF = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BUILD_HANDOFF)


class HandoffProvenanceTests(unittest.TestCase):
    def test_extracts_loaded_package_version_from_session_info(self):
        """Catches regex escaping that reports loaded packages as unavailable."""
        session_info = "otherPkgs:\n[1] airway_1.26.0 DESeq2_1.46.0"

        self.assertEqual(BUILD_HANDOFF.package_version(session_info, "DESeq2"), "1.46.0")
        self.assertEqual(BUILD_HANDOFF.package_version(session_info, "airway"), "1.26.0")


if __name__ == "__main__":
    unittest.main()

