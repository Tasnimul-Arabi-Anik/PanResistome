import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from normalize_fetchm2_output import normalize_fetchm2_output


class FetchM2AdapterTests(unittest.TestCase):
    def test_normalizes_fetchm2_output_to_panresistome_layout(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "fetchm_results"
            metadata_dir = root / "metadata_output"
            metadata_dir.mkdir(parents=True)
            pd.DataFrame(
                [
                    {
                        "Assembly Accession": "GCF_000001.1",
                        "Assembly BioSample Accession": "SAMN000001",
                        "Organism Name": "Klebsiella oxytoca",
                        "Organism Taxonomic ID": 571,
                        "Country": "Bangladesh",
                        "Collection_Year": 2024,
                        "Host_SD": "Homo sapiens",
                        "Isolation_Source_SD": "blood",
                    }
                ]
            ).to_csv(metadata_dir / "fetchm2_clean.csv", index=False)
            (metadata_dir / "fetchm2_clean.tsv").write_text(
                "Assembly Accession\tCountry\nGCF_000001.1\tBangladesh\n",
                encoding="utf-8",
            )
            (root / "metadata_analysis" / "tables").mkdir(parents=True)
            (root / "audit").mkdir()
            (root / "sequence").mkdir()

            sample_dir = normalize_fetchm2_output(root)

            self.assertEqual(sample_dir.name, "Klebsiella_oxytoca")
            self.assertTrue((sample_dir / "metadata_output" / "fetchm2_clean.csv").exists())
            self.assertTrue((sample_dir / "metadata_output" / "ncbi_clean.csv").exists())
            self.assertTrue((sample_dir / "metadata_analysis" / "tables").is_dir())
            self.assertTrue((sample_dir / "sequence").is_dir())
            compat = pd.read_csv(sample_dir / "metadata_output" / "ncbi_clean.csv")
            self.assertEqual(compat.loc[0, "Geographic Location"], "Bangladesh")
            self.assertEqual(str(compat.loc[0, "Collection Date"]), "2024")
            self.assertEqual(compat.loc[0, "Host"], "Homo sapiens")
            self.assertEqual(compat.loc[0, "Isolation Source"], "blood")
            self.assertEqual(compat.loc[0, "Genus"], "Klebsiella")
            self.assertEqual(compat.loc[0, "Species"], "Klebsiella oxytoca")
            self.assertEqual(compat.loc[0, "TaxID"], 571)


if __name__ == "__main__":
    unittest.main()
