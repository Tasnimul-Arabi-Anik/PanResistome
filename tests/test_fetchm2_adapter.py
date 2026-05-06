import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from normalize_fetchm2_output import normalize_fetchm2_output
from export_panr2_inputs import parse_version_line
from panr2_contract import export_contract


class FetchM2AdapterTests(unittest.TestCase):
    def test_parses_tool_named_version_lines(self):
        self.assertEqual(
            parse_version_line("seqkit v2.13.0", "fetchm_env_versions"),
            {"component": "seqkit", "version": "2.13.0"},
        )
        self.assertEqual(
            parse_version_line("abricate 1.4.0", "abricate_env_versions"),
            {"component": "abricate", "version": "1.4.0"},
        )
        self.assertEqual(
            parse_version_line("QUAST v5.3.0", "quast_env_versions"),
            {"component": "QUAST", "version": "5.3.0"},
        )

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
            (metadata_dir / "fetchm2_all_assemblies.csv").write_text(
                "Assembly Accession,Assembly Name\nGCF_000001.1,ASM1\nGCA_000001.1,ASM1\n",
                encoding="utf-8",
            )
            (root / "metadata_analysis" / "tables").mkdir(parents=True)
            (root / "audit").mkdir()
            (root / "sequence").mkdir()

            sample_dir = normalize_fetchm2_output(root)

            self.assertEqual(sample_dir.name, "Klebsiella_oxytoca")
            self.assertTrue((sample_dir / "metadata_output" / "fetchm2_clean.csv").exists())
            self.assertTrue((sample_dir / "metadata_output" / "fetchm2_all_assemblies.csv").exists())
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

    def test_exports_and_validates_panr2_contract_features(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_dir = Path(tmpdir) / "Klebsiella_oxytoca"
            metadata_dir = sample_dir / "metadata_output"
            abricate_dir = sample_dir / "abricate"
            amrfinder_dir = sample_dir / "amrfinderplus" / "raw"
            metadata_dir.mkdir(parents=True)
            abricate_dir.mkdir()
            amrfinder_dir.mkdir(parents=True)

            pd.DataFrame(
                [
                    {"Assembly Accession": "GCF_000000001.1", "Organism Name": "Klebsiella oxytoca"},
                    {"Assembly Accession": "GCF_000000002.1", "Organism Name": "Klebsiella oxytoca"},
                ]
            ).to_csv(metadata_dir / "ncbi_clean.csv", index=False)
            (metadata_dir / "sample_map.csv").write_text(
                "sample_id,Assembly Accession\nsampleA,GCF_000000001.1\nsampleB,GCF_000000002.1\n",
                encoding="utf-8",
            )
            (abricate_dir / "ncbi_results.tab").write_text(
                "#FILE\tSEQUENCE\tSTART\tEND\tGENE\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n"
                "sampleA.fna\tcontig1\t10\t90\tblaTEM-1\t100\t99.5\tncbi\tACC1\tbeta-lactamase\tbeta-lactam\n",
                encoding="utf-8",
            )
            (amrfinder_dir / "sampleB.tsv").write_text(
                "sample_id\tGene symbol\tClass\t% Identity to reference sequence\t% Coverage of reference sequence\tContig id\tStart\tStop\n"
                "sampleB\ttet(A)\ttetracycline\t98.2\t99.0\tcontig2\t5\t80\n",
                encoding="utf-8",
            )

            outputs = export_contract(sample_dir, sample_dir / "panr2_inputs")

            self.assertTrue(Path(outputs["amr"]).exists())
            self.assertTrue(Path(outputs["amrfinderplus"]).exists())
            all_features = pd.read_csv(outputs["all_features"], sep="\t")
            self.assertEqual(set(all_features["database"]), {"amr", "amrfinderplus"})
            self.assertEqual(
                set(all_features["assembly_accession"]),
                {"GCF_000000001.1", "GCF_000000002.1"},
            )
            summary = Path(outputs["schema_validation_summary"]).read_text(encoding="utf-8")
            self.assertIn("feature_rows=2", summary)
            unmatched = pd.read_csv(outputs["unmatched_features"])
            self.assertEqual(len(unmatched), 0)


if __name__ == "__main__":
    unittest.main()
