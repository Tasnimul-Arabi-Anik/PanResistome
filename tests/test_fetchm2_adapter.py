import sys
import tempfile
import unittest
import subprocess
import os
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
            mlst_dir = sample_dir / "mlst" / "merged_output"
            metadata_dir.mkdir(parents=True)
            abricate_dir.mkdir()
            amrfinder_dir.mkdir(parents=True)
            mlst_dir.mkdir(parents=True)

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
            (mlst_dir / "mlst_merged.csv").write_text(
                "Assembly Accession,sample,scheme,st,allele_profile,sequence_type\n"
                "GCF_000000001.1,sampleA.fna,koxytoca,199,gapA(2);infB(2),koxytoca:ST199\n",
                encoding="utf-8",
            )
            (mlst_dir / "mlst_unknown.csv").write_text(
                "Assembly Accession,sample,scheme,st,allele_profile,sequence_type\n"
                "GCF_000000002.1,sampleB.fna,-,-,,-:ST-\n",
                encoding="utf-8",
            )

            outputs = export_contract(sample_dir, sample_dir / "panr2_inputs")

            self.assertTrue(Path(outputs["amr"]).exists())
            self.assertTrue(Path(outputs["amrfinderplus"]).exists())
            self.assertTrue(Path(outputs["mlst"]).exists())
            all_features = pd.read_csv(outputs["all_features"], sep="\t")
            self.assertEqual(set(all_features["database"]), {"amr", "amrfinderplus", "mlst"})
            self.assertEqual(
                set(all_features["assembly_accession"]),
                {"GCF_000000001.1", "GCF_000000002.1"},
            )
            self.assertTrue({"koxytoca:ST199", "ST_199", "gapA_2", "infB_2"}.issubset(set(all_features["feature_id"])))
            self.assertFalse({"-:ST-", "ST_-"}.intersection(set(all_features["feature_id"])))
            summary = Path(outputs["schema_validation_summary"]).read_text(encoding="utf-8")
            self.assertIn("feature_rows=6", summary)
            unmatched = pd.read_csv(outputs["unmatched_features"])
            self.assertEqual(len(unmatched), 0)
            audit = pd.read_csv(outputs["feature_completeness_audit"], sep="\t")
            mlst_status = audit.loc[audit["database"] == "mlst", "status"].iloc[0]
            self.assertEqual(mlst_status, "PASS")
            self.assertTrue(Path(outputs["all_feature_matrix"]).exists())
            self.assertTrue(Path(outputs["feature_cooccurrence"]).exists())

    def test_isfinder_blast_converter_writes_abricate_style_tables(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            blast = root / "sampleA.blast.tsv"
            results = root / "sampleA_results.tab"
            summary = root / "sampleA_summary.tab"
            blast.write_text(
                "contig1\tIS26\t99.5\t800\t1000\t820\t5\t804\t1\t800\t1e-100\t500\n"
                "contig2\tISlow\t85.0\t500\t1000\t900\t10\t509\t1\t500\t1e-20\t100\n",
                encoding="utf-8",
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "isfinder_blast_to_abricate.py"),
                    "--blast",
                    str(blast),
                    "--sample-id",
                    "sampleA",
                    "--out-results",
                    str(results),
                    "--out-summary",
                    str(summary),
                    "--min-identity",
                    "90",
                    "--min-coverage",
                    "80",
                ],
                check=True,
            )

            table = pd.read_csv(results, sep="\t")
            self.assertEqual(len(table), 1)
            self.assertEqual(table.loc[0, "GENE"], "IS26")
            self.assertEqual(table.loc[0, "DATABASE"], "isfinder")
            self.assertAlmostEqual(float(table.loc[0, "%COVERAGE"]), 97.56, places=2)
            summary_text = summary.read_text(encoding="utf-8")
            self.assertIn("IS26", summary_text)

    def test_cross_database_proximity_outputs_same_contig_evidence(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_dir = Path(tmpdir) / "Klebsiella_oxytoca"
            metadata_dir = sample_dir / "metadata_output"
            metadata_dir.mkdir(parents=True)
            pd.DataFrame(
                [
                    {
                        "Assembly Accession": "GCF_000000001.1",
                        "Organism Name": "Klebsiella oxytoca",
                        "Country": "Bangladesh",
                        "Isolation_Source_SD": "blood",
                    }
                ]
            ).to_csv(metadata_dir / "ncbi_clean.csv", index=False)

            result_header = "#FILE\tSEQUENCE\tSTART\tEND\tGENE\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n"
            for directory, filename, database, gene, start, end in [
                ("abricate", "ncbi_results.tab", "ncbi", "blaTEM-1", 100, 500),
                ("isfinder/tables", "isfinder_results.tab", "isfinder", "IS26", 700, 1200),
                ("plasmidfinder", "plasmidfinder_results.tab", "plasmidfinder", "IncFIB", 1400, 1800),
                ("integronfinder", "integronfinder_results.tab", "integronfinder", "complete_integron_intI1", 450, 900),
            ]:
                out_dir = sample_dir / directory
                out_dir.mkdir(parents=True)
                (out_dir / filename).write_text(
                    result_header
                    + f"GCF_000000001.1.fna\tcontigA\t{start}\t{end}\t{gene}\t100\t99\t{database}\tACC\tproduct\tcategory\n",
                    encoding="utf-8",
                )

            outputs = export_contract(sample_dir, sample_dir / "panr2_inputs")
            proximity = pd.read_csv(outputs["feature_proximity"], sep="\t")
            self.assertFalse(proximity.empty)
            self.assertIn("level_3_same_contig_within_10kb", set(proximity["interpretation_level"]))
            self.assertIn("level_4_same_contig_overlapping", set(proximity["interpretation_level"]))
            self.assertIn("evidence_level", proximity.columns)
            self.assertIn("interpretation_warning", proximity.columns)
            amr_mge = pd.read_csv(outputs["amr_mge_same_contig"], sep="\t")
            self.assertTrue({"isfinder", "integronfinder"}.issubset(set(amr_mge["feature_b_database"])))
            context = pd.read_csv(outputs["amr_mge_context"], sep="\t")
            self.assertEqual(context.loc[0, "same_contig_evidence"], "yes")
            self.assertTrue(Path(outputs["cross_database_interpretation_html"]).exists())
            self.assertTrue(Path(outputs["top_findings_html"]).exists())

    def test_bioproject_bias_and_amrfinder_abricate_concordance_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_dir = Path(tmpdir) / "Klebsiella_pneumoniae"
            metadata_dir = sample_dir / "metadata_output"
            abricate_dir = sample_dir / "abricate"
            amrfinder_dir = sample_dir / "amrfinderplus" / "raw"
            metadata_dir.mkdir(parents=True)
            abricate_dir.mkdir()
            amrfinder_dir.mkdir(parents=True)

            metadata_rows = []
            map_lines = ["sample_id,Assembly Accession\n"]
            for idx in range(1, 7):
                accession = f"GCF_00000000{idx}.1"
                country = "CountryA" if idx <= 3 else "CountryB"
                bioproject = "PRJNA_DOMINANT" if idx <= 3 else f"PRJNA_OTHER_{idx}"
                metadata_rows.append(
                    {
                        "Assembly Accession": accession,
                        "Organism Name": "Klebsiella pneumoniae",
                        "Country": country,
                        "Assembly BioProject Accession": bioproject,
                    }
                )
                map_lines.append(f"sample{idx},{accession}\n")
            pd.DataFrame(metadata_rows).to_csv(metadata_dir / "ncbi_clean.csv", index=False)
            (metadata_dir / "sample_map.csv").write_text("".join(map_lines), encoding="utf-8")

            header = "#FILE\tSEQUENCE\tSTART\tEND\tGENE\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n"
            abricate_lines = [header]
            for idx in range(1, 4):
                abricate_lines.append(
                    f"sample{idx}.fna\tcontig{idx}\t10\t90\tblaABC\t100\t99.5\tncbi\tACC{idx}\tbeta-lactamase\tbeta-lactam\n"
                )
            abricate_lines.append(
                "sample4.fna\tcontig4\t10\t90\ttetB\t100\t99.5\tncbi\tACC4\ttetracycline efflux\ttetracycline\n"
            )
            (abricate_dir / "ncbi_results.tab").write_text("".join(abricate_lines), encoding="utf-8")
            (amrfinder_dir / "calls.tsv").write_text(
                "sample_id\tGene symbol\tClass\t% Identity to reference sequence\t% Coverage of reference sequence\tContig id\tStart\tStop\n"
                "sample1\tblaABC\tbeta-lactam\t99.0\t100\tcontig1\t12\t88\n"
                "sample4\ttet(A)\ttetracycline\t98.0\t99\tcontig4\t15\t85\n",
                encoding="utf-8",
            )

            outputs = export_contract(sample_dir, sample_dir / "panr2_inputs")
            top_findings = pd.read_csv(outputs["top_findings"], sep="\t")
            self.assertIn("warning_flags", top_findings.columns)
            self.assertIn("interpretation_label", top_findings.columns)
            self.assertIn("single_bioproject_dominance", ";".join(top_findings["warning_flags"].fillna("")))
            self.assertIn("bioproject_bias_warning", set(top_findings["interpretation_label"]))

            bioproject = pd.read_csv(outputs["bioproject_bias_report"], sep="\t")
            self.assertIn("single_bioproject_dominance", ";".join(bioproject["warning"].fillna("")))
            self.assertTrue(Path(outputs["bioproject_bias_html"]).exists())

            concordance = pd.read_csv(outputs["amrfinder_abricate_concordance"], sep="\t")
            self.assertIn("called_by_both", set(concordance["status"]))
            self.assertIn("possible_class_match", set(concordance["status"]))
            self.assertTrue(Path(outputs["amrfinder_abricate_concordance_html"]).exists())

    def test_database_setup_status_passes_required_abricate_databases(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            sample_dir = root / "Klebsiella_oxytoca"
            (sample_dir / "metadata_output").mkdir(parents=True)
            (sample_dir / "sequence").mkdir()
            (sample_dir / "metadata_output" / "ncbi_clean.csv").write_text(
                "Assembly Accession,Organism Name\nGCF_000000001.1,Klebsiella oxytoca\n",
                encoding="utf-8",
            )
            (sample_dir / "sequence" / "GCF_000000001.1_genomic.fna").write_text(
                ">contig1\nATGC\n",
                encoding="utf-8",
            )
            fake_bin = root / "bin"
            fake_bin.mkdir()
            (fake_bin / "abricate").write_text(
                "#!/bin/sh\n"
                "if [ \"$1\" = \"--list\" ]; then\n"
                "  printf 'DATABASE\\tSEQUENCES\\n'\n"
                "  printf 'ncbi\\t10\\n'\n"
                "  printf 'vfdb\\t10\\n'\n"
                "  printf 'plasmidfinder\\t10\\n'\n"
                "else\n"
                "  printf 'abricate 1.0.1\\n'\n"
                "fi\n",
                encoding="utf-8",
            )
            (fake_bin / "panr").write_text("#!/bin/sh\nprintf 'panr 0.1.3\\n'\n", encoding="utf-8")
            for executable in fake_bin.iterdir():
                executable.chmod(0o755)
            out = root / "database_setup_status.tsv"
            env = os.environ.copy()
            env["PATH"] = f"{fake_bin}:{env.get('PATH', '')}"

            subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "database_setup_status.py"),
                    "--sample-dir",
                    str(sample_dir),
                    "--out",
                    str(out),
                    "--panr2-dbs",
                    "ncbi,vfdb,plasmidfinder",
                    "--run-panr2-comprehensive",
                    "true",
                    "--strict",
                ],
                check=True,
                env=env,
            )

            status = pd.read_csv(out, sep="\t")
            db_rows = status[status["database_or_tool"].str.startswith("abricate_db:")]
            self.assertEqual(set(db_rows["status"]), {"PASS"})
            isfinder = status.loc[status["database_or_tool"] == "isfinder_authorized_fasta", "status"].iloc[0]
            self.assertEqual(isfinder, "SKIPPED")

    def test_database_setup_status_fails_missing_required_abricate_database(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            sample_dir = root / "Klebsiella_oxytoca"
            (sample_dir / "metadata_output").mkdir(parents=True)
            (sample_dir / "sequence").mkdir()
            (sample_dir / "metadata_output" / "ncbi_clean.csv").write_text(
                "Assembly Accession,Organism Name\nGCF_000000001.1,Klebsiella oxytoca\n",
                encoding="utf-8",
            )
            (sample_dir / "sequence" / "GCF_000000001.1_genomic.fna").write_text(
                ">contig1\nATGC\n",
                encoding="utf-8",
            )
            fake_bin = root / "bin"
            fake_bin.mkdir()
            (fake_bin / "abricate").write_text(
                "#!/bin/sh\n"
                "if [ \"$1\" = \"--list\" ]; then\n"
                "  printf 'DATABASE\\tSEQUENCES\\n'\n"
                "  printf 'ncbi\\t10\\n'\n"
                "else\n"
                "  printf 'abricate 1.0.1\\n'\n"
                "fi\n",
                encoding="utf-8",
            )
            (fake_bin / "panr").write_text("#!/bin/sh\nprintf 'panr 0.1.3\\n'\n", encoding="utf-8")
            for executable in fake_bin.iterdir():
                executable.chmod(0o755)
            out = root / "database_setup_status.tsv"
            env = os.environ.copy()
            env["PATH"] = f"{fake_bin}:{env.get('PATH', '')}"

            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "database_setup_status.py"),
                    "--sample-dir",
                    str(sample_dir),
                    "--out",
                    str(out),
                    "--panr2-dbs",
                    "ncbi,vfdb",
                    "--run-panr2-comprehensive",
                    "true",
                    "--strict",
                ],
                check=False,
                env=env,
            )

            self.assertNotEqual(result.returncode, 0)
            status = pd.read_csv(out, sep="\t")
            vfdb_status = status.loc[status["database_or_tool"] == "abricate_db:vfdb", "status"].iloc[0]
            self.assertEqual(vfdb_status, "FAIL")

    def test_validation_summary_script_reports_manifest_metrics(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            sample_dir = root / "Klebsiella_oxytoca"
            manifest = sample_dir / "panr2_inputs" / "manifest"
            features = sample_dir / "panr2_inputs" / "features"
            metadata = sample_dir / "metadata_output"
            qc = sample_dir / "qc"
            report = sample_dir / "report"
            for directory in [manifest, features, metadata, qc, report, sample_dir / "sequence"]:
                directory.mkdir(parents=True)
            (metadata / "ncbi_clean.csv").write_text(
                "Assembly Accession,Organism Name\nGCF_000000001.1,Klebsiella oxytoca\n",
                encoding="utf-8",
            )
            (sample_dir / "sequence" / "GCF_000000001.1_genomic.fna").write_text(">c\nATGC\n", encoding="utf-8")
            (qc / "qc_master_report.csv").write_text(
                "assembly_accession,qc_master_status\nGCF_000000001.1,PASS\n",
                encoding="utf-8",
            )
            (manifest / "schema_validation_summary.txt").write_text(
                "feature_rows=1\nunmatched_feature_rows=0\ninvalid_feature_rows=0\nduplicate_feature_rows=0\n",
                encoding="utf-8",
            )
            (manifest / "database_setup_status.tsv").write_text(
                "database_or_tool\trequired_for_profile\tchecked\tstatus\tsetup_action\tversion_or_path\tmessage\n"
                "abricate_db:ncbi\ttrue\ttrue\tPASS\tpresent\t-\tOK\n",
                encoding="utf-8",
            )
            (manifest / "feature_completeness_audit.tsv").write_text(
                "database\texpected_from_profile\tmodule_enabled\traw_output_found\tfeature_table_found\tfeature_rows\tunique_features\tsamples_with_features\tsamples_processed\tstatus\tmessage\n"
                "amr\ttrue\ttrue\ttrue\ttrue\t1\t1\t1\t1\tPASS\tOK\n",
                encoding="utf-8",
            )
            (features / "amr.features.tsv").write_text(
                "sample_id\tassembly_accession\tdatabase\tfeature_id\tfeature_category\tpresence\tidentity\tcoverage\tcontig\tstart\tend\ttool\ttool_version\tdatabase_version\n"
                "GCF_000000001.1\tGCF_000000001.1\tamr\tblaTEM\tbeta-lactam\t1\t99\t100\tcontig1\t1\t100\tabricate\t1.0\tncbi\n",
                encoding="utf-8",
            )
            (features / "all_features.tsv").write_text((features / "amr.features.tsv").read_text(encoding="utf-8"), encoding="utf-8")
            (report / "index.html").write_text("<html></html>\n", encoding="utf-8")
            nested = sample_dir / "tool_results" / "integronfinder" / "raw" / "sampleA" / "panr2_inputs"
            nested.mkdir(parents=True)
            (nested / "not_a_handoff.txt").write_text("nested copy marker\n", encoding="utf-8")

            subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "summarize_validation_run.py"),
                    "--run-dir",
                    str(root),
                    "--out-dir",
                    str(root),
                ],
                check=True,
            )

            summary = pd.read_csv(root / "validation_summary.csv")
            metrics = set(summary["metric"])
            self.assertIn("schema_feature_rows", metrics)
            self.assertIn("database_setup_required_failures", metrics)
            self.assertIn("feature_table_rows", metrics)
            self.assertFalse(summary["sample_dir"].str.contains("tool_results").any())
            self.assertTrue((root / "validation_summary.md").exists())


if __name__ == "__main__":
    unittest.main()
