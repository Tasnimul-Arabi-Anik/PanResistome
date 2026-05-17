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
from check_genomad_readiness import resolve_database_dir
from check_container_readiness import as_list as container_as_list
from check_container_readiness import image_exec_command
from check_container_readiness import parse_args as parse_container_args
from check_comprehensive_validation_outputs import check_sample_dir


class FetchM2AdapterTests(unittest.TestCase):
    def test_resolves_nested_genomad_database_dir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "genomad_cache"
            nested = root / "genomad_db"
            nested.mkdir(parents=True)
            (nested / "marker.tsv").write_text("ok\n", encoding="utf-8")

            self.assertEqual(resolve_database_dir(root), nested)

    def test_container_readiness_parses_comma_separated_paths(self):
        self.assertEqual(
            container_as_list(" /db/checkm2 , /db/genomad ,, "),
            ["/db/checkm2", "/db/genomad"],
        )

    def test_container_readiness_builds_runtime_pull_test_commands(self):
        self.assertEqual(
            image_exec_command("singularity", "/usr/bin/singularity", "docker://alpine:3.19"),
            ["/usr/bin/singularity", "exec", "docker://alpine:3.19", "true"],
        )
        self.assertEqual(
            image_exec_command("docker", "/usr/bin/docker", "alpine:3.19"),
            ["/usr/bin/docker", "run", "--rm", "--entrypoint", "true", "alpine:3.19"],
        )

    def test_container_readiness_accepts_long_pull_test_timeout(self):
        args = parse_container_args(
            [
                "--runtime",
                "singularity",
                "--image",
                "docker://example/image:latest",
                "--pull-test",
                "--pull-test-timeout",
                "7200",
            ]
        )
        self.assertEqual(args.pull_test_timeout, 7200)

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
                    {"Assembly Accession": "GCF_000000001.1", "Organism Name": "Klebsiella oxytoca", "Collection_Year": 2019, "Country": "Bangladesh"},
                    {"Assembly Accession": "GCF_000000002.1", "Organism Name": "Klebsiella oxytoca", "Collection_Year": 2021, "Country": "United States"},
                    {"Assembly Accession": "GCF_000000003.1", "Organism Name": "Klebsiella oxytoca", "Collection_Year": 2023, "Country": "United Kingdom"},
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
            self.assertTrue(Path(outputs["report_controls"]).exists())
            self.assertTrue(Path(outputs["report_controls_html"]).exists())
            self.assertTrue((sample_dir / "basic" / "enriched_genome_dataset.csv").exists())
            self.assertTrue((sample_dir / "basic" / "enriched_genome_dataset.tsv").exists())
            self.assertTrue((sample_dir / "important" / "results.html").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "geographic_distribution_map.html").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "geographic_distribution_map.png").exists())
            self.assertTrue((sample_dir / "important" / "key_tables" / "geographic_distribution.tsv").exists())
            self.assertTrue((sample_dir / "important" / "key_tables" / "qc_step_summary.tsv").exists())
            self.assertTrue((sample_dir / "important" / "key_tables" / "qc_by_genome.tsv").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "qc_funnel.png").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "qc_status_overview.svg").exists())
            self.assertTrue((sample_dir / "important" / "key_tables" / "feature_prevalence_summary.tsv").exists())
            self.assertTrue((sample_dir / "important" / "key_tables" / "feature_variation_summary.tsv").exists())
            self.assertTrue((sample_dir / "important" / "key_tables" / "feature_variation_hits.tsv").exists())
            self.assertTrue((sample_dir / "important" / "key_tables" / "temporal_database_burden.tsv").exists())
            self.assertTrue((sample_dir / "important" / "key_tables" / "temporal_feature_prevalence.tsv").exists())
            self.assertTrue((sample_dir / "important" / "key_tables" / "temporal_trend_summary.tsv").exists())
            self.assertTrue((sample_dir / "important" / "key_tables" / "temporal_increasing_features.tsv").exists())
            self.assertTrue((sample_dir / "important" / "key_tables" / "temporal_decreasing_features.tsv").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "temporal_database_burden_top20.svg").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "temporal_feature_heatmap_top40.png").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "temporal_trends.html").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "temporal_selected_feature_prevalence.svg").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "temporal_slope_top40.png").exists())
            temporal_summary = pd.read_csv(sample_dir / "important" / "key_tables" / "temporal_trend_summary.tsv", sep="\t")
            self.assertTrue({"trend_label", "support_label", "temporal_pattern_label", "warning_flags"}.issubset(temporal_summary.columns))
            temporal_html = (sample_dir / "important" / "figures" / "temporal_trends.html").read_text(encoding="utf-8")
            for control in ["Database", "Trend", "Support", "Feature", "Selected Feature Prevalence", "First-to-Last Year Slope"]:
                self.assertIn(control, temporal_html)
            report_html = (sample_dir / "important" / "results.html").read_text(encoding="utf-8")
            for section in [
                "Featured Results",
                "Run Overview",
                "QC Summary",
                "Prevalence",
                "Geographic Distribution",
                "Variations",
                "Temporal Trends",
                "Warnings And Limitations",
                "Important Files",
            ]:
                self.assertIn(section, report_html)

    def test_optional_table_inputs_export_contract_features(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_dir = Path(tmpdir) / "Klebsiella_pneumoniae"
            metadata_dir = sample_dir / "metadata_output"
            metadata_dir.mkdir(parents=True)
            pd.DataFrame(
                [
                    {"Assembly Accession": "GCF_000000001.1", "Organism Name": "Klebsiella pneumoniae"},
                    {"Assembly Accession": "GCF_000000002.1", "Organism Name": "Klebsiella pneumoniae"},
                ]
            ).to_csv(metadata_dir / "ncbi_clean.csv", index=False)

            header = "#FILE\tSEQUENCE\tSTART\tEND\tGENE\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n"
            optional_tables = [
                ("mobileelementfinder", "mobileelementfinder_results.tab", "IS26", "mobile_element"),
                ("isfinder/tables", "isfinder_results.tab", "ISEcp1", "insertion_sequence"),
                ("mobsuite/tables", "mobsuite_results.tab", "IncFIB", "replicon"),
                ("prophage/tables", "prophage_results.tab", "region_1", "prophage"),
                ("defensefinder/tables", "defensefinder_results.tab", "RM_Type_I", "defense_system"),
                ("kleborate/tables", "kleborate_results.tab", "K_locus_KL1", "capsule_locus"),
                ("kaptive/tables", "kaptive_results.tab", "KL1", "locus_type"),
                ("ectyper/tables", "ectyper_results.tab", "O1:H7", "serotype"),
                ("serotypefinder/tables", "serotypefinder_results.tab", "O2", "serotype"),
                ("sccmecfinder/tables", "sccmecfinder_results.tab", "SCCmec_IV", "cassette_type"),
            ]
            for index, (directory, filename, feature_id, category) in enumerate(optional_tables, start=1):
                out_dir = sample_dir / directory
                out_dir.mkdir(parents=True)
                sample = "GCF_000000001.1" if index % 2 else "GCF_000000002.1"
                database = directory.split("/")[0]
                (out_dir / filename).write_text(
                    header
                    + f"{sample}.fna\tcontig{index}\t{index * 10}\t{index * 10 + 50}\t{feature_id}\t100\t99\t{database}\tACC{index}\t{category}\t{category}\n",
                    encoding="utf-8",
                )

            outputs = export_contract(sample_dir, sample_dir / "panr2_inputs")
            all_features = pd.read_csv(outputs["all_features"], sep="\t")
            expected_databases = {
                "mobileelementfinder",
                "isfinder",
                "mobsuite",
                "prophage",
                "defensefinder",
                "kleborate",
                "kaptive",
                "ectyper",
                "serotypefinder",
                "sccmecfinder",
            }
            self.assertEqual(set(all_features["database"]), expected_databases)
            for database in expected_databases:
                self.assertTrue(Path(outputs[database]).exists())
                tool_values = set(all_features.loc[all_features["database"] == database, "tool"])
                self.assertEqual(tool_values, {database})
            self.assertEqual(len(pd.read_csv(outputs["unmatched_features"])), 0)
            self.assertEqual(len(pd.read_csv(outputs["invalid_feature_rows"])), 0)
            audit = pd.read_csv(outputs["feature_completeness_audit"], sep="\t")
            present = audit[audit["database"].isin(expected_databases)]
            self.assertEqual(set(present["status"]), {"PASS"})

    def test_optional_runner_empty_tables_create_header_only_feature_tables(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_dir = Path(tmpdir) / "Klebsiella_pneumoniae"
            metadata_dir = sample_dir / "metadata_output"
            metadata_dir.mkdir(parents=True)
            pd.DataFrame(
                [{"Assembly Accession": "GCF_000000001.1", "Organism Name": "Klebsiella pneumoniae"}]
            ).to_csv(metadata_dir / "ncbi_clean.csv", index=False)

            empty_tables = [
                ("mobsuite/tables", "mobsuite.tsv"),
                ("prophage/tables", "prophage.tsv"),
                ("kleborate/tables", "kleborate.tsv"),
                ("kaptive/tables", "kaptive.tsv"),
                ("ectyper/tables", "ectyper.tsv"),
            ]
            for directory, filename in empty_tables:
                out_dir = sample_dir / directory
                out_dir.mkdir(parents=True)
                (out_dir / filename).write_text("sample_id\tstatus\n", encoding="utf-8")

            outputs = export_contract(sample_dir, sample_dir / "panr2_inputs")
            expected_databases = {"mobsuite", "prophage", "kleborate", "kaptive", "ectyper"}
            self.assertTrue(expected_databases.issubset(outputs.keys()))
            for database in expected_databases:
                table = pd.read_csv(outputs[database], sep="\t")
                self.assertEqual(len(table), 0)

            audit = pd.read_csv(outputs["feature_completeness_audit"], sep="\t")
            present = audit[audit["database"].isin(expected_databases)]
            self.assertEqual(set(present["status"]), {"WARNING_EMPTY"})
            self.assertEqual(set(present["feature_table_found"]), {True})

    def test_integronfinder_raw_output_without_features_creates_header_only_table(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_dir = Path(tmpdir) / "Klebsiella_pneumoniae"
            metadata_dir = sample_dir / "metadata_output"
            raw_dir = sample_dir / "tool_results" / "integronfinder" / "panr2_inputs"
            metadata_dir.mkdir(parents=True)
            raw_dir.mkdir(parents=True)
            pd.DataFrame(
                [{"Assembly Accession": "GCF_000000001.1", "Organism Name": "Klebsiella pneumoniae"}]
            ).to_csv(metadata_dir / "ncbi_clean.csv", index=False)
            (raw_dir / "integronfinder_summary.tab").write_text("sample_id\tstatus\n", encoding="utf-8")

            outputs = export_contract(sample_dir, sample_dir / "panr2_inputs")
            self.assertIn("integronfinder", outputs)
            table = pd.read_csv(outputs["integronfinder"], sep="\t")
            self.assertEqual(len(table), 0)

            audit = pd.read_csv(outputs["feature_completeness_audit"], sep="\t")
            status = audit.loc[audit["database"] == "integronfinder", "status"].iloc[0]
            self.assertEqual(status, "WARNING_EMPTY")

    def test_native_integronfinder_handoff_rows_are_exported(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_dir = Path(tmpdir) / "Klebsiella_pneumoniae"
            metadata_dir = sample_dir / "metadata_output"
            raw_dir = sample_dir / "tool_results" / "integronfinder" / "panr2_inputs"
            metadata_dir.mkdir(parents=True)
            raw_dir.mkdir(parents=True)
            pd.DataFrame(
                [{"Assembly Accession": "GCF_000000001.1", "Organism Name": "Klebsiella pneumoniae"}]
            ).to_csv(metadata_dir / "ncbi_clean.csv", index=False)
            (raw_dir / "integronfinder_results.tab").write_text(
                "#FILE\tSEQUENCE\tSTART\tEND\tGENE\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\n"
                "GCF_000000001.1_ASM1_genomic.fna\tcontigA\t100\t200\tcomplete\t100\t100\tintegronfinder\tcomplete\tcomplete\n",
                encoding="utf-8",
            )

            outputs = export_contract(sample_dir, sample_dir / "panr2_inputs")
            table = pd.read_csv(outputs["integronfinder"], sep="\t")
            self.assertEqual(len(table), 1)
            self.assertEqual(table.loc[0, "database"], "integronfinder")
            self.assertEqual(table.loc[0, "feature_id"], "complete")

            audit = pd.read_csv(outputs["feature_completeness_audit"], sep="\t")
            status = audit.loc[audit["database"] == "integronfinder", "status"].iloc[0]
            self.assertEqual(status, "PASS")

    def test_genomad_summary_collection_exports_clean_region_features(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_dir = Path(tmpdir) / "Klebsiella_pneumoniae"
            metadata_dir = sample_dir / "metadata_output"
            raw_dir = sample_dir / "prophage" / "raw" / "GCF_000000001.1_ASM1_genomic" / "sample_summary"
            tables_dir = sample_dir / "prophage" / "tables"
            metadata_dir.mkdir(parents=True)
            raw_dir.mkdir(parents=True)
            tables_dir.mkdir(parents=True)
            pd.DataFrame(
                [{"Assembly Accession": "GCF_000000001.1", "Organism Name": "Klebsiella pneumoniae"}]
            ).to_csv(metadata_dir / "ncbi_clean.csv", index=False)
            (raw_dir / "sample_virus_summary.tsv").write_text(
                "seq_name\tlength\ttopology\tcoordinates\tn_genes\tgenetic_code\tvirus_score\tfdr\tn_hallmarks\tmarker_enrichment\ttaxonomy\n"
                "contig1|provirus_10_100\t91\tProvirus\t10-100\t5\t11\t0.98\tNA\t2\t10\tViruses;Caudoviricetes\n",
                encoding="utf-8",
            )
            (raw_dir / "sample_plasmid_summary.tsv").write_text(
                "seq_name\tlength\ttopology\tn_genes\tgenetic_code\tplasmid_score\tfdr\tn_hallmarks\tmarker_enrichment\tconjugation_genes\tamr_genes\n"
                "contig2\t5000\tNo terminal repeats\t10\t11\t0.91\tNA\t1\t3\tNA\tNA\n",
                encoding="utf-8",
            )
            (raw_dir / "sample_genes.tsv").write_text(
                "gene\tstart\tend\nignored_gene\t1\t9\n",
                encoding="utf-8",
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "collect_optional_tool_tables.py"),
                    "--raw-dir",
                    str(sample_dir / "prophage" / "raw"),
                    "--out",
                    str(tables_dir / "prophage.tsv"),
                    "--tool",
                    "prophage",
                ],
                check=True,
            )
            collected = pd.read_csv(tables_dir / "prophage.tsv", sep="\t")
            self.assertEqual(set(collected["feature_id"]), {"viral_region:contig1|provirus_10_100", "plasmid_region:contig2"})
            self.assertEqual(set(collected["sample_id"]), {"GCF_000000001.1"})
            self.assertEqual(set(collected["category"]), {"viral_region", "plasmid_region"})

            outputs = export_contract(sample_dir, sample_dir / "panr2_inputs")
            prophage = pd.read_csv(outputs["prophage"], sep="\t")
            self.assertEqual(len(prophage), 2)
            self.assertEqual(set(prophage["feature_category"]), {"viral_region", "plasmid_region"})
            self.assertEqual(len(pd.read_csv(outputs["unmatched_features"])), 0)
            self.assertEqual(len(pd.read_csv(outputs["invalid_feature_rows"])), 0)

    def test_kleborate_real_output_exports_biological_features(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_dir = Path(tmpdir) / "Klebsiella_pneumoniae"
            metadata_dir = sample_dir / "metadata_output"
            kleborate_dir = sample_dir / "kleborate" / "tables"
            metadata_dir.mkdir(parents=True)
            kleborate_dir.mkdir(parents=True)
            pd.DataFrame(
                [{"Assembly Accession": "GCA_041085125.2", "Organism Name": "Klebsiella pneumoniae"}]
            ).to_csv(metadata_dir / "ncbi_clean.csv", index=False)
            (metadata_dir / "sample_map.csv").write_text(
                "sample_id,Assembly Accession\nGCA_041085125.2_ASM4108512v2_genomic,GCA_041085125.2\n",
                encoding="utf-8",
            )
            (kleborate_dir / "kleborate.tsv").write_text(
                "sample_id\ttool\tstrain\tklebsiella_pneumo_complex__mlst__ST\t"
                "klebsiella_pneumo_complex__virulence_score__virulence_score\t"
                "klebsiella_pneumo_complex__resistance_score__resistance_score\t"
                "klebsiella__ybst__Yersiniabactin\t"
                "klebsiella_pneumo_complex__kaptive__K_locus\t"
                "klebsiella_pneumo_complex__kaptive__O_type\t"
                "klebsiella_pneumo_complex__amr__Bla_chr\n"
                "GCA_041085125.2_ASM4108512v2_genomic\tkleborate\t"
                "GCA_041085125.2_ASM4108512v2_genomic\tST147\t1\t0\tybt 9; ICEKp?\tKL64\tO2α\tSHV-11^\n",
                encoding="utf-8",
            )

            outputs = export_contract(sample_dir, sample_dir / "panr2_inputs")
            features = pd.read_csv(outputs["kleborate"], sep="\t")
            feature_ids = set(features["feature_id"])
            self.assertTrue({"ST147", "virulence_score_1", "resistance_score_0", "KL64", "O2α", "SHV-11"}.issubset(feature_ids))
            self.assertIn("yersiniabactin_ybt_9", ";".join(feature_ids))
            self.assertEqual(set(features["tool"]), {"kleborate"})
            self.assertEqual(set(features["assembly_accession"]), {"GCA_041085125.2"})

    def test_mobsuite_realistic_output_exports_biological_features(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_dir = Path(tmpdir) / "Klebsiella_pneumoniae"
            metadata_dir = sample_dir / "metadata_output"
            mobsuite_dir = sample_dir / "mobsuite" / "tables"
            metadata_dir.mkdir(parents=True)
            mobsuite_dir.mkdir(parents=True)
            pd.DataFrame(
                [{"Assembly Accession": "GCA_000000001.1", "Organism Name": "Klebsiella pneumoniae"}]
            ).to_csv(metadata_dir / "ncbi_clean.csv", index=False)
            (mobsuite_dir / "mobsuite.tsv").write_text(
                "sample_id\tbiomarker\tqseqid\tpident\tqcovs\tsseqid\tsstart\tsend\tmolecule_type\trep_type(s)\tprimary_cluster_id\tsecondary_cluster_id\tmash_nearest_neighbor\tmge_type\tmge_subtype\n"
                "GCA_000000001.1_genomic\treplicon\t000100__NZ_CP016161_00012|IncFIB\t98.1\t100\tcontig1\t10\t500\tplasmid\tIncFIB\tAC125\tAL185\tCP021940\t\t\n"
                "GCA_000000001.1_genomic.fna:AC125\t\t\t\t\tcontig1\t600\t1800\tplasmid\t\t\t\t\tISKpn26\tIS5\n",
                encoding="utf-8",
            )

            outputs = export_contract(sample_dir, sample_dir / "panr2_inputs")
            features = pd.read_csv(outputs["mobsuite"], sep="\t")
            self.assertEqual(set(features["sample_id"]), {"GCA_000000001.1_genomic"})
            self.assertIn("IncFIB", set(features["feature_id"]))
            self.assertIn("AC125", set(features["feature_id"]))
            self.assertIn("ISKpn26", set(features["feature_id"]))
            self.assertIn("replicon", set(features["feature_category"]))
            self.assertIn("plasmid_cluster", set(features["feature_category"]))
            self.assertIn("is5", set(features["feature_category"]))

    def test_large_dataset_export_controls_limit_matrices_and_top_features(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_dir = Path(tmpdir) / "Klebsiella_pneumoniae"
            metadata_dir = sample_dir / "metadata_output"
            abricate_dir = sample_dir / "abricate"
            plasmid_dir = sample_dir / "plasmidfinder"
            vfdb_dir = sample_dir / "vfdb"
            metadata_dir.mkdir(parents=True)
            abricate_dir.mkdir()
            plasmid_dir.mkdir()
            vfdb_dir.mkdir()

            pd.DataFrame(
                [
                    {"Assembly Accession": "GCF_000000001.1", "Organism Name": "Klebsiella pneumoniae"},
                    {"Assembly Accession": "GCF_000000002.1", "Organism Name": "Klebsiella pneumoniae"},
                    {"Assembly Accession": "GCF_000000003.1", "Organism Name": "Klebsiella pneumoniae"},
                ]
            ).to_csv(metadata_dir / "ncbi_clean.csv", index=False)
            header = "#FILE\tSEQUENCE\tSTART\tEND\tGENE\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n"
            (abricate_dir / "ncbi_results.tab").write_text(
                header
                + "GCF_000000001.1.fna\tcontig1\t1\t100\tblaA\t100\t99\tncbi\tA1\tproduct\tbeta-lactam\n"
                + "GCF_000000002.1.fna\tcontig2\t1\t100\tblaA\t100\t99\tncbi\tA1\tproduct\tbeta-lactam\n"
                + "GCF_000000003.1.fna\tcontig3\t1\t100\tblaB\t100\t99\tncbi\tA2\tproduct\tbeta-lactam\n",
                encoding="utf-8",
            )
            (vfdb_dir / "vfdb_results.tab").write_text(
                header
                + "GCF_000000001.1.fna\tcontig1\t120\t220\tfimH\t100\t99\tvfdb\tV1\tadhesin\tvirulence\n"
                + "GCF_000000002.1.fna\tcontig2\t120\t220\tybtS\t100\t99\tvfdb\tV2\tsiderophore\tvirulence\n",
                encoding="utf-8",
            )
            (plasmid_dir / "plasmidfinder_results.tab").write_text(
                header
                + "GCF_000000001.1.fna\tcontig1\t150\t250\tIncFIB\t100\t99\tplasmidfinder\tP1\treplicon\tplasmid\n"
                + "GCF_000000002.1.fna\tcontig2\t150\t250\tIncFIB\t100\t99\tplasmidfinder\tP1\treplicon\tplasmid\n"
                + "GCF_000000003.1.fna\tcontig3\t150\t250\tIncX\t100\t99\tplasmidfinder\tP2\treplicon\tplasmid\n",
                encoding="utf-8",
            )

            outputs = export_contract(
                sample_dir,
                sample_dir / "panr2_inputs",
                large_dataset=True,
                report_mode="compact",
                max_features_heatmap=2,
                max_features_network=2,
                max_metadata_columns=5,
                top_n_features_per_database=1,
                skip_heavy_interactive_plots=True,
            )

            matrix = pd.read_csv(outputs["all_feature_matrix"], sep="\t")
            self.assertLessEqual(len(matrix.columns), 3)
            top_features = pd.read_csv(outputs["top_features_by_database"], sep="\t")
            self.assertLessEqual(top_features.groupby("database").size().max(), 1)
            controls = pd.read_csv(outputs["report_controls"], sep="\t")
            control_values = dict(zip(controls["setting"], controls["value"]))
            self.assertEqual(control_values["large_dataset"], "true")
            self.assertEqual(control_values["report_mode"], "compact")
            self.assertEqual(control_values["max_features_heatmap"], "2")
            self.assertEqual(control_values["skip_heavy_interactive_plots"], "true")
            proximity = pd.read_csv(outputs["feature_proximity"], sep="\t")
            proximity_features = set(zip(proximity["feature_a_database"], proximity["feature_a_id"]))
            proximity_features.update(zip(proximity["feature_b_database"], proximity["feature_b_id"]))
            self.assertLessEqual(len(proximity_features), 2)
            complete_proximity = pd.read_csv(outputs["feature_proximity_all"], sep="\t")
            complete_features = set(zip(complete_proximity["feature_a_database"], complete_proximity["feature_a_id"]))
            complete_features.update(zip(complete_proximity["feature_b_database"], complete_proximity["feature_b_id"]))
            self.assertGreater(len(complete_features), len(proximity_features))
            self.assertTrue(Path(outputs["report_controls_html"]).exists())

    def test_basic_output_mode_skips_important_user_bundle(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_dir = Path(tmpdir) / "Klebsiella_oxytoca"
            metadata_dir = sample_dir / "metadata_output"
            abricate_dir = sample_dir / "abricate"
            metadata_dir.mkdir(parents=True)
            abricate_dir.mkdir()
            pd.DataFrame(
                [
                    {
                        "Assembly Accession": "GCF_000000001.1",
                        "Organism Name": "Klebsiella oxytoca",
                        "Country": "Bangladesh",
                        "Collection_Year": "2024",
                        "combined_qc_status": "PASS",
                    }
                ]
            ).to_csv(metadata_dir / "ncbi_clean.csv", index=False)
            (abricate_dir / "ncbi_results.tab").write_text(
                "#FILE\tSEQUENCE\tSTART\tEND\tGENE\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n"
                "GCF_000000001.1.fna\tcontig1\t10\t90\tblaTEM-1\t100\t99.5\tncbi\tACC1\tbeta-lactamase\tbeta-lactam\n",
                encoding="utf-8",
            )

            export_contract(sample_dir, sample_dir / "panr2_inputs", output_mode="basic")

            self.assertTrue((sample_dir / "basic" / "enriched_genome_dataset.csv").exists())
            self.assertTrue((sample_dir / "basic" / "enriched_genome_dataset.tsv").exists())
            self.assertFalse((sample_dir / "important").exists())

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

    def test_lineage_diversity_and_statistical_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sample_dir = Path(tmpdir) / "Klebsiella_pneumoniae"
            metadata_dir = sample_dir / "metadata_output"
            abricate_dir = sample_dir / "abricate"
            mlst_dir = sample_dir / "tool_results" / "mlst" / "raw"
            ani_dir = sample_dir / "ani" / "analysis"
            metadata_dir.mkdir(parents=True)
            abricate_dir.mkdir()
            mlst_dir.mkdir(parents=True)
            ani_dir.mkdir(parents=True)

            metadata_rows = []
            sample_map = ["sample_id,Assembly Accession\n"]
            mlst_lines = []
            ani_lines = ["ani_cluster,representative,genome,cluster_size,duplicate_threshold\n"]
            abricate_lines = ["#FILE\tSEQUENCE\tSTART\tEND\tGENE\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n"]
            for idx in range(1, 7):
                accession = f"GCF_00000000{idx}.1"
                country = "CountryA" if idx <= 3 else "CountryB"
                st = "11" if idx <= 3 else f"{20 + idx}"
                ani_cluster = "ANI_CLUSTER_0001" if idx <= 3 else f"ANI_CLUSTER_000{idx}"
                metadata_rows.append(
                    {
                        "Assembly Accession": accession,
                        "Organism Name": "Klebsiella pneumoniae",
                        "Country": country,
                        "Assembly BioProject Accession": f"PRJNA_{idx}",
                    }
                )
                sample_map.append(f"sample{idx},{accession}\n")
                mlst_lines.append(f"sample{idx}.fna\tklebsiella\t{st}\tgapA({idx})\n")
                ani_lines.append(f"{ani_cluster}\t{accession}\t{accession}\t3\t99.9\n")
                abricate_lines.append(
                    f"sample{idx}.fna\tcontig{idx}\t1\t80\tcoreGene\t100\t99\tncbi\tCORE\tcore product\tcore\n"
                )
                if idx <= 3:
                    abricate_lines.append(
                        f"sample{idx}.fna\tcontig{idx}\t100\t180\tlineageGene\t100\t99\tncbi\tLIN\tlineage product\tlineage\n"
                    )
                if idx == 6:
                    abricate_lines.append(
                        "sample6.fna\tcontig6\t200\t260\trareGene\t100\t99\tncbi\tRARE\trare product\trare\n"
                    )
            pd.DataFrame(metadata_rows).to_csv(metadata_dir / "ncbi_clean.csv", index=False)
            (metadata_dir / "sample_map.csv").write_text("".join(sample_map), encoding="utf-8")
            (mlst_dir / "mlst.tsv").write_text("".join(mlst_lines), encoding="utf-8")
            (ani_dir / "duplicate_clusters.csv").write_text("".join(ani_lines), encoding="utf-8")
            (abricate_dir / "ncbi_results.tab").write_text("".join(abricate_lines), encoding="utf-8")

            outputs = export_contract(
                sample_dir,
                sample_dir / "panr2_inputs",
                rare_feature_threshold=0.2,
            )

            top_findings = pd.read_csv(outputs["top_findings"], sep="\t")
            self.assertIn("dominant_ST", top_findings.columns)
            self.assertIn("lineage_warning_flags", top_findings.columns)
            self.assertIn("single_ST_dominance", ";".join(top_findings["lineage_warning_flags"].fillna("")))

            lineage = pd.read_csv(outputs["lineage_summary"], sep="\t")
            self.assertIn("ST_11", set(lineage["mlst_ST"]))
            self.assertTrue(Path(outputs["lineage_context_html"]).exists())

            core = pd.read_csv(outputs["core_accessory_rare_features"], sep="\t")
            class_by_feature = dict(zip(core["feature_id"], core["feature_class"]))
            self.assertEqual(class_by_feature["coreGene"], "core")
            self.assertEqual(class_by_feature["lineageGene"], "accessory")
            self.assertEqual(class_by_feature["rareGene"], "rare")

            jaccard = pd.read_csv(outputs["jaccard_distance_matrix"], sep="\t")
            pair_ab = jaccard[(jaccard["sample_a"] == "GCF_000000001.1") & (jaccard["sample_b"] == "GCF_000000006.1")]["jaccard_distance"].iloc[0]
            pair_ba = jaccard[(jaccard["sample_a"] == "GCF_000000006.1") & (jaccard["sample_b"] == "GCF_000000001.1")]["jaccard_distance"].iloc[0]
            self.assertEqual(pair_ab, pair_ba)
            self.assertTrue(Path(outputs["diversity_summary_html"]).exists())

            stats = pd.read_csv(outputs["statistical_summary"], sep="\t")
            self.assertIn("lineage_warnings", set(stats["metric"]))
            self.assertTrue(Path(outputs["statistical_summary_html"]).exists())

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
            manifest = sample_dir / "panr2_inputs" / "manifest"
            manifest.mkdir(parents=True)
            (manifest / "abricate_database_setup_status.tsv").write_text(
                "database\trequested\tpresent_before\tsetup_requested\tupdate_requested\tsetup_status\tupdate_status\tpresent_after\tstatus\tmessage\n"
                "ncbi\ttrue\tfalse\ttrue\tfalse\tPASS\tSKIPPED\ttrue\tPASS\tok\n"
                "vfdb\ttrue\tfalse\ttrue\tfalse\tPASS\tSKIPPED\ttrue\tPASS\tok\n"
                "plasmidfinder\ttrue\tfalse\ttrue\tfalse\tPASS\tSKIPPED\ttrue\tPASS\tok\n",
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
            manifest = sample_dir / "panr2_inputs" / "manifest"
            manifest.mkdir(parents=True)
            (manifest / "abricate_database_setup_status.tsv").write_text(
                "database\trequested\tpresent_before\tsetup_requested\tupdate_requested\tsetup_status\tupdate_status\tpresent_after\tstatus\tmessage\n"
                "ncbi\ttrue\tfalse\ttrue\tfalse\tPASS\tSKIPPED\ttrue\tPASS\tok\n"
                "vfdb\ttrue\tfalse\ttrue\tfalse\tPASS\tSKIPPED\ttrue\tPASS\tok\n",
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

    def test_database_setup_status_marks_single_genome_pairwise_modules_inapplicable(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            sample_dir = root / "Klebsiella_oxytoca"
            (sample_dir / "metadata_output").mkdir(parents=True)
            (sample_dir / "ani" / "analysis").mkdir(parents=True)
            (sample_dir / "mash" / "analysis").mkdir(parents=True)
            (sample_dir / "metadata_output" / "ncbi_clean.csv").write_text(
                "Assembly Accession,Organism Name\nGCF_000000001.1,Klebsiella oxytoca\n",
                encoding="utf-8",
            )
            (sample_dir / "ani" / "analysis" / "panr2_ani_summary.csv").write_text(
                "sample_id,closest_genome,closest_ani,ani_cluster\n",
                encoding="utf-8",
            )
            (sample_dir / "ani" / "analysis" / "ani_run_status.tsv").write_text(
                "tool\tgenome_count\testimated_comparisons\tstrategy\tmax_all_vs_all_genomes\tlarge_dataset\tdecision\tstatus\tmessage\n"
                "fastani\t1\t1\tauto\t200\tfalse\tinsufficient_genomes\tSKIPPED_INAPPLICABLE\tANI requires at least two genomes; pairwise ANI was skipped.\n",
                encoding="utf-8",
            )
            (sample_dir / "mash" / "analysis" / "mash_distance_long.csv").write_text(
                "query,reference,mash_distance,p_value,matching_hashes\n",
                encoding="utf-8",
            )
            (sample_dir / "mash" / "analysis" / "mash_run_status.tsv").write_text(
                "tool\tgenome_count\tpair_rows\tnonself_pair_rows\tdecision\tstatus\tmessage\n"
                "mash\t1\t0\t0\tinsufficient_genomes\tSKIPPED_INAPPLICABLE\tMash requires at least two genomes; pairwise Mash screening was skipped.\n",
                encoding="utf-8",
            )
            out = root / "database_setup_status.tsv"

            subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "database_setup_status.py"),
                    "--sample-dir",
                    str(sample_dir),
                    "--out",
                    str(out),
                    "--run-ani",
                    "true",
                    "--run-mash",
                    "true",
                ],
                check=True,
            )

            status = pd.read_csv(out, sep="\t")
            ani_status = status.loc[status["database_or_tool"] == "ani", "status"].iloc[0]
            mash_status = status.loc[status["database_or_tool"] == "mash", "status"].iloc[0]
            self.assertEqual(ani_status, "SKIPPED_INAPPLICABLE")
            self.assertEqual(mash_status, "SKIPPED_INAPPLICABLE")

    def test_setup_abricate_databases_runs_optional_force_update(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            fake_bin = root / "bin"
            fake_bin.mkdir()
            log = root / "commands.log"
            (fake_bin / "panr").write_text(
                "#!/bin/sh\n"
                f"printf 'panr %s\\n' \"$*\" >> {log}\n"
                "printf 'setup complete\\n'\n",
                encoding="utf-8",
            )
            (fake_bin / "abricate-get_db").write_text(
                "#!/bin/sh\n"
                f"printf 'abricate-get_db %s\\n' \"$*\" >> {log}\n"
                "printf 'updated\\n'\n",
                encoding="utf-8",
            )
            (fake_bin / "abricate").write_text(
                "#!/bin/sh\n"
                "if [ \"$1\" = \"--list\" ]; then\n"
                "  printf 'DATABASE\\tSEQUENCES\\n'\n"
                "  printf 'ncbi\\t10\\n'\n"
                "  printf 'vfdb\\t10\\n'\n"
                "elif [ \"$1\" = \"--setupdb\" ]; then\n"
                f"  printf 'abricate --setupdb\\n' >> {log}\n"
                "  printf 'indexed\\n'\n"
                "else\n"
                "  printf 'abricate 1.0.1\\n'\n"
                "fi\n",
                encoding="utf-8",
            )
            for executable in fake_bin.iterdir():
                executable.chmod(0o755)

            out = root / "abricate_database_setup_status.tsv"
            env = os.environ.copy()
            env["PATH"] = f"{fake_bin}:{env.get('PATH', '')}"

            subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "setup_abricate_databases.py"),
                    "--dbs",
                    "ncbi,vfdb",
                    "--out",
                    str(out),
                    "--update",
                ],
                check=True,
                env=env,
            )

            status = pd.read_csv(out, sep="\t")
            self.assertEqual(set(status["status"]), {"PASS"})
            self.assertTrue(status["update_requested"].astype(str).str.lower().eq("true").all())
            commands = log.read_text(encoding="utf-8")
            self.assertIn("panr setup-db --dbs ncbi,vfdb", commands)
            self.assertIn("abricate-get_db --db ncbi --force", commands)
            self.assertIn("abricate-get_db --db vfdb --force", commands)
            self.assertIn("abricate --setupdb", commands)

    def test_setup_mobsuite_database_runs_mob_init_cache(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            fake_bin = root / "bin"
            fake_bin.mkdir()
            db_dir = root / "mobsuite_db"
            required_files = [
                "clusters.txt",
                "host_range_literature_plasmidDB.txt",
                "mob.proteins.faa",
                "mpf.proteins.faa",
                "ncbi_plasmid_full_seqs.fas",
                "ncbi_plasmid_full_seqs.fas.msh",
                "orit.fas",
                "rep.dna.fas",
                "repetitive.dna.fas",
                "ncbi_plasmid_full_seqs.fas.nhr",
                "repetitive.dna.fas.nhr",
            ]
            create_files = "\n".join([f"printf 'x\\n' > \"$db/{name}\"" for name in required_files])
            (fake_bin / "mob_init").write_text(
                "#!/bin/sh\n"
                "db=''\n"
                "while [ $# -gt 0 ]; do\n"
                "  case \"$1\" in\n"
                "    -d|--database_directory) shift; db=\"$1\" ;;\n"
                "  esac\n"
                "  shift\n"
                "done\n"
                "mkdir -p \"$db\"\n"
                f"{create_files}\n"
                "printf 'mob init ok\\n'\n",
                encoding="utf-8",
            )
            (fake_bin / "mob_init").chmod(0o755)
            out = root / "mobsuite_database_setup_status.tsv"
            env = os.environ.copy()
            env["PATH"] = f"{fake_bin}:{env.get('PATH', '')}"

            subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "setup_mobsuite_database.py"),
                    "--db-dir",
                    str(db_dir),
                    "--out",
                    str(out),
                    "--auto-init",
                    "true",
                    "--auto-init-taxa",
                    "false",
                ],
                check=True,
                env=env,
            )

            status = pd.read_csv(out, sep="\t").iloc[0]
            self.assertEqual(status["mob_init_status"], "PASS")
            self.assertEqual(status["core_status"], "PASS")
            self.assertEqual(status["status"], "WARNING_TAXA_MISSING")

    def test_setup_genomad_database_downloads_to_cache(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            fake_bin = root / "bin"
            fake_bin.mkdir()
            db_dir = root / "genomad_cache"
            (fake_bin / "genomad").write_text(
                "#!/bin/sh\n"
                "if [ \"$1\" = \"download-database\" ]; then\n"
                "  mkdir -p \"$2/genomad_db\"\n"
                "  printf 'db\\n' > \"$2/genomad_db/version.txt\"\n"
                "  printf 'downloaded\\n'\n"
                "else\n"
                "  printf 'genomad 1.0\\n'\n"
                "fi\n",
                encoding="utf-8",
            )
            (fake_bin / "genomad").chmod(0o755)
            out = root / "genomad_database_setup_status.tsv"
            env = os.environ.copy()
            env["PATH"] = f"{fake_bin}:{env.get('PATH', '')}"

            subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "setup_genomad_database.py"),
                    "--db-dir",
                    str(db_dir),
                    "--out",
                    str(out),
                    "--auto-download",
                    "true",
                ],
                check=True,
                env=env,
            )

            status = pd.read_csv(out, sep="\t").iloc[0]
            self.assertEqual(status["download_status"], "PASS")
            self.assertEqual(status["status"], "PASS")
            self.assertTrue(str(status["resolved_database_dir"]).endswith("genomad_db"))

    def test_database_setup_status_warns_when_mobsuite_taxa_sqlite_missing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            sample_dir = root / "Klebsiella_oxytoca"
            (sample_dir / "metadata_output").mkdir(parents=True)
            (sample_dir / "sequence").mkdir()
            (sample_dir / "mobsuite" / "tables").mkdir(parents=True)
            (sample_dir / "metadata_output" / "ncbi_clean.csv").write_text(
                "Assembly Accession,Organism Name\nGCF_000000001.1,Klebsiella oxytoca\n",
                encoding="utf-8",
            )
            (sample_dir / "sequence" / "GCF_000000001.1_genomic.fna").write_text(
                ">contig1\nATGC\n",
                encoding="utf-8",
            )
            mobsuite_db = root / "mobsuite_db"
            mobsuite_db.mkdir()
            for name in [
                "clusters.txt",
                "host_range_literature_plasmidDB.txt",
                "mob.proteins.faa",
                "mpf.proteins.faa",
                "ncbi_plasmid_full_seqs.fas",
                "ncbi_plasmid_full_seqs.fas.msh",
                "orit.fas",
                "rep.dna.fas",
                "repetitive.dna.fas",
                "ncbi_plasmid_full_seqs.fas.nhr",
                "repetitive.dna.fas.nhr",
            ]:
                (mobsuite_db / name).write_text("placeholder\n", encoding="utf-8")
            out = root / "database_setup_status.tsv"

            subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "database_setup_status.py"),
                    "--sample-dir",
                    str(sample_dir),
                    "--out",
                    str(out),
                    "--run-mobsuite",
                    "true",
                    "--mobsuite-db",
                    str(mobsuite_db),
                    "--strict",
                ],
                check=True,
            )

            status = pd.read_csv(out, sep="\t")
            row = status.loc[status["database_or_tool"] == "mobsuite_database"].iloc[0]
            self.assertEqual(row["status"], "WARNING")
            self.assertIn("taxa.sqlite", row["message"])

    def test_database_setup_status_does_not_strict_fail_experimental_optional_tools(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            sample_dir = root / "Acinetobacter_pittii"
            (sample_dir / "metadata_output").mkdir(parents=True)
            (sample_dir / "sequence").mkdir()
            (sample_dir / "mobsuite" / "tables").mkdir(parents=True)
            (sample_dir / "mobsuite" / "mobsuite_database_setup_status.tsv").write_text(
                "database_dir\tauto_init_requested\tauto_init_taxa_requested\tmob_init_status\ttaxa_init_status\tcore_status\ttaxa_status\tstatus\tmessage\n"
                f"{root / 'mobsuite_db'}\ttrue\ttrue\tFAIL\tSKIPPED\tFAIL\tWARNING_TAXA_MISSING\tFAIL\tcore database files missing\n",
                encoding="utf-8",
            )
            (sample_dir / "metadata_output" / "ncbi_clean.csv").write_text(
                "Assembly Accession,Organism Name\nGCF_000000001.1,Acinetobacter pittii\n",
                encoding="utf-8",
            )
            (sample_dir / "sequence" / "GCF_000000001.1_genomic.fna").write_text(
                ">contig1\nATGC\n",
                encoding="utf-8",
            )
            fake_bin = root / "bin"
            fake_bin.mkdir()
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
                    "--run-panr2-comprehensive",
                    "true",
                    "--run-defensefinder",
                    "true",
                    "--run-mobsuite",
                    "true",
                    "--mobsuite-db",
                    str(root / "mobsuite_db"),
                    "--strict",
                ],
                check=True,
                env=env,
            )

            status = pd.read_csv(out, sep="\t")
            defensefinder = status.loc[status["database_or_tool"] == "defensefinder"].iloc[0]
            mobsuite_setup = status.loc[status["database_or_tool"] == "mobsuite_database_setup"].iloc[0]
            mobsuite_db = status.loc[status["database_or_tool"] == "mobsuite_database"].iloc[0]
            self.assertEqual(defensefinder["status"], "WARNING_MISSING")
            self.assertEqual(str(defensefinder["required_for_profile"]).lower(), "false")
            self.assertEqual(mobsuite_setup["status"], "WARNING_FAILED")
            self.assertEqual(mobsuite_db["status"], "WARNING_FAILED")

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

    def test_comprehensive_validation_checker_accepts_nonfatal_mobileelementfinder_warning(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            sample_dir = root / "Acinetobacter_pittii"
            manifest = sample_dir / "panr2_inputs" / "manifest"
            features = sample_dir / "panr2_inputs" / "features"
            for directory in [
                manifest,
                features,
                sample_dir / "report",
                sample_dir / "checkm2",
                sample_dir / "metadata_output",
                sample_dir / "sequence_qc",
                sample_dir / "prophage",
            ]:
                directory.mkdir(parents=True)

            (sample_dir / "report" / "index.html").write_text("<html></html>\n", encoding="utf-8")
            (features / "all_features.tsv").write_text(
                "sample_id\tassembly_accession\tdatabase\tfeature_id\n"
                "GCF_000000001.1\tGCF_000000001.1\tamr\tblaTEM\n",
                encoding="utf-8",
            )
            (features / "prophage.features.tsv").write_text(
                "sample_id\tassembly_accession\tdatabase\tfeature_id\n",
                encoding="utf-8",
            )
            (manifest / "schema_validation_summary.txt").write_text(
                "feature_rows=1\nunmatched_feature_rows=0\ninvalid_feature_rows=0\nduplicate_feature_rows=0\n",
                encoding="utf-8",
            )
            (manifest / "database_setup_status.tsv").write_text(
                "database_or_tool\trequired_for_profile\tstatus\tmessage\n"
                "checkm2\ttrue\tPASS\tOK\n"
                "mobileelementfinder\tfalse\tWARNING_MISSING\toptional\n",
                encoding="utf-8",
            )
            (manifest / "native_runner_merge_audit.tsv").write_text(
                "module\tstatus\tmessage\n"
                "abricate\tPASS\tOK\n"
                "integronfinder\tPASS\tOK\n"
                "mlst\tPASS\tOK\n"
                "mobileelementfinder\tWARNING_FAILED\toptional runner failed nonfatally\n",
                encoding="utf-8",
            )
            (sample_dir / "checkm2" / "quality_report.tsv").write_text(
                "Name\tCompleteness\tContamination\nGCF_000000001.1\t99.0\t0.1\n",
                encoding="utf-8",
            )
            (sample_dir / "sequence_qc" / "qc_decisions.tsv").write_text(
                "sequence_file\tcombined_qc_status\nGCF_000000001.1_genomic.fna\tPASS\n",
                encoding="utf-8",
            )
            (sample_dir / "metadata_output" / "ncbi_enriched.csv").write_text(
                "Assembly Accession,checkm2_completeness\nGCF_000000001.1,99.0\n",
                encoding="utf-8",
            )
            (sample_dir / "prophage" / "module_status.tsv").write_text(
                "module\tstatus\tmessage\ngenomad\tPASS\tOK\n",
                encoding="utf-8",
            )

            failures = check_sample_dir(
                sample_dir,
                require_checkm2=True,
                require_genomad=True,
                allow_mobileelementfinder_warning=True,
                expect_zero_schema_errors=True,
            )

            self.assertEqual(failures, [])


if __name__ == "__main__":
    unittest.main()
