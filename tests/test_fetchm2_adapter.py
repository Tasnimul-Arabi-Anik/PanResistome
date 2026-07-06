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
from panr2_contract import (
    REPORT_FIGURE_REGISTRY,
    export_contract,
    _figure_caption_for,
    _figure_visibility_metadata,
    _figure_render_quality,
    _human_figure_title,
    _is_same_feature_pair,
    _kruskal_wallis,
    _svg_geographic_map,
    _write_cooccurrence_network_svg,
    _write_variation_scatter_svg,
    write_important_diversity_outputs,
    write_important_geographic_outputs,
    write_important_lineage_outputs,
    write_important_prevalence_outputs,
    write_important_temporal_outputs,
)
from check_genomad_readiness import resolve_database_dir
from check_container_readiness import as_list as container_as_list
from check_container_readiness import image_exec_command
from check_container_readiness import parse_args as parse_container_args
from check_comprehensive_validation_outputs import check_sample_dir
from check_important_report_outputs import (
    _check_highlight_quality,
    _check_prevalence_consistency,
    _check_temporal_placeholders,
    _check_visual_quality,
)


class FetchM2AdapterTests(unittest.TestCase):
    def test_important_geographic_distribution_counts_unique_genomes_and_flags_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            sample_dir = root / "Sample"
            metadata_dir = sample_dir / "metadata_output"
            features_dir = root / "panr2_inputs" / "features"
            metadata_dir.mkdir(parents=True)
            features_dir.mkdir(parents=True)
            pd.DataFrame(
                [
                    {"Assembly Accession": "G1", "Country": "Bangladesh", "Collection_Year": "2020", "Assembly BioProject Accession": "PRJ1"},
                    {"Assembly Accession": "G2", "Country": "Bangladesh", "Collection_Year": "2021", "Assembly BioProject Accession": "PRJ2"},
                    {"Assembly Accession": "G3", "Country": "", "Collection_Year": "2021", "Assembly BioProject Accession": "PRJ3"},
                ]
            ).to_csv(metadata_dir / "ncbi_clean.csv", index=False)
            pd.DataFrame(
                [
                    {"assembly_accession": "G1", "sample_id": "s1", "database": "amr", "feature_id": "blaA", "feature_name": "blaA", "feature_category": "beta_lactam", "presence": "1", "identity": "99", "coverage": "100", "contig": "c1", "start": "1", "end": "100", "tool": "abricate", "database_version": "v1"},
                    {"assembly_accession": "G1", "sample_id": "s1", "database": "amr", "feature_id": "blaA", "feature_name": "blaA", "feature_category": "beta_lactam", "presence": "1", "identity": "98", "coverage": "95", "contig": "c2", "start": "1", "end": "90", "tool": "abricate", "database_version": "v1"},
                    {"assembly_accession": "G2", "sample_id": "s2", "database": "vfdb", "feature_id": "fimH", "feature_name": "fimH", "feature_category": "adhesion", "presence": "1", "identity": "97", "coverage": "90", "contig": "c3", "start": "1", "end": "80", "tool": "abricate", "database_version": "v2"},
                ]
            ).to_csv(features_dir / "all_features.tsv", sep="\t", index=False)

            outputs = write_important_geographic_outputs(sample_dir, root / "panr2_inputs", sample_dir / "important", top_n=2)
            self.assertTrue(Path(outputs["important_geographic_database_burden"]).exists())
            self.assertTrue((sample_dir / "important" / "figures" / "geographic_distribution.html").exists())
            self.assertTrue((sample_dir / "important" / "geographic_tables.zip").exists())

            burden = pd.read_csv(sample_dir / "important" / "tables" / "geographic_database_burden.tsv", sep="\t")
            bd = burden[(burden["database"] == "amr") & (burden["geo_level"] == "country") & (burden["group_name"] == "Bangladesh")].iloc[0]
            self.assertEqual(int(bd["total_genomes"]), 2)
            self.assertEqual(int(bd["positive_genomes"]), 1)
            self.assertEqual(int(bd["total_feature_rows"]), 2)
            self.assertAlmostEqual(float(bd["prevalence_percent"]), 50.0)
            self.assertEqual(str(bd["continent"]), "Asia")
            self.assertEqual(str(bd["subcontinent"]), "South Asia")
            self.assertIn("small_group_warning", str(bd["warning_flags"]))

            features = pd.read_csv(sample_dir / "important" / "tables" / "geographic_feature_distribution.tsv", sep="\t")
            bla = features[(features["database"] == "amr") & (features["feature_id"] == "blaA") & (features["geo_level"] == "country") & (features["group_name"] == "Bangladesh")].iloc[0]
            self.assertEqual(int(bla["positive_genomes"]), 1)
            self.assertEqual(int(bla["feature_rows"]), 2)
            self.assertAlmostEqual(float(bla["mean_hits_per_positive_genome"]), 2.0)

            warnings = pd.read_csv(sample_dir / "important" / "tables" / "geographic_warning_summary.tsv", sep="\t")
            missing = warnings[(warnings["geo_level"] == "country") & (warnings["group_name"] == "missing")].iloc[0]
            self.assertIn("missing_country_metadata", str(missing["warning_flags"]))
            html = (sample_dir / "important" / "figures" / "geographic_distribution.html").read_text(encoding="utf-8")
            for control in ["Database", "Mode", "Feature", "Geographic level", "Minimum group size", "Warning filter"]:
                self.assertIn(control, html)
            self.assertIn("Gene map", html)
            self.assertIn("modeSelect.value = 'individual_feature'", html)
            self.assertIn("blaA", html)
            self.assertIn("Search feature / gene", html)
            self.assertIn("featureSearchInput", html)
            self.assertIn("Search detected features", html)
            self.assertIn("Feature / gene", html)
            self.assertIn("plotly.min.js", html)
            self.assertIn("Plotly.react", html)
            self.assertIn("type: 'choropleth'", html)
            self.assertIn("locationmode: 'country names'", html)
            self.assertIn("natural earth", html)
            self.assertIn("filled-country map", html)
            self.assertIn("map-workspace", html)
            self.assertIn("mapPlot", html)
            self.assertIn("overflow-x: hidden", html)
            self.assertIn("Map reading guide", html)
            self.assertIn("Higher selected-gene prevalence", html)
            self.assertIn("positive / total genomes", html)
            self.assertIn("small group", html)
            self.assertIn("Dataset-specific summary", html)
            self.assertNotIn("country tile", html.lower())
            self.assertNotIn("bubble", html.lower())
            self.assertTrue((sample_dir / "important" / "figures" / "plotly.min.js").exists())
            self.assertGreater((sample_dir / "important" / "figures" / "plotly.min.js").stat().st_size, 1_000_000)
            svg = (sample_dir / "important" / "figures" / "geographic_distribution_map.svg").read_text(encoding="utf-8")
            self.assertIn("Color shows dataset prevalence/burden", svg)
            self.assertIn("Map reading guide", svg)
            self.assertIn("positive / total genomes", svg)
            self.assertIn("higher selected-gene prevalence", svg)
            self.assertIn("not a regional or global prevalence estimate", svg)
            self.assertIn("small group", svg)
            self.assertNotIn("bubble", svg.lower())

    def test_important_prevalence_counts_unique_genomes_and_duplicate_rows(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            sample_dir = root / "Sample"
            basic_dir = sample_dir / "basic"
            features_dir = root / "panr2_inputs" / "features"
            basic_dir.mkdir(parents=True)
            features_dir.mkdir(parents=True)
            pd.DataFrame(
                [
                    {"assembly_accession": "G1", "sample_id": "s1"},
                    {"assembly_accession": "G2", "sample_id": "s2"},
                    {"assembly_accession": "G3", "sample_id": "s3"},
                ]
            ).to_csv(basic_dir / "enriched_genome_dataset.csv", index=False)
            pd.DataFrame(
                [
                    {"assembly_accession": "G1", "sample_id": "s1", "database": "amr", "feature_id": "blaA", "feature_name": "blaA", "feature_category": "beta_lactam", "feature_subcategory": "", "presence": "1", "identity": "99", "coverage": "100", "contig": "c1", "start": "1", "end": "100", "tool": "abricate", "database_version": "v1"},
                    {"assembly_accession": "G1", "sample_id": "s1", "database": "amr", "feature_id": "blaA", "feature_name": "blaA", "feature_category": "beta_lactam", "feature_subcategory": "", "presence": "1", "identity": "98", "coverage": "95", "contig": "c2", "start": "1", "end": "90", "tool": "abricate", "database_version": "v1"},
                    {"assembly_accession": "G2", "sample_id": "s2", "database": "amr", "feature_id": "blaA", "feature_name": "blaA", "feature_category": "beta_lactam", "feature_subcategory": "", "presence": "1", "identity": "97", "coverage": "90", "contig": "c1", "start": "1", "end": "80", "tool": "abricate", "database_version": "v1"},
                    {"assembly_accession": "G3", "sample_id": "s3", "database": "vfdb", "feature_id": "fimH", "feature_name": "fimH", "feature_category": "adhesion", "feature_subcategory": "", "presence": "1", "identity": "96", "coverage": "88", "contig": "", "start": "", "end": "", "tool": "abricate", "database_version": "v2"},
                ]
            ).to_csv(features_dir / "all_features.tsv", sep="\t", index=False)

            outputs = write_important_prevalence_outputs(sample_dir, root / "panr2_inputs", sample_dir / "important", top_n=2)
            self.assertTrue(Path(outputs["important_feature_prevalence"]).exists())
            self.assertTrue((sample_dir / "important" / "figures" / "prevalence_analysis.html").exists())

            prevalence = pd.read_csv(sample_dir / "important" / "tables" / "feature_prevalence.tsv", sep="\t")
            bla = prevalence.loc[prevalence["feature_id"] == "blaA"].iloc[0]
            self.assertEqual(int(bla["positive_genomes"]), 2)
            self.assertEqual(int(bla["feature_rows"]), 3)
            self.assertAlmostEqual(float(bla["prevalence_percent"]), 66.7, places=1)
            self.assertAlmostEqual(float(bla["mean_hits_per_positive_genome"]), 1.5, places=2)
            self.assertIn("duplicate_feature_rows_detected", str(bla["warning_flags"]))

            summary = pd.read_csv(sample_dir / "important" / "tables" / "prevalence_summary_by_database.tsv", sep="\t")
            self.assertTrue({"total_feature_rows", "unique_features", "median_features_per_genome", "top_feature_id"}.issubset(summary.columns))
            amr_summary = summary[summary["database"] == "amr"].iloc[0]
            self.assertEqual(int(amr_summary["positive_genomes"]), 2)
            self.assertGreaterEqual(int(amr_summary["positive_genomes"]), int(amr_summary["top_feature_positive_genomes"]))
            self.assertGreaterEqual(float(amr_summary["genomes_positive_percent"]), float(amr_summary["top_feature_prevalence_percent"]))
            html = (sample_dir / "important" / "figures" / "prevalence_analysis.html").read_text(encoding="utf-8")
            for control in ["Database", "Top 20", "Complete", "Genome prevalence %", "Positive genome count", "Feature row count", "Minimum prevalence %"]:
                self.assertIn(control, html)

    def test_temporal_outputs_drop_placeholder_years(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            sample_dir = root / "Sample"
            basic_dir = sample_dir / "basic"
            features_dir = root / "panr2_inputs" / "features"
            basic_dir.mkdir(parents=True)
            features_dir.mkdir(parents=True)
            pd.DataFrame(
                [
                    {"assembly_accession": "G1", "sample_id": "s1", "collection_year": "1-01-01", "bioproject": "PRJ1"},
                    {"assembly_accession": "G2", "sample_id": "s2", "collection_year": "1900", "bioproject": "PRJ1"},
                    {"assembly_accession": "G3", "sample_id": "s3", "collection_year": "2020", "bioproject": "PRJ2"},
                    {"assembly_accession": "G4", "sample_id": "s4", "collection_year": "2021", "bioproject": "PRJ3"},
                    {"assembly_accession": "G5", "sample_id": "s5", "collection_year": "2022", "bioproject": "PRJ4"},
                ]
            ).to_csv(basic_dir / "enriched_genome_dataset.csv", index=False)
            pd.DataFrame(
                [
                    {"assembly_accession": "G1", "sample_id": "s1", "database": "amr", "feature_id": "blaA", "presence": "1"},
                    {"assembly_accession": "G2", "sample_id": "s2", "database": "amr", "feature_id": "blaA", "presence": "1"},
                    {"assembly_accession": "G3", "sample_id": "s3", "database": "amr", "feature_id": "blaA", "presence": "1"},
                    {"assembly_accession": "G4", "sample_id": "s4", "database": "amr", "feature_id": "blaA", "presence": "1"},
                    {"assembly_accession": "G5", "sample_id": "s5", "database": "amr", "feature_id": "blaA", "presence": "1"},
                ]
            ).to_csv(features_dir / "all_features.tsv", sep="\t", index=False)

            write_important_temporal_outputs(sample_dir, root / "panr2_inputs", sample_dir / "important")
            burden = pd.read_csv(sample_dir / "important" / "key_tables" / "temporal_database_burden.tsv", sep="\t")
            self.assertNotIn(1900, set(burden["collection_year"].astype(int)))
            self.assertEqual(set(burden["collection_year"].astype(int)), {2020, 2021, 2022})

    def test_human_figure_title_conversions(self):
        self.assertEqual(_human_figure_title("geographic_map_vfdb_adeH"), "VFDB adeH geographic distribution")
        self.assertEqual(_human_figure_title("diversity_richness_by_metadata_isolation_source"), "Feature richness by Isolation Source")
        self.assertEqual(_human_figure_title("feature_profile_pcoa_by_bioproject"), "Feature-profile PCoA by BioProject")
        self.assertEqual(_human_figure_title("geographic_map_amr_ant_3_-IIa"), "AMR ant(3)-IIa geographic distribution")
        self.assertEqual(_human_figure_title("geographic_country_bar_amr_ant_3_-IIa"), "AMR ant(3)-IIa prevalence by country")
        self.assertEqual(_human_figure_title("variation_top_variable_amr_top20"), "AMR most variable features")
        self.assertEqual(_human_figure_title("variation_identity_coverage_vfdb_top20"), "VFDB identity vs coverage")

    def test_figure_caption_conversions_are_specific(self):
        caption_cases = {
            "variation_identity_coverage_amr_top20": "Identity-versus-coverage",
            "lineage_metadata_overlap_mlst_ST_isolation_source": "metadata groups are dominated by one lineage",
            "diversity_core_common_accessory_rare_by_database": "core, common, accessory, or rare",
            "temporal_slope_top40": "first and last valid-year prevalence",
            "metadata_volcano_amr_isolation_source_missing": "effect size and FDR support",
            "cooccurrence_heatmap_vfdb_vs_mlst": "feature-pair association direction",
            "prevalence_genomes_positive_by_database": "at least one feature from each database",
            "geographic_map_amr_burden": "Filled-country geographic map",
        }
        for stem, expected_phrase in caption_cases.items():
            with self.subTest(stem=stem):
                self.assertIn(expected_phrase, _figure_caption_for("", stem))
        self.assertTrue(_is_same_feature_pair("aph(3'')-Ib", "aph(3'')-Ib"))

    def test_geographic_map_preview_has_basemap_and_gene_map_guidance(self):
        svg = _svg_geographic_map(
            [
                {
                    "country": "Bangladesh",
                    "total_genomes": "10",
                    "positive_genomes": "7",
                    "prevalence_percent": "70.0",
                    "prevalence_display": "70.0% (7/10)",
                    "warning_flags": "",
                }
            ],
            "Geographic gene map",
        )

        self.assertIn("<polygon", svg)
        self.assertIn("Geographic gene map", svg)
        self.assertIn("higher selected-gene prevalence", svg)
        self.assertIn("Positive / total genomes: 7/10", svg)

    def test_empty_report_figures_render_explicit_unavailable_state(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            important_dir = Path(tmpdir) / "important"
            figures = important_dir / "figures"
            figures.mkdir(parents=True)

            scatter_stem = "variation_identity_coverage_amrfinderplus_top20"
            (figures / f"{scatter_stem}.data.tsv").write_text(
                "database\tfeature_id\tidentity\tcoverage\n"
                "amrfinderplus\tblaA\t\t\n",
                encoding="utf-8",
            )
            _write_variation_scatter_svg(figures / f"{scatter_stem}.svg", [{"database": "amrfinderplus", "feature_id": "blaA", "identity": "", "coverage": ""}], "AMRFinderPlus Identity vs Coverage")
            scatter_svg = (figures / f"{scatter_stem}.svg").read_text(encoding="utf-8")
            self.assertIn("No plottable data", scatter_svg)
            render_quality, axis_status, _action = _figure_render_quality(important_dir, scatter_stem)
            self.assertEqual(render_quality, "missing_required_numeric_metrics")
            self.assertEqual(axis_status, "axis_labels_present")

            network_stem = "cooccurrence_network_amr_vs_plasmidfinder"
            (figures / f"{network_stem}.data.tsv").write_text(
                "source_feature\ttarget_feature\tphi_correlation\n",
                encoding="utf-8",
            )
            _write_cooccurrence_network_svg(figures / f"{network_stem}.svg", [], [], "AMR vs PlasmidFinder Co-occurrence Network")
            network_svg = (figures / f"{network_stem}.svg").read_text(encoding="utf-8")
            self.assertIn("No plottable data", network_svg)
            render_quality, axis_status, _action = _figure_render_quality(important_dir, network_stem)
            self.assertEqual(render_quality, "empty_network")
            self.assertEqual(axis_status, "not_applicable")

    def test_figure_visibility_metadata_classifies_public_and_technical_figures(self):
        self.assertIn("notable_genomes_ranked", REPORT_FIGURE_REGISTRY)
        featured = _figure_visibility_metadata(
            "diversity_core_common_accessory_rare_by_database",
            "Diversity / Pan-feature Summary",
            "descriptive",
            "none",
            "asset_ready",
            "interpretation_ready",
            "human_readable_title",
            "specific_caption",
        )
        self.assertEqual(featured["default_visibility"], "featured")
        self.assertEqual(featured["publication_candidate"], "true")
        technical = _figure_visibility_metadata(
            "cooccurrence_network_amr_vs_vfdb",
            "Co-occurrence / Genomic Context",
            "exploratory",
            "none",
            "asset_ready",
            "exploratory_interpretation",
            "human_readable_title",
            "specific_caption",
        )
        self.assertEqual(technical["default_visibility"], "technical")
        self.assertEqual(technical["recommended_audience"], "technical users")

    def test_important_report_qa_catches_release_blockers(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            important_dir = Path(tmpdir) / "important"
            tables = important_dir / "tables"
            key_tables = important_dir / "key_tables"
            tables.mkdir(parents=True)
            key_tables.mkdir(parents=True)
            (tables / "prevalence_summary_by_database.tsv").write_text(
                "database\tpositive_genomes\tgenomes_positive_percent\ttop_feature_positive_genomes\ttop_feature_prevalence_percent\n"
                "amr\t1\t0.9\t50\t45.0\n",
                encoding="utf-8",
            )
            (tables / "report_highlights.tsv").write_text(
                "highlight_type\tprimary_feature\tsecondary_feature\n"
                "informative_cooccurrence\taph(3'')-Ib\taph(3'')-Ib\n",
                encoding="utf-8",
            )
            (tables / "report_highlights_by_section.tsv").write_text(
                "section\thighlight_type\n"
                + "".join(f"Co-occurrence / Genomic Context\trow{i}\n" for i in range(30)),
                encoding="utf-8",
            )
            (key_tables / "temporal_trend_summary.tsv").write_text(
                "database\tfeature_id\tfirst_year\tlast_year\n"
                "amr\tblaA\t1900\t2025\n",
                encoding="utf-8",
            )
            (tables / "report_visual_quality.tsv").write_text(
                "figure_stem\tsection\tasset_quality_label\trender_quality_label\taxis_label_status\tinterpretation_quality_label\ttitle_quality_label\tcaption_quality_label\tfinal_publication_label\tdefault_visibility\tpublication_candidate\n"
                "geographic_map_amr_burden\tGeographic Distribution\tasset_ready\trendered\tnot_applicable\tinterpretation_ready\thuman_readable_title\tspecific_caption\tpublication_ready\tstandard\ttrue\n"
                "variation_identity_coverage_amrfinderplus_top20\tVariations\tasset_ready\tmissing_required_numeric_metrics\tmissing_axis_labels\tinterpretation_ready\thuman_readable_title\tspecific_caption\tsupporting_only\tstandard\ttrue\n",
                encoding="utf-8",
            )

            errors: list[str] = []
            _check_prevalence_consistency(important_dir, errors)
            _check_highlight_quality(important_dir, errors)
            _check_temporal_placeholders(important_dir, errors)
            _check_visual_quality(important_dir, errors)

            combined = "\n".join(errors)
            self.assertIn("Database prevalence inconsistency", combined)
            self.assertIn("identical normalized features", combined)
            self.assertIn("only 1 distinct sections", combined)
            self.assertIn("placeholder first_year", combined)
            self.assertIn("Interpretation-sensitive figure marked publication_ready", combined)
            self.assertIn("render_quality_label=missing_required_numeric_metrics", combined)
            self.assertIn("missing axis labels", combined)

    def test_important_lineage_report_counts_and_warnings(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            sample_dir = root / "Sample"
            basic_dir = sample_dir / "basic"
            features_dir = root / "panr2_inputs" / "features"
            basic_dir.mkdir(parents=True)
            features_dir.mkdir(parents=True)
            pd.DataFrame(
                [
                    {"assembly_accession": f"G{i}", "sample_id": f"s{i}", "country": "CountryA" if i <= 3 else "CountryB", "isolation_source": "clinical" if i <= 3 else "environmental", "bioproject": f"PRJ{i}", "ani_cluster": "ANI_1" if i <= 3 else f"ANI_{i}"}
                    for i in range(1, 7)
                ]
            ).to_csv(basic_dir / "enriched_genome_dataset.csv", index=False)
            feature_rows = []
            for i in range(1, 7):
                st = "ST_11" if i <= 3 else f"ST_{20 + i}"
                feature_rows.append({"assembly_accession": f"G{i}", "sample_id": f"s{i}", "database": "mlst", "feature_id": st, "feature_name": st, "feature_category": "sequence_type", "presence": "1"})
                feature_rows.append({"assembly_accession": f"G{i}", "sample_id": f"s{i}", "database": "amr", "feature_id": "coreGene", "feature_name": "coreGene", "feature_category": "core", "presence": "1"})
                if i <= 3:
                    feature_rows.append({"assembly_accession": f"G{i}", "sample_id": f"s{i}", "database": "amr", "feature_id": "lineageGene", "feature_name": "lineageGene", "feature_category": "lineage", "presence": "1"})
            feature_rows.append({"assembly_accession": "G1", "sample_id": "s1", "database": "amr", "feature_id": "lineageGene", "feature_name": "lineageGene", "feature_category": "lineage", "presence": "1"})
            pd.DataFrame(feature_rows).to_csv(features_dir / "all_features.tsv", sep="\t", index=False)

            outputs = write_important_lineage_outputs(sample_dir, root / "panr2_inputs", sample_dir / "important", top_n=3)
            self.assertTrue(Path(outputs["important_lineage_summary"]).exists())
            self.assertTrue(Path(outputs["important_lineage_written_summaries"]).exists())
            self.assertTrue((sample_dir / "important" / "figures" / "lineage_clonal_structure.html").exists())
            self.assertTrue((sample_dir / "important" / "lineage_tables.zip").exists())

            distribution = pd.read_csv(sample_dir / "important" / "tables" / "lineage_distribution.tsv", sep="\t")
            st11 = distribution[(distribution["lineage_type"] == "mlst_ST") & (distribution["lineage_id"] == "ST_11")].iloc[0]
            self.assertEqual(int(st11["total_genomes"]), 3)
            self.assertAlmostEqual(float(st11["fraction_of_dataset"]), 0.5)

            overlap = pd.read_csv(sample_dir / "important" / "tables" / "lineage_metadata_overlap.tsv", sep="\t")
            country_a = overlap[(overlap["lineage_type"] == "mlst_ST") & (overlap["metadata_column"] == "country") & (overlap["metadata_group"] == "CountryA") & (overlap["lineage_id"] == "ST_11")].iloc[0]
            self.assertAlmostEqual(float(country_a["lineage_fraction_in_group"]), 1.0)
            self.assertEqual(country_a["dominance_label"], "severe_lineage_confounding")
            self.assertIn("metadata_lineage_confounding", str(country_a["warning_flags"]))

            presence = pd.read_csv(sample_dir / "important" / "tables" / "lineage_feature_presence.tsv", sep="\t")
            lin = presence[(presence["database"] == "amr") & (presence["feature_id"] == "lineageGene") & (presence["lineage_type"] == "mlst_ST") & (presence["lineage_id"] == "ST_11")].iloc[0]
            self.assertEqual(int(lin["positive_genomes"]), 3)
            self.assertEqual(int(lin["feature_rows"]), 4)
            self.assertAlmostEqual(float(lin["prevalence_percent"]), 100.0)

            enrichment = pd.read_csv(sample_dir / "important" / "tables" / "lineage_feature_enrichment.tsv", sep="\t")
            self.assertTrue({"odds_ratio", "p_value", "q_value", "support_label", "interpretation_label", "warning_flags"}.issubset(enrichment.columns))
            written = pd.read_csv(sample_dir / "important" / "tables" / "lineage_written_summaries.tsv", sep="\t")
            self.assertIn("overall_lineage_summary", set(written["section"]))
            html = (sample_dir / "important" / "figures" / "lineage_clonal_structure.html").read_text(encoding="utf-8")
            for control in ["Lineage type", "Database", "Feature mode", "Search feature", "Feature", "Metadata overlay", "Minimum lineage size", "custom", "Display", "Written Summaries"]:
                self.assertIn(control, html)

    def test_important_diversity_counts_unique_features_and_writes_complete_tables(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            sample_dir = root / "Sample"
            basic_dir = sample_dir / "basic"
            features_dir = root / "panr2_inputs" / "features"
            basic_dir.mkdir(parents=True)
            features_dir.mkdir(parents=True)
            pd.DataFrame(
                [
                    {"assembly_accession": "G1", "sample_id": "s1", "country": "Bangladesh", "isolation_source": "clinical", "bioproject": "PRJ1", "mlst_ST": "ST_11", "ani_cluster": "ANI_1"},
                    {"assembly_accession": "G2", "sample_id": "s2", "country": "Bangladesh", "isolation_source": "clinical", "bioproject": "PRJ1", "mlst_ST": "ST_11", "ani_cluster": "ANI_1"},
                    {"assembly_accession": "G3", "sample_id": "s3", "country": "India", "isolation_source": "environmental", "bioproject": "PRJ2", "mlst_ST": "ST_15", "ani_cluster": "ANI_2"},
                    {"assembly_accession": "G4", "sample_id": "s4", "country": "India", "isolation_source": "environmental", "bioproject": "PRJ2", "mlst_ST": "ST_15", "ani_cluster": "ANI_2"},
                ]
            ).to_csv(basic_dir / "enriched_genome_dataset.csv", index=False)
            pd.DataFrame(
                [
                    {"assembly_accession": "G1", "sample_id": "s1", "database": "amr", "feature_id": "coreA", "feature_name": "coreA", "feature_category": "beta_lactam", "presence": "1", "identity": "99", "coverage": "100"},
                    {"assembly_accession": "G2", "sample_id": "s2", "database": "amr", "feature_id": "coreA", "feature_name": "coreA", "feature_category": "beta_lactam", "presence": "1", "identity": "98", "coverage": "99"},
                    {"assembly_accession": "G3", "sample_id": "s3", "database": "amr", "feature_id": "coreA", "feature_name": "coreA", "feature_category": "beta_lactam", "presence": "1", "identity": "97", "coverage": "98"},
                    {"assembly_accession": "G4", "sample_id": "s4", "database": "amr", "feature_id": "coreA", "feature_name": "coreA", "feature_category": "beta_lactam", "presence": "1", "identity": "96", "coverage": "97"},
                    {"assembly_accession": "G1", "sample_id": "s1", "database": "vfdb", "feature_id": "rareV", "feature_name": "rareV", "feature_category": "adhesion", "presence": "1", "identity": "95", "coverage": "88"},
                    {"assembly_accession": "G1", "sample_id": "s1", "database": "vfdb", "feature_id": "rareV", "feature_name": "rareV", "feature_category": "adhesion", "presence": "1", "identity": "94", "coverage": "87"},
                    {"assembly_accession": "G2", "sample_id": "s2", "database": "plasmidfinder", "feature_id": "IncF", "feature_name": "IncF", "feature_category": "replicon", "presence": "1", "identity": "100", "coverage": "100"},
                ]
            ).to_csv(features_dir / "all_features.tsv", sep="\t", index=False)

            outputs = write_important_diversity_outputs(sample_dir, root / "panr2_inputs", sample_dir / "important", top_n=2)
            self.assertTrue(Path(outputs["important_diversity_feature_richness"]).exists())
            self.assertTrue((sample_dir / "important" / "figures" / "diversity_analysis.html").exists())
            self.assertTrue((sample_dir / "important" / "diversity_tables.zip").exists())
            self.assertTrue((sample_dir / "important" / "diversity_figures.zip").exists())

            richness = pd.read_csv(sample_dir / "important" / "tables" / "diversity_feature_richness_by_sample.tsv", sep="\t")
            g1 = richness[richness["assembly_accession"] == "G1"].iloc[0]
            self.assertEqual(int(g1["total_unique_features"]), 2)
            self.assertEqual(int(g1["total_feature_rows"]), 3)
            self.assertEqual(int(g1["vfdb_richness"]), 1)

            classes = pd.read_csv(sample_dir / "important" / "tables" / "diversity_core_common_accessory_rare_features.tsv", sep="\t")
            self.assertEqual(classes.loc[classes["feature_id"] == "coreA", "feature_class"].iloc[0], "core")
            self.assertEqual(int(classes.loc[classes["feature_id"] == "rareV", "feature_rows"].iloc[0]), 2)

            jaccard = pd.read_csv(sample_dir / "important" / "tables" / "diversity_jaccard_distance_matrix.tsv", sep="\t")
            diag = jaccard[(jaccard["sample_a"] == "G1") & (jaccard["sample_b"] == "G1")].iloc[0]
            self.assertAlmostEqual(float(diag["jaccard_distance"]), 0.0)
            self.assertTrue((jaccard["jaccard_distance"].astype(float) >= 0).all())

            accumulation = pd.read_csv(sample_dir / "important" / "tables" / "diversity_pan_feature_accumulation.tsv", sep="\t")
            self.assertTrue((accumulation["new_features_added"].astype(int) >= 0).all())
            self.assertTrue(accumulation["cumulative_unique_features"].is_monotonic_increasing)

            html = (sample_dir / "important" / "figures" / "diversity_analysis.html").read_text(encoding="utf-8")
            for control in ["Diversity scope", "Diversity view", "Metadata color/group", "Display", "Sort", "Pan-feature accumulation", "Jaccard similarity/distance"]:
                self.assertIn(control, html)
            for metric in ["total_feature_rows", "max_features_in_one_genome", "databases_represented", "jaccard_matrix_available", "pan_feature_curve_available"]:
                self.assertIn(metric, html)
            report_summary = pd.read_csv(sample_dir / "important" / "tables" / "diversity_report_summary.tsv", sep="\t")
            self.assertTrue({"total_feature_rows", "max_features_in_one_genome", "databases_represented", "jaccard_matrix_available", "pan_feature_curve_available"}.issubset(set(report_summary["metric"])))

    def test_kruskal_wallis_detects_multi_group_burden_difference(self):
        result = _kruskal_wallis([[1, 1, 2, 2], [5, 5, 6, 6], [9, 9, 10, 10]])
        self.assertIsNotNone(result)
        statistic, p_value = result
        self.assertGreater(statistic, 6.0)
        self.assertLess(p_value, 0.05)

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
            self.assertTrue((sample_dir / "important" / "figures" / "geographic_distribution.html").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "geographic_distribution_map.png").exists())
            self.assertTrue((sample_dir / "important" / "key_tables" / "geographic_distribution.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "geographic_distribution_summary.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "geographic_feature_distribution.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "geographic_database_burden.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "geographic_warning_summary.tsv").exists())
            self.assertTrue((sample_dir / "important" / "geographic_tables.zip").exists())
            self.assertTrue((sample_dir / "important" / "geographic_figures.zip").exists())
            self.assertTrue(any((sample_dir / "important" / "figures").glob("geographic_country_bar_*.data.tsv")))
            self.assertTrue(any((sample_dir / "important" / "figures").glob("geographic_map_*.svg")))
            self.assertTrue((sample_dir / "important" / "key_tables" / "qc_step_summary.tsv").exists())
            self.assertTrue((sample_dir / "important" / "key_tables" / "qc_by_genome.tsv").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "qc_funnel.png").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "qc_status_overview.svg").exists())
            self.assertTrue((sample_dir / "important" / "key_tables" / "feature_prevalence_summary.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "feature_prevalence.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "feature_prevalence_top.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "prevalence_summary_by_database.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "prevalence_core_accessory_rare_summary.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "prevalence_database_burden_by_sample.tsv").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "prevalence_analysis.html").exists())
            self.assertTrue((sample_dir / "important" / "prevalence_tables.zip").exists())
            self.assertTrue((sample_dir / "important" / "prevalence_figures.zip").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "prevalence_feature_counts_by_database.svg").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "prevalence_genomes_positive_by_database.png").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "prevalence_core_accessory_rare_by_database.data.tsv").exists())
            self.assertTrue(any((sample_dir / "important" / "figures").glob("prevalence_top_features_*.pdf")))
            self.assertTrue((sample_dir / "important" / "key_tables" / "feature_variation_summary.tsv").exists())
            self.assertTrue((sample_dir / "important" / "key_tables" / "feature_variation_hits.tsv").exists())
            self.assertTrue((sample_dir / "important" / "key_tables" / "feature_variation_database_summary.tsv").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "variation_analysis.html").exists())
            self.assertTrue((sample_dir / "important" / "variation_figures.zip").exists())
            self.assertTrue(any((sample_dir / "important" / "figures").glob("variation_identity_*_top20.pdf")))
            self.assertTrue(any((sample_dir / "important" / "figures").glob("variation_coverage_*_top20.svg")))
            self.assertTrue(any((sample_dir / "important" / "figures").glob("variation_identity_coverage_*_top20.png")))
            self.assertTrue(any((sample_dir / "important" / "figures").glob("variation_top_variable_*_top20.data.tsv")))
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
            self.assertTrue((sample_dir / "important" / "tables" / "cooccurrence_pair_summary.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "cooccurrence_heatmap_matrix.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "cooccurrence_network_edges.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "cooccurrence_network_nodes.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "genomic_context_evidence.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "contig_neighborhoods.tsv").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "cooccurrence_context.html").exists())
            self.assertTrue((sample_dir / "important" / "cooccurrence_tables.zip").exists())
            self.assertTrue((sample_dir / "important" / "cooccurrence_figures.zip").exists())
            self.assertTrue(any((sample_dir / "important" / "figures").glob("cooccurrence_heatmap_*_vs_*.svg")))
            self.assertTrue(any((sample_dir / "important" / "figures").glob("cooccurrence_heatmap_*_vs_*.pdf")))
            self.assertTrue(any((sample_dir / "important" / "figures").glob("cooccurrence_network_*_vs_*.png")))
            self.assertTrue((sample_dir / "important" / "tables" / "metadata_feature_enrichment.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "metadata_burden_associations.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "metadata_category_enrichment.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "metadata_association_summary.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "metadata_usability_summary.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "metadata_burden_omnibus.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "metadata_category_omnibus.tsv").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "metadata_associations.html").exists())
            self.assertTrue((sample_dir / "important" / "metadata_association_tables.zip").exists())
            self.assertTrue((sample_dir / "important" / "metadata_association_figures.zip").exists())
            self.assertTrue(any((sample_dir / "important" / "figures").glob("metadata_volcano_*_*.svg")))
            self.assertTrue(any((sample_dir / "important" / "figures").glob("metadata_enrichment_heatmap_*_*.pdf")))
            self.assertTrue(any((sample_dir / "important" / "figures").glob("metadata_burden_boxplot_*_*.png")))
            self.assertTrue((sample_dir / "important" / "tables" / "lineage_summary.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "lineage_distribution.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "lineage_metadata_overlap.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "lineage_feature_burden.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "lineage_feature_enrichment.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "lineage_adjusted_top_findings.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "lineage_feature_presence.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "lineage_written_summaries.tsv").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "lineage_clonal_structure.html").exists())
            self.assertTrue((sample_dir / "important" / "lineage_tables.zip").exists())
            self.assertTrue((sample_dir / "important" / "lineage_figures.zip").exists())
            self.assertTrue(any((sample_dir / "important" / "figures").glob("lineage_distribution_*.svg")))
            self.assertTrue(any((sample_dir / "important" / "figures").glob("lineage_metadata_overlap_*_*.png")))
            self.assertTrue(any((sample_dir / "important" / "figures").glob("lineage_feature_enrichment_*_*.pdf")))
            self.assertTrue((sample_dir / "important" / "tables" / "diversity_feature_richness_by_sample.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "diversity_database_by_sample.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "diversity_database_by_sample_wide.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "diversity_core_common_accessory_rare_features.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "diversity_core_accessory_summary_by_database.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "diversity_pan_feature_accumulation.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "diversity_jaccard_distance_matrix.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "diversity_jaccard_pairs.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "diversity_by_metadata_group.tsv").exists())
            self.assertTrue((sample_dir / "important" / "tables" / "diversity_written_summaries.tsv").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "diversity_analysis.html").exists())
            self.assertTrue((sample_dir / "important" / "diversity_tables.zip").exists())
            self.assertTrue((sample_dir / "important" / "diversity_figures.zip").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "diversity_feature_richness_by_sample.svg").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "diversity_database_by_sample_heatmap.png").exists())
            self.assertTrue((sample_dir / "important" / "figures" / "diversity_pan_feature_accumulation.data.tsv").exists())
            for table_name in [
                "notable_genomes.tsv",
                "notable_genome_score_components.tsv",
                "feature_profile_ordination.tsv",
                "database_concordance_summary.tsv",
                "amr_concordance_feature_level.tsv",
                "amr_concordance_by_sample.tsv",
                "evidence_summary.tsv",
                "finding_confidence_summary.tsv",
                "evidence_by_section.tsv",
                "report_highlights.tsv",
                "report_highlights_by_section.tsv",
                "warning_priority_summary.tsv",
                "report_visual_index.tsv",
                "report_visual_quality.tsv",
                "warnings_and_limitations_summary.tsv",
                "warnings_and_limitations.tsv",
                "warnings_by_section.tsv",
                "module_warning_summary.tsv",
                "report_cap_summary.tsv",
                "important_file_index.tsv",
                "download_manifest.tsv",
            ]:
                self.assertTrue((sample_dir / "important" / "tables" / table_name).exists())
            self.assertTrue((sample_dir / "important" / "downloads" / "important_summary_tables.zip").exists())
            self.assertTrue((sample_dir / "important" / "downloads" / "important_tables.zip").exists())
            self.assertTrue((sample_dir / "important" / "downloads" / "important_figures.zip").exists())
            self.assertTrue((sample_dir / "important" / "downloads" / "publication_candidate_figures.zip").exists())
            self.assertTrue((sample_dir / "important" / "downloads" / "important_report_assets.zip").exists())
            for figure_name in [
                "notable_genomes_ranked.svg",
                "notable_genome_score_heatmap.png",
                "feature_profile_pcoa_by_lineage.svg",
                "amr_concordance_summary.svg",
                "evidence_confidence_summary.svg",
                "warnings_summary.svg",
            ]:
                self.assertTrue((sample_dir / "important" / "figures" / figure_name).exists())
            prevalence_svg = (sample_dir / "important" / "figures" / "prevalence_genomes_positive_by_database.svg").read_text(encoding="utf-8")
            self.assertIn("role='img'", prevalence_svg)
            self.assertIn("<desc id=", prevalence_svg)
            self.assertIn("positive and total genome counts", prevalence_svg)
            self.assertIn("#dc2626", prevalence_svg)
            temporal_slope_svg = (sample_dir / "important" / "figures" / "temporal_slope_top40.svg").read_text(encoding="utf-8")
            self.assertIn("Red indicates increasing prevalence", temporal_slope_svg)
            self.assertIn("#dc2626", temporal_slope_svg)
            self.assertIn("#2563eb", temporal_slope_svg)
            variation_scatter = next((sample_dir / "important" / "figures").glob("variation_identity_coverage_*_top20.svg")).read_text(encoding="utf-8")
            self.assertIn("90% identity", variation_scatter)
            self.assertIn("80% coverage", variation_scatter)
            metadata_volcano = next((sample_dir / "important" / "figures").glob("metadata_volcano_*_*.svg")).read_text(encoding="utf-8").lower()
            self.assertNotIn("upregulated", metadata_volcano)
            self.assertNotIn("downregulated", metadata_volcano)
            temporal_summary = pd.read_csv(sample_dir / "important" / "key_tables" / "temporal_trend_summary.tsv", sep="\t")
            self.assertTrue({"trend_label", "support_label", "temporal_pattern_label", "warning_flags"}.issubset(temporal_summary.columns))
            prevalence = pd.read_csv(sample_dir / "important" / "tables" / "feature_prevalence.tsv", sep="\t")
            self.assertTrue({"positive_genomes", "total_genomes", "prevalence_percent", "feature_rows", "mean_hits_per_positive_genome", "prevalence_label", "warning_flags"}.issubset(prevalence.columns))
            prevalence_db = pd.read_csv(sample_dir / "important" / "tables" / "prevalence_summary_by_database.tsv", sep="\t")
            for _, db_row in prevalence_db.iterrows():
                self.assertGreaterEqual(int(db_row["positive_genomes"]), int(db_row["top_feature_positive_genomes"]))
                self.assertGreaterEqual(float(db_row["genomes_positive_percent"]), float(db_row["top_feature_prevalence_percent"]))
            prevalence_html = (sample_dir / "important" / "figures" / "prevalence_analysis.html").read_text(encoding="utf-8")
            for control in ["Database", "Top 10", "Top 20", "Top 50", "Complete", "Genome prevalence %", "Minimum positive genomes"]:
                self.assertIn(control, prevalence_html)
            geographic_html = (sample_dir / "important" / "figures" / "geographic_distribution.html").read_text(encoding="utf-8")
            for control in ["Database", "Mode", "Feature", "Geographic level", "Metric", "Minimum group size", "Display", "Warning filter"]:
                self.assertIn(control, geographic_html)
            self.assertIn("Search feature / gene", geographic_html)
            self.assertIn("featureSearchInput", geographic_html)
            for guide_text in ["Map reading guide", "Higher selected-gene prevalence", "positive / total genomes", "small group"]:
                self.assertIn(guide_text, geographic_html)
            self.assertIn("plotly.min.js", geographic_html)
            self.assertIn("Plotly.react", geographic_html)
            self.assertIn("type: 'choropleth'", geographic_html)
            self.assertIn("locationmode: 'country names'", geographic_html)
            self.assertIn("natural earth", geographic_html)
            self.assertIn("map-workspace", geographic_html)
            self.assertIn("mapPlot", geographic_html)
            self.assertNotIn("country tile", geographic_html.lower())
            geographic_burden = pd.read_csv(sample_dir / "important" / "tables" / "geographic_database_burden.tsv", sep="\t")
            self.assertTrue({"database", "geo_level", "group_name", "positive_genomes", "prevalence_percent", "warning_flags"}.issubset(geographic_burden.columns))
            country_rows = geographic_burden[(geographic_burden["database"] == "amr") & (geographic_burden["geo_level"] == "country")]
            bd_row = country_rows[country_rows["group_name"] == "Bangladesh"].iloc[0]
            self.assertEqual(int(bd_row["positive_genomes"]), 1)
            self.assertEqual(int(bd_row["total_genomes"]), 1)
            self.assertAlmostEqual(float(bd_row["prevalence_percent"]), 100.0)
            self.assertIn("small_group_warning", str(bd_row["warning_flags"]))
            variation_summary = pd.read_csv(sample_dir / "important" / "key_tables" / "feature_variation_summary.tsv", sep="\t")
            self.assertTrue({"feature_name", "prevalence_percent", "mean_hits_per_positive_genome", "median_alignment_length", "iqr_alignment_length", "variation_score"}.issubset(variation_summary.columns))
            variation_hits = pd.read_csv(sample_dir / "important" / "key_tables" / "feature_variation_hits.tsv", sep="\t")
            self.assertIn("alignment_length", variation_hits.columns)
            variation_html = (sample_dir / "important" / "figures" / "variation_analysis.html").read_text(encoding="utf-8")
            for control in ["Database", "Metric", "View", "Display", "Sort by", "Identity", "Coverage", "Alignment length", "Feature count per genome"]:
                self.assertIn(control, variation_html)
            temporal_html = (sample_dir / "important" / "figures" / "temporal_trends.html").read_text(encoding="utf-8")
            for control in ["Database", "Trend", "Support", "Feature", "Selected Feature Prevalence", "First-to-Last Year Slope"]:
                self.assertIn(control, temporal_html)
            cooccurrence_html = (sample_dir / "important" / "figures" / "cooccurrence_context.html").read_text(encoding="utf-8")
            for control in [
                "Analysis mode",
                "X database",
                "Y database",
                "Feature set",
                "Minimum sample support",
                "Minimum feature prevalence",
                "Significance",
                "Effect size",
                "Evidence level",
                "Download all co-occurrence tables ZIP",
                "Download all co-occurrence figures ZIP",
            ]:
                self.assertIn(control, cooccurrence_html)
            cooccurrence_pairs = pd.read_csv(sample_dir / "important" / "tables" / "cooccurrence_pair_summary.tsv", sep="\t")
            self.assertTrue({"phi_correlation", "q_value", "significance_label", "evidence_level", "warning_flags"}.issubset(cooccurrence_pairs.columns))
            metadata_html = (sample_dir / "important" / "figures" / "metadata_associations.html").read_text(encoding="utf-8")
            for control in ["Database", "Association type", "Metadata variable", "Group", "Minimum group size", "Significance", "Effect size", "Warning filter", "Display"]:
                self.assertIn(control, metadata_html)
            metadata_enrichment = pd.read_csv(sample_dir / "important" / "tables" / "metadata_feature_enrichment.tsv", sep="\t")
            self.assertTrue({"odds_ratio", "p_value", "q_value", "support_label", "interpretation_label", "warning_flags"}.issubset(metadata_enrichment.columns))
            metadata_usability = pd.read_csv(sample_dir / "important" / "tables" / "metadata_usability_summary.tsv", sep="\t")
            self.assertTrue({"metadata_column", "non_missing_count", "missing_fraction", "eligible_for_testing", "recommended_use"}.issubset(metadata_usability.columns))
            metadata_omnibus = pd.read_csv(sample_dir / "important" / "tables" / "metadata_burden_omnibus.tsv", sep="\t")
            self.assertTrue({"test_name", "test_statistic", "p_value", "q_value", "support_label", "interpretation_label"}.issubset(metadata_omnibus.columns))
            lineage_html = (sample_dir / "important" / "figures" / "lineage_clonal_structure.html").read_text(encoding="utf-8")
            for control in ["Lineage type", "Database", "Feature mode", "Search feature", "Feature", "Metadata overlay", "Minimum lineage size", "custom", "Display", "Written Summaries"]:
                self.assertIn(control, lineage_html)
            lineage_distribution = pd.read_csv(sample_dir / "important" / "tables" / "lineage_distribution.tsv", sep="\t")
            self.assertTrue({"lineage_type", "lineage_id", "total_genomes", "fraction_of_dataset", "warning_flags"}.issubset(lineage_distribution.columns))
            lineage_enrichment = pd.read_csv(sample_dir / "important" / "tables" / "lineage_feature_enrichment.tsv", sep="\t")
            self.assertTrue({"odds_ratio", "p_value", "q_value", "support_label", "interpretation_label", "warning_flags"}.issubset(lineage_enrichment.columns))
            lineage_written = pd.read_csv(sample_dir / "important" / "tables" / "lineage_written_summaries.tsv", sep="\t")
            self.assertIn("lineage_adjusted_top_findings_summary", set(lineage_written["section"]))
            lineage_summary_text = " ".join(lineage_written["summary"].fillna("").astype(str))
            self.assertNotIn("not available (not available)", lineage_summary_text)
            self.assertNotIn("among 0 genome(s)", lineage_summary_text)
            self.assertNotIn("report-facing burden row", lineage_summary_text)
            diversity_html = (sample_dir / "important" / "figures" / "diversity_analysis.html").read_text(encoding="utf-8")
            for control in ["Diversity scope", "Diversity view", "Metadata color/group", "Display", "Sort", "Core/common/accessory/rare", "Pan-feature accumulation", "Jaccard similarity/distance"]:
                self.assertIn(control, diversity_html)
            diversity_richness = pd.read_csv(sample_dir / "important" / "tables" / "diversity_feature_richness_by_sample.tsv", sep="\t")
            self.assertTrue({"assembly_accession", "total_unique_features", "total_feature_rows", "richness_label", "warning_flags"}.issubset(diversity_richness.columns))
            diversity_classes = pd.read_csv(sample_dir / "important" / "tables" / "diversity_core_common_accessory_rare_features.tsv", sep="\t")
            self.assertTrue({"database", "feature_id", "positive_genomes", "prevalence_percent", "feature_class", "feature_rows"}.issubset(diversity_classes.columns))
            diversity_jaccard = pd.read_csv(sample_dir / "important" / "tables" / "diversity_jaccard_distance_matrix.tsv", sep="\t")
            self.assertTrue({"sample_a", "sample_b", "jaccard_distance", "jaccard_similarity"}.issubset(diversity_jaccard.columns))
            diversity_summary = pd.read_csv(sample_dir / "important" / "tables" / "diversity_report_summary.tsv", sep="\t")
            self.assertTrue({"total_feature_rows", "max_features_in_one_genome", "databases_represented", "jaccard_matrix_available", "pan_feature_curve_available"}.issubset(set(diversity_summary["metric"])))
            notable = pd.read_csv(sample_dir / "important" / "tables" / "notable_genomes.tsv", sep="\t")
            self.assertTrue({"notable_genome_score", "notable_label", "score_explanation", "warning_flags"}.issubset(notable.columns))
            self.assertNotIn("clinical risk", " ".join(notable["score_explanation"].fillna("").astype(str)).lower())
            components = pd.read_csv(sample_dir / "important" / "tables" / "notable_genome_score_components.tsv", sep="\t")
            top_sample = notable.sort_values("rank").iloc[0]["assembly_accession"]
            component_sum = components[components["assembly_accession"] == top_sample]["component_score"].astype(float).sum()
            top_score = float(notable[notable["assembly_accession"] == top_sample]["notable_genome_score"].iloc[0])
            self.assertAlmostEqual(component_sum, top_score, places=2)
            temporal_components = components[components["component"] == "temporal_increasing_feature_score"]
            self.assertTrue((temporal_components["component_score"].astype(float) >= 0).all())
            ordination = pd.read_csv(sample_dir / "important" / "tables" / "feature_profile_ordination.tsv", sep="\t")
            self.assertTrue({"PCoA1", "PCoA2", "explained_variance_PCoA1", "explained_variance_PCoA2"}.issubset(ordination.columns))
            concordance = pd.read_csv(sample_dir / "important" / "tables" / "amr_concordance_feature_level.tsv", sep="\t")
            self.assertTrue({"called_by_both", "abricate_only", "amrfinderplus_only", "concordance_label"}.issubset(concordance.columns))
            confidence = pd.read_csv(sample_dir / "important" / "tables" / "finding_confidence_summary.tsv", sep="\t")
            self.assertTrue({"confidence_label", "recommended_interpretation"}.issubset(confidence.columns))
            highlights = pd.read_csv(sample_dir / "important" / "tables" / "report_highlights.tsv", sep="\t")
            self.assertTrue({"rank", "section", "highlight_type", "triage_score", "confidence_label", "recommended_action"}.issubset(highlights.columns))
            if not highlights.empty:
                self.assertTrue(highlights["confidence_label"].fillna("").astype(str).str.len().gt(0).all())
            self.assertNotIn("database_burden", set(highlights.head(20)["primary_feature"].fillna("").astype(str)))
            informative_pairs = highlights[highlights["highlight_type"] == "informative_cooccurrence"]
            for _, row in informative_pairs.iterrows():
                self.assertFalse(_is_same_feature_pair(str(row.get("primary_feature", "")), str(row.get("secondary_feature", ""))))
            highlight_sections = set(highlights.head(20)["section"].fillna("").astype(str))
            self.assertGreaterEqual(len(highlight_sections), 2)
            all_highlight_sections = set(highlights["section"].fillna("").astype(str))
            self.assertIn("Diversity / Pan-feature Summary", all_highlight_sections)
            self.assertIn("Warnings & Limitations", all_highlight_sections)
            self.assertIn("Concordance / Database Agreement", all_highlight_sections)
            if confidence["confidence_label"].isin(["high_confidence", "moderate_confidence", "warning_heavy"]).any():
                self.assertIn("Evidence & Confidence", all_highlight_sections)
            by_section = pd.read_csv(sample_dir / "important" / "tables" / "report_highlights_by_section.tsv", sep="\t")
            self.assertTrue({"section_rank", "section", "highlight_type", "triage_score", "confidence_label"}.issubset(by_section.columns))
            if not by_section.empty:
                self.assertTrue(by_section["confidence_label"].fillna("").astype(str).str.len().gt(0).all())
            warning_priorities = pd.read_csv(sample_dir / "important" / "tables" / "warning_priority_summary.tsv", sep="\t")
            self.assertTrue({"rank", "section", "severity", "warning_type", "priority_score", "why_it_matters"}.issubset(warning_priorities.columns))
            visual_index = pd.read_csv(sample_dir / "important" / "tables" / "report_visual_index.tsv", sep="\t")
            self.assertTrue({
                "figure_stem",
                "section",
                "interpretation_type",
                "svg_path",
                "recommended_use",
                "title",
                "default_visibility",
                "display_reason",
                "recommended_audience",
                "main_report_priority",
                "publication_candidate",
            }.issubset(visual_index.columns))
            self.assertTrue(set(visual_index["default_visibility"].dropna().astype(str)).issubset({"featured", "standard", "supporting", "technical"}))
            self.assertIn("featured", set(visual_index["default_visibility"].fillna("").astype(str)))
            self.assertIn("true", {value.lower() for value in visual_index["publication_candidate"].fillna("").astype(str)})
            network_rows = visual_index[visual_index["figure_stem"].fillna("").astype(str).str.startswith("cooccurrence_network_")]
            if not network_rows.empty:
                self.assertTrue(set(network_rows["default_visibility"].fillna("").astype(str)).issubset({"technical"}))
            pcoa_rows = visual_index[visual_index["figure_stem"] == "feature_profile_pcoa_by_bioproject"]
            if not pcoa_rows.empty:
                self.assertEqual(pcoa_rows.iloc[0]["title"], "Feature-profile PCoA by BioProject")
            concordance_rows = visual_index[visual_index["figure_stem"] == "amr_concordance_summary"]
            if not concordance_rows.empty:
                self.assertEqual(concordance_rows.iloc[0]["title"], "AMR tool concordance summary")
            geographic_rows = visual_index[visual_index["figure_stem"] == "geographic_distribution_map"]
            if not geographic_rows.empty:
                self.assertEqual(geographic_rows.iloc[0]["title"], "Geographic gene map")
            visual_quality = pd.read_csv(sample_dir / "important" / "tables" / "report_visual_quality.tsv", sep="\t")
            self.assertTrue({
                "figure_stem",
                "quality_label",
                "svg_available",
                "data_tsv_available",
                "asset_quality_label",
                "interpretation_quality_label",
                "title_quality_label",
                "caption_quality_label",
                "final_publication_label",
                "default_visibility",
                "display_reason",
                "recommended_audience",
                "main_report_priority",
                "publication_candidate",
            }.issubset(visual_quality.columns))
            exploratory_quality = visual_quality[visual_quality["section"].isin(["Geographic Distribution", "Co-occurrence / Genomic Context", "Lineage / Clonal Structure"])]
            if not exploratory_quality.empty:
                self.assertNotIn("publication_ready", set(exploratory_quality["final_publication_label"].fillna("").astype(str)))
            variation_quality = visual_quality[visual_quality["figure_stem"].fillna("").astype(str).str.startswith("variation_")]
            if not variation_quality.empty:
                self.assertNotIn("generic_caption", set(variation_quality["caption_quality_label"].fillna("").astype(str)))
            qc_quality = visual_quality[visual_quality["figure_stem"].fillna("").astype(str).str.startswith("qc_")]
            if not qc_quality.empty:
                self.assertNotIn("generic_caption", set(qc_quality["caption_quality_label"].fillna("").astype(str)))
                self.assertNotIn("asset_incomplete", set(qc_quality["asset_quality_label"].fillna("").astype(str)))
            temporal_quality = visual_quality[visual_quality["figure_stem"].fillna("").astype(str).str.startswith("temporal_")]
            if not temporal_quality.empty:
                self.assertNotIn("missing_axis_labels", set(temporal_quality["axis_label_status"].fillna("").astype(str)))
                self.assertNotIn("asset_incomplete", set(temporal_quality["asset_quality_label"].fillna("").astype(str)))
            geographic_map_quality = visual_quality[visual_quality["figure_stem"].fillna("").astype(str) == "geographic_distribution_map"]
            if not geographic_map_quality.empty:
                self.assertEqual(geographic_map_quality.iloc[0]["asset_quality_label"], "asset_ready")
            warnings = pd.read_csv(sample_dir / "important" / "tables" / "warnings_and_limitations.tsv", sep="\t")
            self.assertTrue({"warning_id", "section", "severity", "warning_type", "recommended_action"}.issubset(warnings.columns))
            self.assertNotIn("Warning flags reported in", "\n".join(warnings["description"].fillna("").astype(str)))
            self.assertTrue(warnings["recommended_action"].fillna("").astype(str).str.len().gt(20).all())
            warning_summary = pd.read_csv(sample_dir / "important" / "tables" / "warnings_and_limitations_summary.tsv", sep="\t")
            self.assertTrue({"section", "severity", "warning_type", "warning_count", "recommended_action"}.issubset(warning_summary.columns))
            warning_priorities = pd.read_csv(sample_dir / "important" / "tables" / "warning_priority_summary.tsv", sep="\t")
            if not warning_priorities.empty:
                self.assertNotEqual(
                    set(warning_priorities["why_it_matters"].fillna("").astype(str)),
                    {"This warning can change whether a finding is broad, lineage-driven, BioProject-driven, or under-supported."},
                )
            download_manifest = pd.read_csv(sample_dir / "important" / "tables" / "download_manifest.tsv", sep="\t")
            self.assertIn("basic/enriched_genome_dataset.csv", set(download_manifest["file_path"]))
            all_index_rows = download_manifest[download_manifest["file_path"] == "all/file_index.tsv"]
            if not all_index_rows.empty:
                self.assertEqual(all_index_rows.iloc[0]["complete_or_capped"], "not_generated_in_important_mode")
            report_controls = pd.read_csv(outputs["report_controls"], sep="\t")
            self.assertIn("important_lineage_feature_cap_per_database", set(report_controls["setting"]))
            self.assertIn("important_diversity_jaccard_heatmap_cap", set(report_controls["setting"]))
            report_html = (sample_dir / "important" / "results.html").read_text(encoding="utf-8")
            self.assertNotIn("Report-facing figure with PNG", report_html)
            self.assertNotIn("Warning rows represent", report_html)
            self.assertNotIn("warning rows:", report_html)
            self.assertIn("Warning flags average", report_html)
            self.assertIn("one comparison may carry several warnings", report_html)
            self.assertIn("Look for:", report_html)
            self.assertIn("Caution:", report_html)
            self.assertIn("sidebar-links", report_html)
            self.assertIn("analysis-card", report_html)
            self.assertIn("grid-template-columns: repeat(auto-fit, minmax(150px, 1fr))", report_html)
            self.assertIn("Diagnostics and unavailable plots", report_html)
            self.assertIn("Report storyline", report_html)
            self.assertIn("Best figure to start with", report_html)
            self.assertIn("interactive-explorer-card", report_html)
            self.assertIn("explorer-frame", report_html)
            self.assertIn("Load embedded explorer", report_html)
            self.assertIn("Open filled-country gene map", report_html)
            self.assertGreaterEqual(report_html.count("loading='lazy'"), 8)
            if by_section["section"].nunique() >= 5:
                balanced_preview = report_html.split("<h3>Balanced highlights by section</h3>", 1)[1].split('<div class="downloads"', 1)[0]
                preview_sections = [
                    section
                    for section in sorted(set(by_section["section"].fillna("").astype(str)))
                    if section and section in balanced_preview
                ]
                self.assertGreaterEqual(len(preview_sections), 5)
            for section in [
                "Featured Results",
                "Run Overview",
                "QC Summary",
                "Enriched Dataset",
                "Prevalence",
                "Geographic Distribution",
                "Variations",
                "Temporal Trends",
                "Co-occurrence / Genomic Context",
                "Metadata Associations",
                "Lineage / Clonal Structure",
                "Diversity / Pan-feature Summary",
                "Notable Genomes / Genome Prioritization",
                "Feature-profile Ordination",
                "Concordance / Database Agreement",
                "Evidence & Confidence",
                "Warnings & Limitations",
                "Downloads / Important Files",
            ]:
                self.assertIn(section, report_html)
            for anchor in [
                'id="featured"',
                'id="overview"',
                'id="qc"',
                'id="enriched-dataset"',
                'id="prevalence"',
                'id="geography"',
                'id="variations"',
                'id="temporal"',
                'id="cooccurrence"',
                'id="metadata-associations"',
                'id="lineage"',
                'id="diversity"',
                'id="notable-genomes"',
                'id="ordination"',
                'id="concordance"',
                'id="evidence"',
                'id="warnings"',
                'id="downloads"',
            ]:
                self.assertIn(anchor, report_html)
            for ui_class in [
                "report-header",
                "header-meta",
                "sidebar",
                "section-header",
                "summary-card",
                "figure-card",
                "figure-guidance",
                "section-focus",
                "best-figure-note",
                "interactive-explorer-card",
                "explorer-frame",
                "report-storyline",
                "table-card",
                "table-actions",
                "table-download-link",
                "table-search",
                "warning-box",
                "download-card",
                "back-to-top",
                "details-block",
            ]:
                self.assertIn(ui_class, report_html)
            for css_token in [
                "--primary: #0f766e",
                "--red: #dc2626",
                "--green: #16a34a",
                "overflow-x: hidden",
                "@media (max-width: 920px)",
                "@media (max-width: 520px)",
                "overflow-wrap: anywhere",
                "main > .report-header, main > .section",
                "width: auto !important",
                ".figure-guidance",
                ".section-focus",
                ".best-figure-note",
                ".interactive-explorer-card",
                ".explorer-frame",
                ".report-storyline",
            ]:
                self.assertIn(css_token, report_html)
            self.assertIn("Featured figure gallery", report_html)
            self.assertIn("How to read this figure", report_html)
            self.assertIn("Next table:", report_html)
            self.assertIn("Section interpretation guide", report_html)
            self.assertIn("What this section answers", report_html)
            self.assertIn("Open next", report_html)
            self.assertIn("Featured results are triage shortcuts", report_html)
            self.assertIn(
                "Lineage summaries are exploratory and do not replace phylogenetic analysis",
                report_html,
            )
            self.assertIn("download_manifest.tsv and important_file_index.tsv", report_html)
            self.assertIn("metadata_feature_enrichment.tsv", report_html)
            self.assertIn("lineage_adjusted_top_findings.tsv", report_html)
            self.assertIn("genomic_context_evidence.tsv", report_html)
            self.assertIn("feature_profile_ordination.tsv", report_html)
            self.assertIn("Featured Results reading guide", report_html)
            self.assertIn("Executive finding cards", report_html)
            self.assertIn("Trust, caution, and review queues", report_html)
            self.assertIn("Balanced highlight table by section", report_html)
            self.assertIn("Download full TSV", report_html)
            self.assertIn("tables/feature_prevalence.tsv", report_html)
            self.assertIn("tables/report_visual_index.tsv", report_html)
            self.assertIn("What to review first", report_html)
            self.assertNotIn(">1-01-01<", report_html)
            self.assertNotIn(">0001-01-01<", report_html)
            self.assertIn("What to trust first", report_html)
            self.assertIn("What needs caution", report_html)
            self.assertIn("Balanced highlights by section", report_html)
            for report_phrase in [
                "Warning triage cards",
                "Prevalence reading guide",
                "Detailed prevalence tables",
                "Geographic distribution reading guide",
                "Detailed geographic tables",
                "Variation reading guide",
                "Detailed variation tables",
                "Temporal trends reading guide",
                "Detailed temporal tables",
                "Co-occurrence and context reading guide",
                "Detailed co-occurrence and context tables",
                "Metadata association reading guide",
                "Detailed metadata association tables",
                "Lineage reading guide",
                "Detailed lineage tables",
                "Detailed warning tables",
                "Main files",
                "Review and interpretation tables",
                "Reproducibility and audit files",
                "Visual index and quality previews",
                "Important file index preview",
            ]:
                self.assertIn(report_phrase, report_html)
            geography_section = report_html.split('id="geography"', 1)[1].split('id="variations"', 1)[0]
            self.assertIn("Filled-country geographic gene map", geography_section)
            self.assertIn("figures/geographic_distribution.html", geography_section)
            self.assertIn("Embedded explorer", geography_section)
            self.assertIn("geographic_country_bar", geography_section)
            self.assertLess(
                geography_section.index("figures/geographic_distribution.html"),
                geography_section.index("Geographic denominator, missingness, and reading guide"),
            )
            self.assertLess(
                geography_section.index("figures/geographic_distribution.html"),
                geography_section.index("geographic_country_bar"),
            )
            metadata_section = report_html.split('id="metadata-associations"', 1)[1].split('id="lineage"', 1)[0]
            self.assertIn("metadata_volcano_", metadata_section)
            self.assertIn("metadata_enrichment_heatmap_", metadata_section)
            self.assertIn("metadata_burden_boxplot_", metadata_section)
            self.assertLess(
                metadata_section.index("metadata_volcano_"),
                metadata_section.index("metadata_burden_boxplot_"),
            )
            lineage_section = report_html.split('id="lineage"', 1)[1].split('id="diversity"', 1)[0]
            self.assertIn("lineage_distribution_", lineage_section)
            self.assertIn("lineage_metadata_overlap_", lineage_section)
            self.assertIn("lineage_feature_heatmap_", lineage_section)
            self.assertIn("lineage_confounding_top_findings", lineage_section)
            self.assertLess(
                lineage_section.index("lineage_feature_heatmap_"),
                lineage_section.index("lineage_database_burden_"),
            )
            self.assertIn("Download report highlights", report_html)
            self.assertIn("Visual index", report_html)
            self.assertIn("Visual quality", report_html)
            self.assertIn("Download enriched dataset", report_html)
            self.assertIn("Download summary tables ZIP", report_html)
            self.assertIn("Download complete tables ZIP", report_html)
            self.assertIn("Download important figures ZIP", report_html)
            self.assertIn("Download publication candidates ZIP", report_html)
            self.assertIn("Download report assets ZIP", report_html)
            self.assertIn("More supporting", report_html)
            self.assertIn("Interpret these patterns as dataset-specific sampling summaries", report_html)
            self.assertIn("Warning flags average", report_html)
            self.assertIn("PNG</a>", report_html)
            self.assertIn("SVG</a>", report_html)
            self.assertIn("Data TSV</a>", report_html)
            self.assertIn("alt='", report_html)
            self.assertIn("Showing ", report_html)
            self.assertIn("download the full TSV", report_html)
            for zip_name in ["important_summary_tables.zip", "important_tables.zip", "important_figures.zip", "publication_candidate_figures.zip", "important_report_assets.zip"]:
                self.assertGreater((sample_dir / "important" / "downloads" / zip_name).stat().st_size, 0)
            qa = subprocess.run(
                [sys.executable, str(REPO_ROOT / "scripts" / "check_important_report_outputs.py"), str(sample_dir)],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(qa.returncode, 0, qa.stderr + qa.stdout)
            visual_qa = subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "check_important_report_visual_layout.py"),
                    str(sample_dir),
                    "--browser",
                    "skip",
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(visual_qa.returncode, 0, visual_qa.stderr + visual_qa.stdout)

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
            important_context = pd.read_csv(sample_dir / "important" / "tables" / "genomic_context_evidence.tsv", sep="\t")
            self.assertFalse(important_context.empty)
            self.assertIn("within_10kb", set(important_context["evidence_level"]))
            neighborhoods = pd.read_csv(sample_dir / "important" / "tables" / "contig_neighborhoods.tsv", sep="\t")
            self.assertFalse(neighborhoods.empty)
            self.assertTrue((sample_dir / "important" / "figures" / "cooccurrence_context.html").exists())
            self.assertTrue(any((sample_dir / "important" / "figures").glob("genomic_context_evidence_ladder_*.svg")))

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
            important_lineage = pd.read_csv(sample_dir / "important" / "tables" / "lineage_distribution.tsv", sep="\t")
            self.assertIn("ST_11", set(important_lineage["lineage_id"]))
            important_presence = pd.read_csv(sample_dir / "important" / "tables" / "lineage_feature_presence.tsv", sep="\t")
            lineage_gene = important_presence[(important_presence["feature_id"] == "lineageGene") & (important_presence["lineage_id"] == "ST_11")].iloc[0]
            self.assertEqual(int(lineage_gene["positive_genomes"]), 3)
            self.assertIn("Lineage type", (sample_dir / "important" / "figures" / "lineage_clonal_structure.html").read_text(encoding="utf-8"))

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

    def test_setup_abricate_databases_uses_external_datadir_without_force_update(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            fake_bin = root / "bin"
            fake_bin.mkdir()
            db_dir = root / "shared_abricate_db"
            db_dir.mkdir()
            log = root / "commands.log"
            (fake_bin / "panr").write_text(
                "#!/bin/sh\n"
                f"printf 'panr %s\\n' \"$*\" >> {log}\n"
                "printf 'unexpected panr call\\n'\n",
                encoding="utf-8",
            )
            (fake_bin / "abricate-get_db").write_text(
                "#!/bin/sh\n"
                f"printf 'abricate-get_db %s\\n' \"$*\" >> {log}\n"
                "printf 'unexpected update call\\n'\n",
                encoding="utf-8",
            )
            (fake_bin / "abricate").write_text(
                "#!/bin/sh\n"
                f"printf 'abricate %s\\n' \"$*\" >> {log}\n"
                "if [ \"$1\" = \"--datadir\" ]; then shift 2; fi\n"
                "if [ \"$1\" = \"--list\" ]; then\n"
                "  printf 'DATABASE\\tSEQUENCES\\n'\n"
                "  printf 'ncbi\\t10\\n'\n"
                "  printf 'vfdb\\t10\\n'\n"
                "elif [ \"$1\" = \"--setupdb\" ]; then\n"
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
                    "--datadir",
                    str(db_dir),
                    "--update",
                ],
                check=True,
                env=env,
            )

            status = pd.read_csv(out, sep="\t")
            self.assertEqual(set(status["status"]), {"PASS"})
            self.assertTrue(status["database_dir"].astype(str).str.contains(str(db_dir)).all())
            self.assertEqual(set(status["update_status"]), {"SKIPPED"})
            commands = log.read_text(encoding="utf-8")
            self.assertIn(f"abricate --datadir {db_dir} --list", commands)
            self.assertIn(f"abricate --datadir {db_dir} --setupdb", commands)
            self.assertNotIn("panr setup-db", commands)
            self.assertNotIn("abricate-get_db", commands)

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
