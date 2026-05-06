# PanResistome: Scalable Pipeline for Global Antimicrobial Resistance Analysis

Current pipeline version: `0.2.1`

## Overview

**PanResistome** is a scalable, modular, and reproducible bioinformatics pipeline built using [Nextflow](https://www.nextflow.io/). It automates the end-to-end analysis of global antimicrobial resistance (AMR) patterns in bacterial populations using genome assemblies. The pipeline is designed for researchers working in microbial genomics, resistome surveillance, and public health, enabling large-scale comparative analysis of resistance gene profiles across time and geography.

PanResistome integrates several state-of-the-art tools including:

* [**FetchM2**](https://github.com/Tasnimul-Arabi-Anik/FetchM2): for standardized genome metadata, host/source/environment normalization, metadata audit reports, and download-ready assembly tables
* [**CheckM2**](https://github.com/chklovski/CheckM2): for genome completeness and contamination quality assessment
* [**GTDB-Tk**](https://github.com/Ecogenomics/GTDBTk): for taxonomy classification and genus/species consistency checks
* [**QUAST**](https://github.com/ablab/quast): for assembly structure metrics such as N50, contig count, GC, and total length
* [**FastANI/skani**](https://github.com/ParBLiSS/FastANI): for pairwise ANI, species-consistency screening, and near-duplicate clustering
* [**Mash**](https://github.com/marbl/Mash): for optional fast sketch-based distance pre-screening
* [**ABRicate**](https://github.com/tseemann/abricate): for AMR, virulence, and plasmid database annotation through the PanR2 comprehensive mode
* [**MobileElementFinder**](https://bitbucket.org/genomicepidemiology/mobileelementfinder): for mobile genetic element annotation in comprehensive mode
* [**IntegronFinder**](https://github.com/gem-pasteur/Integron_Finder): for integron annotation in comprehensive mode
* [**MLST**](https://github.com/tseemann/mlst): for sequence-type context in comprehensive mode
* [**MOB-suite**](https://github.com/phac-nml/mob-suite): optional plasmid reconstruction/typing before PanR2 analysis
* [**geNomad**](https://github.com/apcamargo/genomad): optional prophage/viral-region annotation when a geNomad database is provided
* **Organism-specific typing modules**: optional Kleborate/Kaptive/ECTyper table generation or ingestion for organism-focused comparative genomics
* [**PanR2**](https://github.com/Tasnimul-Arabi-Anik/PanR2): for downstream statistical analysis, cross-database associations, reports, and interactive visualization

## Key Features

* 🔄 **Fully automated** end-to-end pipeline from genome download to visualization
* 🧬 **Panresistome analysis** using resistance gene profiling from ABRicate
* 🧫 **Comprehensive feature mode** for NCBI AMR, VFDB, PlasmidFinder, MobileElementFinder, IntegronFinder, MLST, optional MOB-suite, prophage/geNomad, and organism-specific typing tables through PanR2
* 📊 **Visualization-ready outputs** including heatmaps, barplots, boxplots, and interactive HTML figures
* 📈 **Statistical summaries** and correlation-based insights on resistance gene distribution
* 🌍 **Geospatial & temporal comparison** of AMR gene prevalence
* 💡 **Epidemic signal detection** by comparing ARG prevalence across time and location
* ⚙️ **Nextflow-based** for reproducibility, scalability, and cloud/HPC compatibility
* 🧾 **PanR2 handoff exports** so heavy tool execution stays in PanResistome while PanR2 stays lightweight

---

## Workflow Overview

```
+-------------+   +-------------+   +----------------------+   +----------------------+   +-------+
|   FetchM2   |-->| Sequence QC |-->| CheckM2/QUAST/ANI/QC |-->| PanR2 comprehensive  |-->| Report|
+-------------+   +-------------+   +----------------------+   +----------------------+   +-------+
| Assemblies  |   | Assembly    |   | Heavy comparative    |   | ABRicate databases, |   | plots |
| & metadata  |   | stats       |   | genomics + filtering |   | MGE, integron, MLST |   | HTML  |
+-------------+   +-------------+   +----------------------+   +----------------------+   +-------+

Optional: add `--run_gtdbtk true`, `--run_quast true`, `--run_ani true`, or `--run_mash true` to insert heavier comparative-genomics checks before ABRicate and PanR2.
```

Each run also writes Conda environment version reports under `pipeline_versions/` so analyses can be traced back to the exact tool versions used.

Sequence QC filtering is optional. By default, the pipeline reports QC metrics but keeps all assemblies for downstream analysis. Add `--qc_filter true` with one or more thresholds to exclude failed assemblies from ABRicate, PanR2, and later tools.

CheckM2 is enabled by default. If `--checkm2_db` is not provided, PanResistome now attempts to download the CheckM2 database automatically under `<outdir>/databases/checkm2` unless `--checkm2_auto_download_db false` is set. This improves the one-command user path, but the database is large, so users on restricted networks can still pre-download it and pass `--checkm2_db /path/to/checkm2_database.dmnd`.

For a lighter validation run on modest hardware, use `--stop_after_qc true`. This runs FetchM2 metadata processing, sequence QC, and CheckM2, then collects QC outputs without running GTDB-Tk, ABRicate, or PanR2.

GTDB-Tk is disabled by default because it is resource-intensive. If enabled, it requires its reference data to be available in the run environment. If `GTDBTK_DATA_PATH` is not configured globally, pass the location with `--gtdbtk_data_path /path/to/gtdbtk_data`. Taxonomy matching compares GTDB-Tk classification against the organism/species metadata at genus rank by default; use `--taxonomy_match_rank species` for stricter matching.

QUAST, FastANI/skani, and Mash are also optional. They are part of PanResistome because they require external tools and can be expensive on large genome sets. Their outputs are summarized into standardized tables and exported under `panr2_inputs/` for PanR2 to analyze without inheriting the heavy dependencies.

For the broader database/tool workflow, add `--run_panr2_comprehensive true`. This makes PanResistome call PanR2's integrated runners from the pinned comprehensive Conda environment. The tested default comprehensive mode runs ABRicate `ncbi`, `vfdb`, and `plasmidfinder`, plus MobileElementFinder, IntegronFinder, and MLST. It also writes PanR2 database-specific folders, cross-database associations, temporal summaries, a top-level dashboard at `report/index.html`, citations, software versions, and a standardized `panr2_inputs/` handoff bundle.

DefenseFinder remains available as `--panr2_run_defensefinder true`, but it is not part of the default comprehensive mode until its Conda dependency stack is stable across fresh installs. GTDB-Tk remains off by default because it requires a large external reference database.

Optional heavy modules can be inserted before PanR2:

```bash
--run_mobsuite true
--run_genomad true --genomad_db /path/to/genomad_db
--run_kleborate true
--run_ectyper true
--run_kaptive true --kaptive_db /path/to/kaptive_db
```

The same feature families can also be supplied as precomputed tables without running the external tools inside PanResistome:

```bash
--mobsuite_dir /path/to/mobsuite_tables
--defensefinder_dir /path/to/defensefinder_tables
--prophage_dir /path/to/prophage_tables
--kleborate_dir /path/to/kleborate_tables
--kaptive_dir /path/to/kaptive_tables
--ectyper_dir /path/to/ectyper_tables
--serotypefinder_dir /path/to/serotypefinder_tables
--sccmecfinder_dir /path/to/sccmecfinder_tables
```

This table-input path is the preferred stable route for organism-specific typing and CGE outputs until their database setup is reproducible across fresh machines. SerotypeFinder and SCCmecFinder are currently supported as PanR2-compatible table inputs; their CGE database-driven runners are intentionally kept outside the default environment.

---

## 🚀 Getting Started

### ✅ Prerequisites

* [Nextflow](https://www.nextflow.io/)
* [Conda](https://docs.conda.io/en/latest/)
* [Mamba](https://mamba.readthedocs.io/) is strongly recommended for optional heavy-tool environments
* Git

### 🔧 Installation

```bash
git clone https://github.com/Tasnimul-Arabi-Anik/PanResistome.git
cd PanResistome
```

### ✅ Preflight Setup

Run the bootstrap script before launching a long analysis. It checks Java, Nextflow, Conda, environment files, ABRicate database files, CheckM2 database path, and Nextflow syntax.

If you already have the CheckM2 database:

```bash
scripts/bootstrap.sh \
  --checkm2-db /path/to/CheckM2_database/uniref100.KO.1.dmnd \
  --abricate-db ./db
```

To download the CheckM2 database:

```bash
scripts/bootstrap.sh \
  --download-checkm2-db "$HOME/databases" \
  --abricate-db ./db
```

The script prints ready-to-run validation and full pipeline commands using the checked paths.

For the most practical fresh-install experience, use `-profile conda,mamba` when mamba is available. Plain `-profile conda` remains supported, but optional heavy environments such as MOB-suite, geNomad, Kleborate, Kaptive, and ECTyper can solve slowly with classic conda on some systems.

### 🧪 Offline CI Test

The repository includes tiny local fixtures for syntax-safe CI and development checks. This profile does not download from NCBI and does not require CheckM2, ABRicate, PanR2, or GTDB-Tk databases.

```bash
nextflow run main.nf -profile test
```

It validates sequence QC, metadata enrichment, pass-only metadata generation, and result collection using `tests/fixtures/local_samples/`.

## Input
Download ncbi_dataset.tsv of your target organism(s) from the [NCBI genome database](https://www.ncbi.nlm.nih.gov/datasets/genome/).
-**ncbi_dataset.tsv**

## 🧪 Running the Pipeline

First validate the download and QC path:

```bash
nextflow run main.nf \
  --input test_small.tsv \
  --outdir results_small \
  -profile conda \
  --stop_after_qc true \
  --run_gtdbtk false \
  --threads 8 \
  --checkm2_db /path/to/CheckM2_database/uniref100.KO.1.dmnd
```

Then run the downstream analysis with QC filtering:

```bash
nextflow run main.nf \
  --input test_small.tsv \
  --outdir results_small \
  -profile conda \
  --run_gtdbtk false \
  --qc_filter true \
  --threads 8 \
  --checkm2_db /path/to/CheckM2_database/uniref100.KO.1.dmnd

```

## 📦 Repository Structure

```
PanResistome/
├── main.nf
├── nextflow.config
├── envs/
│   ├── fetchm.yaml
│   └── abricate.yml
├── results/
├── figures/ 
├── LICENSE
├── test.tsv
└── README.md
```

---

## 🧾 Command-line Options

### ✅ Required Arguments

| Argument   | Description                              |
| ---------- | ---------------------------------------- |
| `--input`  | Input TSV file listing genome accessions |
| `--outdir` | Output directory for results             |

### ⚙️ Optional Arguments for FetchM2

| Argument    | Type   | Default | Description                                           |
| ----------- | ------ | ------- | ----------------------------------------------------- |
| `--checkm`  | float  | -       | Minimum CheckM completeness threshold (e.g., 90)      |
| `--ani`     | str    | all     | ANI filter status: OK, Inconclusive, Failed, or all   |
| `--sleep`   | float  | 0.5     | Time to wait between fetch requests (in seconds)      |
| `--metadata_engine` | str | fetchm2 | Metadata engine: `fetchm2` or reversible legacy `legacy_fetchm` |
| `--fetchm2_offline` | bool | false | Use FetchM2 offline metadata mode |
| `--fetchm2_no_analysis` | bool | false | Skip FetchM2 metadata analysis figures/tables |
| `--fetchm2_download` | bool | true | Download selected assemblies from FetchM2 metadata |
| `--fetchm2_download_engine` | str | native | Sequence downloader: `native` uses `fetchm2 seq`; `panresistome` uses the fallback downloader |
| `--fetchm2_workers` | int | 3 | FetchM2 BioSample metadata workers |
| `--fetchm2_download_workers` | int | 1 | Sequence-download workers; increase cautiously for larger runs |
| `--fetchm2_max_genomes` | int | - | Maximum genomes selected for sequence download after metadata filtering |
| `--fetchm2_keep_assembly_duplicates` | bool | false | Keep paired GCA/GCF assembly rows in `fetchm2_clean.csv`; by default FetchM2 keeps one representative row per Assembly Name |
| `--host`    | str\[] | -       | Host species (e.g., "Homo sapiens", "Bos taurus")     |
| `--year`    | str\[] | -       | Filter by year or range (e.g., "2015" or "2015-2023") |
| `--country` | str\[] | -       | Filter by country (e.g., "Bangladesh", "USA")         |
| `--cont`    | str\[] | -       | Filter by continent (e.g., "Asia", "Africa")          |
| `--subcont` | str\[] | -       | Filter by subcontinent (e.g., "Southern Asia")        |
| `--sample_type` | str\[] | - | FetchM2 `Sample_Type_SD` filter |
| `--isolation_source` | str\[] | - | FetchM2 `Isolation_Source_SD` filter |
| `--environment_medium` | str\[] | - | FetchM2 `Environment_Medium_SD` filter |
| `--year_from` | int | - | Minimum FetchM2 `Collection_Year` |
| `--year_to` | int | - | Maximum FetchM2 `Collection_Year` |

### 🧬 Optional Arguments for PanR2

| Argument   | Type  | Default | Description                                                |
| ---------- | ----- | ------- | ---------------------------------------------------------- |
| `--genep`  | float | -       | Minimum % gene presence to include in heatmaps             |
| `--nseq`   | int   | -       | Minimum number of sequences required per group in heatmaps |
| `--format` | str   | png     | Output format for figures (tiff, svg, png, pdf)            |

### 🧪 Optional Arguments for Sequence QC Filtering

| Argument                | Type  | Default | Description                                      |
| ----------------------- | ----- | ------- | ------------------------------------------------ |
| `--qc_filter`           | bool  | false   | Use only QC-passing assemblies downstream        |
| `--min_total_length`    | int   | -       | Minimum assembly length required to pass QC      |
| `--max_contigs`         | int   | -       | Maximum contig count allowed to pass QC          |
| `--min_n50`             | int   | -       | Minimum N50 required to pass QC                  |
| `--min_gc`              | float | -       | Minimum GC percentage allowed to pass QC         |
| `--max_gc`              | float | -       | Maximum GC percentage allowed to pass QC         |
| `--max_ambiguous_bases` | int   | -       | Maximum ambiguous/gap bases allowed to pass QC   |
| `--checkm2_db`          | path  | -       | Optional CheckM2 database path                   |
| `--checkm2_auto_download_db` | bool | true | Download the CheckM2 database automatically when no database path is provided |
| `--checkm2_db_dir`      | path  | `<outdir>/databases/checkm2` | CheckM2 database download/cache directory |
| `--min_completeness`    | float | -       | Minimum CheckM2 completeness required to pass QC |
| `--max_contamination`   | float | -       | Maximum CheckM2 contamination allowed to pass QC |
| `--checkm2_lowmem`      | bool  | true    | Run CheckM2 in low-memory mode                   |
| `--stop_after_qc`       | bool  | false   | Stop after sequence QC and CheckM2               |
| `--run_gtdbtk`          | bool  | false   | Enable GTDB-Tk taxonomy QC                       |
| `--gtdbtk_data_path`    | path  | -       | Optional GTDB-Tk reference data path             |
| `--taxonomy_match_rank` | str   | genus   | Compare expected taxonomy at genus or species    |
| `--run_quast`           | bool  | false   | Enable QUAST assembly-structure QC               |
| `--run_ani`             | bool  | false   | Enable FastANI/skani pairwise ANI analysis       |
| `--ani_tool`            | str   | fastani | ANI engine: `fastani` or `skani`                 |
| `--ani_duplicate_threshold` | float | 99.9 | ANI threshold for near-duplicate clusters        |
| `--ani_species_threshold` | float | 95.0 | ANI warning threshold for species consistency    |
| `--run_mash`            | bool  | false   | Enable Mash sketch/distance pre-screening        |
| `--representative_only` | bool  | false   | Keep one representative per near-duplicate ANI cluster when filtering |
| `--export_panr2_inputs` | bool  | true    | Export standardized `panr2_inputs/` handoff directory |
| `--run_panr2_comprehensive` | bool | false | Run PanR2 integrated ABRicate NCBI/VFDB/PlasmidFinder, MobileElementFinder, IntegronFinder, and MLST |
| `--panr2_setup_abricate_db` | bool | true | Run `panr setup-db` before comprehensive analysis |
| `--panr2_abricate_dbs` | str | ncbi,vfdb,plasmidfinder | ABRicate databases used in comprehensive mode; add `isfinder` only if installed |
| `--panr2_run_defensefinder` | bool | false | Add DefenseFinder when a working installation is available |
| `--defensefinder_dir` | path | - | Existing DefenseFinder table directory to pass into PanR2 |
| `--panr2_min_identity` | float | 90 | Minimum identity threshold for PanR2 feature calls |
| `--panr2_plot_style` | str | publication | PanR2 plot preset: `publication`, `dashboard`, or `compact` |
| `--panr2_label_max_length` | int | 40 | Maximum feature label length in crowded plots |
| `--panr2_sample_map` | path | - | Optional `sample_id` to `Assembly Accession` map for external PanR2 table inputs; FetchM2 `sample_map.csv` is used automatically when present |
| `--run_mobsuite` | bool | false | Run MOB-suite and pass plasmid reconstruction/typing tables into PanR2 |
| `--mobsuite_dir` | path | - | Existing MOB-suite table directory to pass into PanR2 |
| `--run_genomad` | bool | false | Run geNomad and pass prophage/viral-region tables into PanR2 |
| `--prophage_dir` | path | - | Existing prophage/viral-region table directory to pass into PanR2 |
| `--genomad_db` | path | - | geNomad database directory required for `--run_genomad true` |
| `--run_organism_specific_typing` | bool | false | Run available organism-specific typing helpers |
| `--run_kleborate` | bool | false | Run Kleborate when available |
| `--kleborate_dir` | path | - | Existing Kleborate table directory to pass into PanR2 |
| `--run_kaptive` | bool | false | Run Kaptive when `--kaptive_db` is provided |
| `--kaptive_dir` | path | - | Existing Kaptive table directory to pass into PanR2 |
| `--kaptive_db` | path | - | Kaptive database path |
| `--run_ectyper` | bool | false | Run ECTyper when available |
| `--ectyper_dir` | path | - | Existing ECTyper table directory to pass into PanR2 |
| `--serotypefinder_dir` | path | - | Existing SerotypeFinder table directory to pass into PanR2 |
| `--sccmecfinder_dir` | path | - | Existing SCCmecFinder table directory to pass into PanR2 |

### 🔧 Other Options

| Argument    | Type | Default | Description                             |
| ----------- | ---- | ------- | --------------------------------------- |
| `--threads` | int  | 8       | Number of threads for CheckM2, GTDB-Tk, QUAST, ANI, and abricate |
| `--db`      | str  | ./db    | Directory containing abricate databases |
| `--help`    | flag | -       | Show help message and exit              |


## 📂 Output Structure

```bash
results/
└── <organism>/
    ├── abricate/            # Raw and summary AMR annotation results
    ├── figures/             # PNG, TIFF, and interactive HTML visualizations
    │   ├── heatmap/
    │   ├── mean_ARG/
    │   ├── html_files/
    │   ├── index.html       # Navigation page to generated reports
    │   └── Stat_analysis/
    ├── merged_output/       # Cleaned, joined resistance tables
    ├── metadata_output/     # FetchM2 assembly, annotation, and standardized metadata summary
    │   ├── fetchm2_clean.csv
    │   ├── fetchm2_clean.tsv
    │   ├── fetchm2_all_assemblies.csv
    │   ├── sample_map.csv
    │   ├── metadata_completeness.csv
    │   ├── metadata_bias_warning.txt
    │   ├── fetchm2_manifest.json
    │   ├── fetchm2_report.md
    │   ├── metadata_engine.txt
    │   ├── ncbi_clean.csv
    │   ├── ncbi_clean_unfiltered.csv # original metadata when --qc_filter true
    │   ├── ncbi_clean_qc_pass.csv    # metadata rows passing enabled QC checks
    │   └── ncbi_enriched.csv # ncbi_clean.csv plus sequence QC, CheckM2, and optional GTDB-Tk columns
    ├── sequence_qc/
    │   ├── assembly_stats.tsv # seqkit assembly statistics
    │   └── qc_decisions.tsv   # sequence, CheckM2, GTDB-Tk, and combined QC decisions
    ├── checkm2/
    │   └── quality_report.tsv # CheckM2 completeness and contamination report
    ├── gtdbtk/                # only when --run_gtdbtk true
    │   └── *.summary.tsv      # GTDB-Tk taxonomy classification summaries
    ├── quast/                 # only when --run_quast true
    │   └── analysis/assembly_qc.csv
    ├── ani/                   # only when --run_ani true
    │   └── analysis/
    │       ├── pairwise_ani_long.csv
    │       ├── ani_matrix.csv
    │       ├── closest_genome.csv
    │       ├── duplicate_clusters.csv
    │       └── ani_outliers.csv
    ├── mash/                  # only when --run_mash true
    │   └── analysis/
    ├── vfdb/                  # only when --run_panr2_comprehensive true
    ├── plasmidfinder/         # only when --run_panr2_comprehensive true
    ├── mobileelementfinder/   # only when --run_panr2_comprehensive true
    ├── integronfinder/        # only when --run_panr2_comprehensive true
    ├── mlst/                  # only when --run_panr2_comprehensive true
    ├── mobsuite/              # only when --run_mobsuite true
    ├── prophage/              # only when --run_genomad true or prophage tables are supplied
    ├── kleborate/             # only when organism typing is enabled/supplied
    ├── kaptive/
    ├── ectyper/
    ├── serotypefinder/
    ├── sccmecfinder/
    ├── cross_database/        # PanR2 cross-feature associations
    ├── temporal/              # PanR2 temporal trend summaries
    ├── report/                # PanR2 report, dashboard, citations
    │   └── index.html
    ├── qc/
    │   ├── qc_master_report.csv
    │   ├── qc_pass_samples.txt
    │   ├── qc_fail_samples.txt
    │   ├── qc_warning_samples.txt
    │   └── excluded_for_panr2.csv
    ├── panr2_inputs/          # standardized handoff bundle for PanR2
    │   ├── metadata/
    │   ├── metadata_analysis/
    │   ├── metadata_audit/
    │   ├── sequence/
    │   ├── amr/
    │   ├── ani/
    │   ├── assembly_qc/
    │   ├── qc/
    │   └── manifest/
    ├── sequence_filtered/     # pass-only FASTA files used when --qc_filter true
    └── sequence/            # Downloaded genome FASTA files        
pipeline_versions/
├── fetchm_env_versions.txt   # Python, FetchM2/legacy FetchM, PanR2, seqkit versions
├── checkm2_env_versions.txt  # CheckM2 version
├── gtdbtk_env_versions.txt   # GTDB-Tk version, only when --run_gtdbtk true
├── ani_env_versions.txt      # FastANI/skani versions, only when --run_ani true
├── quast_env_versions.txt    # QUAST version, only when --run_quast true
├── mash_env_versions.txt     # Mash version, only when --run_mash true
└── abricate_env_versions.txt # ABRicate and Perl versions
```

You can open `index.html` in a browser for easy navigation of visual outputs.

## PanResistome And PanR2 Responsibilities

PanResistome should run heavy tools, manage Conda environments/databases, capture versions, filter genomes, and export standardized tables. FetchM2 is the default metadata engine because it provides richer standardized host, source, environment, geography, year, disease, and metadata-audit fields than the legacy FetchM path. PanR2 should remain a lightweight comparative analysis and reporting tool that reads standardized outputs.

Every new PanResistome module should export PanR2-compatible records when possible. The formal schema is documented in [`docs/panr2_input_contract.md`](docs/panr2_input_contract.md).

---

## 🧼 Work Directory

By default, the `.nextflow` and `work/` directories are preserved for reproducibility. To remove intermediate files after a successful run:

```bash
nextflow clean -f
rm -rf work/
```

Only do this once you've verified the results.

---

## 📑 Example Report

You can view an example output here:

➡️ [**Interactive HTML Report**](https://tasnimul-arabi-anik.github.io/PanR2/)

Here’s an updated **📑 Example Report** section that includes the figures you mentioned and links to the interactive HTML dashboard:

---

## 📑 Example Report

You can explore an example output showcasing PanResistome's capabilities:

➡️ [**Interactive HTML Dashboard**](https://tasnimul-arabi-anik.github.io/PanR2/)

### Few Examples from index.html:

* ![Geographic Distribution of ARGs](figures/Geographic_distribution.png)
  *Figure: Regional prevalence of antimicrobial resistance genes across sampled isolates.*

* ![Heatmap of ARG Presence](figures/heatmap.png)
  *Figure: Heatmap showing ARG distribution across genomes.*

* ![Mean ARG Count](figures/mean_ARG.png)
  *Figure: Average number of detected ARGs per genome grouped by region.*

* ![Resistance Gene Frequency](figures/Resistance_gene_frequency.png)
  *Figure: Frequency of specific ARGs in the dataset.*

* ![Resistance Gene Variation](figures/Resistance_gene_variation.png)
  *Figure: Variation in resistance gene presence across samples.*

* ![Correlation Between ARGs](figures/Correlation_plot.png)
  *Figure: Correlation matrix of co-occurring ARGs.*

---

## 🧾 Citations

If you use this pipeline in your research, please cite the following tools:
Recommended citation entries:

---

* **Nextflow:** Di Tommaso *et al.* (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*. [https://doi.org/10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)
* **ABRicate:** Seemann T. ABRicate: Mass screening of contigs for antimicrobial and virulence genes. GitHub repository: [https://github.com/tseemann/abricate](https://github.com/tseemann/abricate)
* **FetchM2:** FetchM2 metadata standardization, audit, and sequence-download workflow. GitHub repository: [https://github.com/Tasnimul-Arabi-Anik/FetchM2](https://github.com/Tasnimul-Arabi-Anik/FetchM2)
* **PanR2:** Anik TA. *PanR2: Panresistome Analysis Tool*. DOI: 10.1101/2025.04.08.647722

---

## 💬 Contact

For suggestions, bug reports, or collaboration:

📧 Email: arabianik987@gmail.com
🌐 GitHub: [@Tasnimul-Arabi-Anik](https://github.com/Tasnimul-Arabi-Anik)

---

## 📜 License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---
