# PanResistome: Scalable Pipeline for Global Antimicrobial Resistance Analysis

Current pipeline version: `0.2.1-dev`

## Overview

**PanResistome** is a scalable, modular, and reproducible bioinformatics pipeline built using [Nextflow](https://www.nextflow.io/). It automates the end-to-end analysis of global antimicrobial resistance (AMR) patterns in bacterial populations using genome assemblies. The pipeline is designed for researchers working in microbial genomics, resistome surveillance, and public health, enabling large-scale comparative analysis of resistance gene profiles across time and geography.

PanResistome integrates several state-of-the-art tools including:

* [**FetchM**](https://github.com/Tasnimul-Arabi-Anik/FetchM): for fetching genome assemblies and standardized NCBI metadata
* [**CheckM2**](https://github.com/chklovski/CheckM2): for genome completeness and contamination quality assessment
* [**GTDB-Tk**](https://github.com/Ecogenomics/GTDBTk): for taxonomy classification and genus/species consistency checks
* [**ABRicate**](https://github.com/tseemann/abricate): for resistance gene annotation using curated ANCBI databases
* [**PanR2**](https://github.com/Tasnimul-Arabi-Anik/PanR2): for downstream statistical analysis and interactive visualization of resistome data

## Key Features

* 🔄 **Fully automated** end-to-end pipeline from genome download to visualization
* 🧬 **Panresistome analysis** using resistance gene profiling from ABRicate
* 📊 **Visualization-ready outputs** including heatmaps, barplots, boxplots, and interactive HTML figures
* 📈 **Statistical summaries** and correlation-based insights on resistance gene distribution
* 🌍 **Geospatial & temporal comparison** of AMR gene prevalence
* 💡 **Epidemic signal detection** by comparing ARG prevalence across time and location
* ⚙️ **Nextflow-based** for reproducibility, scalability, and cloud/HPC compatibility

---

## Workflow Overview

```
+-------------+   +-------------+   +---------+   +----------+   +-------+   +--------+
|   FetchM    |-->| Sequence QC |-->| CheckM2 |-->| ABRicate |-->| PanR2 |-->| Output |
+-------------+   +-------------+   +---------+   +----------+   +-------+   +--------+
| Assemblies  |   | Assembly    |   | Quality |   | AMR gene |   | Stats |   | Reports|
| & metadata  |   | stats       |   | metrics |   | calls    |   | plots |   | tables |
+-------------+   +-------------+   +---------+   +----------+   +-------+   +--------+

Optional: add `--run_gtdbtk true` to insert GTDB-Tk taxonomy matching between CheckM2 and ABRicate.
```

Each run also writes Conda environment version reports under `pipeline_versions/` so analyses can be traced back to the exact tool versions used.

Sequence QC filtering is optional. By default, the pipeline reports QC metrics but keeps all assemblies for downstream analysis. Add `--qc_filter true` with one or more thresholds to exclude failed assemblies from ABRicate, PanR2, and later tools.

CheckM2 requires its genome-quality database to be available in the run environment. If it is not configured globally, pass the database location with `--checkm2_db /path/to/checkm2_database`.

For a lighter validation run on modest hardware, use `--stop_after_qc true`. This runs FetchM, sequence QC, and CheckM2, then collects QC outputs without running GTDB-Tk, ABRicate, or PanR2.

GTDB-Tk is disabled by default because it is resource-intensive. If enabled, it requires its reference data to be available in the run environment. If `GTDBTK_DATA_PATH` is not configured globally, pass the location with `--gtdbtk_data_path /path/to/gtdbtk_data`. Taxonomy matching compares GTDB-Tk classification against the organism/species metadata at genus rank by default; use `--taxonomy_match_rank species` for stricter matching.

---

## 🚀 Getting Started

### ✅ Prerequisites

* [Nextflow](https://www.nextflow.io/)
* [Conda](https://docs.conda.io/en/latest/)
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
│   ├── fetchm.yml
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

### ⚙️ Optional Arguments for FetchM

| Argument    | Type   | Default | Description                                           |
| ----------- | ------ | ------- | ----------------------------------------------------- |
| `--checkm`  | float  | -       | Minimum CheckM completeness threshold (e.g., 90)      |
| `--ani`     | str    | all     | ANI filter status: OK, Inconclusive, Failed, or all   |
| `--sleep`   | float  | 0.5     | Time to wait between fetch requests (in seconds)      |
| `--host`    | str\[] | -       | Host species (e.g., "Homo sapiens", "Bos taurus")     |
| `--year`    | str\[] | -       | Filter by year or range (e.g., "2015" or "2015-2023") |
| `--country` | str\[] | -       | Filter by country (e.g., "Bangladesh", "USA")         |
| `--cont`    | str\[] | -       | Filter by continent (e.g., "Asia", "Africa")          |
| `--subcont` | str\[] | -       | Filter by subcontinent (e.g., "Southern Asia")        |

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
| `--min_completeness`    | float | -       | Minimum CheckM2 completeness required to pass QC |
| `--max_contamination`   | float | -       | Maximum CheckM2 contamination allowed to pass QC |
| `--checkm2_lowmem`      | bool  | true    | Run CheckM2 in low-memory mode                   |
| `--stop_after_qc`       | bool  | false   | Stop after sequence QC and CheckM2               |
| `--run_gtdbtk`          | bool  | false   | Enable GTDB-Tk taxonomy QC                       |
| `--gtdbtk_data_path`    | path  | -       | Optional GTDB-Tk reference data path             |
| `--taxonomy_match_rank` | str   | genus   | Compare expected taxonomy at genus or species    |

### 🔧 Other Options

| Argument    | Type | Default | Description                             |
| ----------- | ---- | ------- | --------------------------------------- |
| `--threads` | int  | 8       | Number of threads for CheckM2, GTDB-Tk, and abricate |
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
    ├── metadata_output/     # Assembly, annotation and metadata summary
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
    ├── sequence_filtered/     # pass-only FASTA files used when --qc_filter true
    └── sequence/            # Downloaded genome FASTA files        
pipeline_versions/
├── fetchm_env_versions.txt   # Python, FetchM, PanR2, seqkit versions
├── checkm2_env_versions.txt  # CheckM2 version
├── gtdbtk_env_versions.txt   # GTDB-Tk version, only when --run_gtdbtk true
└── abricate_env_versions.txt # ABRicate and Perl versions
```

You can open `index.html` in a browser for easy navigation of visual outputs.

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
Here are updated citation entries with proper formatting, pointing to the shared preprint for both FetchM and PanR2:

---

* **Nextflow:** Di Tommaso *et al.* (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*. [https://doi.org/10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)
* **ABRicate:** Seemann T. ABRicate: Mass screening of contigs for antimicrobial and virulence genes. GitHub repository: [https://github.com/tseemann/abricate](https://github.com/tseemann/abricate)
* **FetchM & PanR2:**
  Anik TA. *FetchM: Streamlining Genome and Metadata Integration for Microbial Comparative Genomics* (2025). Preprint available via ResearchGate/bioRxiv, DOI: 10.1101/2025.04.08.647722 ([researchgate.net](https://www.researchgate.net/publication/390754932_FetchM_Streamlining_Genome_and_Metadata_Integration_for_Microbial_Comparative_Genomics?utm_source=chatgpt.com))

---

## 💬 Contact

For suggestions, bug reports, or collaboration:

📧 Email: arabianik987@gmail.com
🌐 GitHub: [@Tasnimul-Arabi-Anik](https://github.com/Tasnimul-Arabi-Anik)

---

## 📜 License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---
