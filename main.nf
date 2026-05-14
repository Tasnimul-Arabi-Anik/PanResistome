#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def pipelineStartMillis = System.currentTimeMillis()

// Parameters
params.pipeline_version = '0.6.0'
params.checkm = null
params.ani = 'all'
params.sleep = 0.5
params.host = []
params.year = []
params.country = []
params.cont = []
params.subcont = []
params.metadata_engine = 'fetchm2'
params.fetchm2_offline = false
params.fetchm2_no_analysis = false
params.fetchm2_download = true
params.fetchm2_download_engine = 'native'
params.fetchm2_workers = 3
params.fetchm2_download_workers = 1
params.fetchm2_retries = 3
params.fetchm2_retry_delay = 5.0
params.fetchm2_keep_gz = false
params.fetchm2_max_genomes = null
params.fetchm2_keep_assembly_duplicates = false
params.sample_type = []
params.isolation_source = []
params.environment_medium = []
params.year_from = null
params.year_to = null

params.input = "test.tsv"
params.outdir = "results"
params.threads = 8
params.checkm2_threads = null
params.db = "$baseDir/db"
params.help = false
params.format = "png"
params.genep   = null
params.nseq    = null
params.local_samples = null
params.capture_versions = true
params.sequence_qc_engine = 'seqkit'
params.qc_filter = false
params.min_total_length = null
params.max_contigs = null
params.min_n50 = null
params.min_gc = null
params.max_gc = null
params.max_ambiguous_bases = null
params.checkm2_db = null
params.checkm2_auto_download_db = true
params.checkm2_db_dir = null
params.checkm2_db_download_retries = 3
params.min_completeness = null
params.max_contamination = null
params.run_checkm2 = true
params.checkm2_lowmem = true
params.stop_after_qc = false
params.run_gtdbtk = false
params.gtdbtk_data_path = null
params.taxonomy_match_rank = 'genus'
params.run_quast = false
params.run_ani = false
params.ani_tool = 'fastani'
params.ani_duplicate_threshold = 99.9
params.ani_species_threshold = 95.0
params.ani_large_run_strategy = 'auto'
params.ani_max_all_vs_all_genomes = 200
params.run_mash = false
params.representative_only = false
params.run_abricate = true
params.export_panr2_inputs = true
params.run_panr2_comprehensive = false
params.panr2_setup_abricate_db = true
params.panr2_update_abricate_db = true
params.panr2_abricate_dbs = 'ncbi,vfdb,plasmidfinder'
params.panr2_min_identity = 90
params.panr2_plot_style = 'publication'
params.panr2_label_max_length = 40
params.panr2_cross_database_max_features = 300
params.large_dataset = false
params.report_mode = 'publication'
params.max_features_heatmap = null
params.max_features_network = null
params.max_metadata_columns = null
params.top_n_features_per_database = null
params.skip_heavy_interactive_plots = false
params.panr2_force_tool_run = false
params.panr2_native_feature_runners = true
params.panr2_native_feature_runner_mode = 'serial'
params.panr2_run_mobileelementfinder = false
params.panr2_mobileelementfinder_allow_failure = true
params.panr2_run_defensefinder = false
params.panr2_sample_map = null
params.defensefinder_dir = null
params.run_isfinder = false
params.isfinder_db_fasta = null
params.isfinder_dir = null
params.isfinder_min_identity = 90
params.isfinder_min_coverage = 80
params.isfinder_database_version = null
params.analysis_profile = 'custom'
params.run_amrfinderplus = false
params.amrfinderplus_dir = null
params.amrfinderplus_organism = null
params.amrfinderplus_update_db = true
params.amrfinderplus_jobs = null
params.amrfinderplus_threads_per_sample = 1
params.amrfinderplus_reuse_existing = true
params.amrfinderplus_progress_every = 10
params.run_mobsuite = false
params.mobsuite_dir = null
params.mobsuite_db = null
params.mobsuite_db_dir = null
params.mobsuite_auto_init_db = true
params.mobsuite_auto_init_taxa = true
params.mobsuite_jobs = null
params.mobsuite_threads_per_sample = 1
params.mobsuite_reuse_existing = true
params.run_genomad = false
params.prophage_dir = null
params.genomad_db = null
params.genomad_db_dir = null
params.genomad_auto_download_db = true
params.genomad_use_host_env = false
params.genomad_splits = null
params.genomad_sensitivity = null
params.genomad_jobs = null
params.genomad_threads_per_sample = 1
params.genomad_reuse_existing = true
params.run_organism_specific_typing = false
params.run_kleborate = false
params.kleborate_dir = null
params.kleborate_jobs = null
params.kleborate_reuse_existing = true
params.run_kaptive = false
params.kaptive_dir = null
params.kaptive_db = null
params.run_ectyper = false
params.ectyper_dir = null
params.serotypefinder_dir = null
params.sccmecfinder_dir = null
params.core_feature_threshold = 0.95
params.rare_feature_threshold = 0.05
params.slurm_queue = null
params.slurm_account = null
params.slurm_cluster_options = ''
params.slurm_queue_size = 50


def paramList(value) {
    if (value == null || value == false) {
        return []
    }
    if (value instanceof List) {
        return value.findAll { it != null && it.toString().trim() }
    }
    return value.toString().split(',').collect { it.trim() }.findAll { it }
}

def shellQuote(value) {
    return "'" + value.toString().replace("'", "'\"'\"'") + "'"
}

def launchPath(value) {
    if (!value) {
        return ""
    }
    def text = value.toString()
    return text.startsWith("/") ? text : "${launchDir}/${text}"
}

def joinedOption(name, value) {
    def values = paramList(value)
    return values ? "${name} ${values.collect { shellQuote(it) }.join(' ')}" : ""
}

def yearRangeOptions(value) {
    def values = paramList(value)
    if (!values) {
        return ""
    }
    def options = []
    values.each { item ->
        def text = item.toString().trim()
        if (text.contains('-')) {
            def parts = text.split('-', 2).collect { it.trim() }
            if (parts[0]) {
                options << "--year-from ${shellQuote(parts[0])}"
            }
            if (parts.size() > 1 && parts[1]) {
                options << "--year-to ${shellQuote(parts[1])}"
            }
        } else if (text) {
            options << "--year-from ${shellQuote(text)}"
            options << "--year-to ${shellQuote(text)}"
        }
    }
    return options.join(' ')
}

def analysisProfile() {
    return (params.analysis_profile ?: 'custom').toString().trim().toLowerCase()
}

def effectiveStopAfterQc() {
    return params.stop_after_qc || analysisProfile() == 'qc_only'
}

def effectiveRunPanr2Comprehensive() {
    return params.run_panr2_comprehensive || analysisProfile() in ['amr_basic', 'amr_vp', 'amr_vp_mge', 'comprehensive']
}

def effectivePanr2Dbs() {
    def profile = analysisProfile()
    if (params.panr2_abricate_dbs != 'ncbi,vfdb,plasmidfinder') {
        return params.panr2_abricate_dbs
    }
    if (profile == 'amr_basic') {
        return 'ncbi'
    }
    return params.panr2_abricate_dbs
}

def abricateSetupCommand(panr2Dbs, sampleDir) {
    def command = "python ${baseDir}/scripts/setup_abricate_databases.py --dbs ${panr2Dbs} --out ${sampleDir}/panr2_inputs/manifest/abricate_database_setup_status.tsv"
    if (!params.panr2_setup_abricate_db) {
        command += " --check-only"
    }
    if (params.panr2_update_abricate_db) {
        command += " --update"
    }
    return command
}

def outdirDatabasePath(name) {
    return params.outdir.toString().startsWith("/") ? "${params.outdir}/databases/${name}" : "${launchDir}/${params.outdir}/databases/${name}"
}

def effectiveMobsuiteDbPath() {
    if (params.mobsuite_db) {
        return launchPath(params.mobsuite_db)
    }
    if (params.mobsuite_db_dir) {
        return launchPath(params.mobsuite_db_dir)
    }
    return outdirDatabasePath("mobsuite")
}

def effectiveGenomadDbPath() {
    if (params.genomad_db) {
        return launchPath(params.genomad_db)
    }
    if (params.genomad_db_dir) {
        return launchPath(params.genomad_db_dir)
    }
    return outdirDatabasePath("genomad")
}

def effectiveRunMobileElementFinder() {
    return params.panr2_run_mobileelementfinder
}

def effectiveRunIntegronFinder() {
    def profile = analysisProfile()
    return (params.run_panr2_comprehensive && profile == 'custom') || profile in ['amr_vp_mge', 'comprehensive']
}

def effectiveRunMlst() {
    def profile = analysisProfile()
    return (params.run_panr2_comprehensive && profile == 'custom') || profile == 'comprehensive'
}

def effectiveRunIsfinder() {
    def profile = analysisProfile()
    return params.run_isfinder || (params.isfinder_db_fasta && profile in ['amr_vp_mge', 'comprehensive'])
}

def truthyParam(value) {
    return value != null && value.toString().trim().toLowerCase() in ['true', '1', 'yes', 'y']
}

def effectiveReportMode() {
    def mode = (params.report_mode ?: '').toString().trim().toLowerCase()
    if (!mode || mode == 'publication') {
        return truthyParam(params.large_dataset) ? 'compact' : 'publication'
    }
    return mode
}

def effectiveMaxFeaturesHeatmap() {
    return params.max_features_heatmap ? params.max_features_heatmap as int : (truthyParam(params.large_dataset) ? 150 : 300)
}

def effectiveMaxFeaturesNetwork() {
    return params.max_features_network ? params.max_features_network as int : (truthyParam(params.large_dataset) ? 150 : params.panr2_cross_database_max_features as int)
}

def effectiveMaxMetadataColumns() {
    return params.max_metadata_columns ? params.max_metadata_columns as int : (truthyParam(params.large_dataset) ? 20 : 80)
}

def effectiveTopNFeaturesPerDatabase() {
    return params.top_n_features_per_database ? params.top_n_features_per_database as int : (truthyParam(params.large_dataset) ? 50 : 25)
}

def effectiveSkipHeavyInteractivePlots() {
    return truthyParam(params.skip_heavy_interactive_plots) || truthyParam(params.large_dataset)
}

def effectivePanr2PlotStyle() {
    if (truthyParam(params.large_dataset) && params.panr2_plot_style == 'publication') {
        return 'compact'
    }
    return params.panr2_plot_style
}


// Help message
def helpMessage() {
    log.info """
📦 PanResistome Pipeline v${params.pipeline_version}
    ------------------------
    A Nextflow pipeline for downloading genomic data, identifying resistance genes, and visualizing pangenome resistome profiles.

    ▶️ Usage:
      nextflow run main.nf --input <input.tsv> --outdir <output_dir> [options]

    ✅ Required arguments:
      --input            Input TSV file listing genome accessions
      --outdir           Output directory for results

    ⚙️ Optional arguments for FetchM2 metadata/download:
      --metadata_engine   Metadata engine: fetchm2 or legacy_fetchm [default: fetchm2]
      --checkm           Minimum CheckM completeness threshold (e.g. 90. Default: null)
      --ani              ANI filter status (Choices: OK, Inconclusive, Failed, all. Default: all)
      --sleep            Time to wait between fetch requests (default: 0.5s)
      --fetchm2_offline   Use FetchM2 offline metadata mode [default: false]
      --fetchm2_no_analysis Skip FetchM2 metadata analysis figures/tables [default: false]
      --fetchm2_download Download assemblies from FetchM2 metadata [default: true]
      --fetchm2_download_engine Sequence downloader: native or panresistome [default: native]
      --fetchm2_workers   FetchM2 BioSample fetch workers [default: 3]
      --fetchm2_download_workers FetchM2 sequence download workers [default: 1]
      --fetchm2_max_genomes Maximum genomes selected for FetchM2 sequence download
      --fetchm2_keep_assembly_duplicates Keep paired GCA/GCF rows in fetchm2_clean.csv [default: false]

        🧬 Instead of global resistance analysis, you may do specific analysis by providing: 
      --host             Host species (e.g. "Homo sapiens" "Bos taurus")
      --year             Filter by year or year range (e.g. "2015" or "2015-2023")
      --country          Country filter (e.g. "Bangladesh" "USA")
      --cont             Continent filter (e.g. "Asia", "Africa")
      --subcont          Subcontinent filter (e.g. "Southern Asia")
      --sample_type      FetchM2 Sample_Type_SD filter
      --isolation_source FetchM2 Isolation_Source_SD filter
      --environment_medium FetchM2 Environment_Medium_SD filter
      --year_from        FetchM2 minimum Collection_Year filter
      --year_to          FetchM2 maximum Collection_Year filter

    🧬 Optional arguments for PanR2:
      --genep            Minimum % gene presence to include in heatmap (float)
      --nseq             Minimum number of sequences per group in heatmaps (int)
      --format           Output format for plots (tiff, svg, png, pdf) [default: png]

    🧪 Sequence QC:
      After assemblies are downloaded from FetchM2 metadata, seqkit generates assembly stats and the pipeline writes
      metadata_output/ncbi_enriched.csv with sequence QC columns appended to ncbi_clean.csv.
      --qc_filter              Exclude failed assemblies from downstream tools [default: false]
      --min_total_length       Minimum assembly length required to pass QC
      --max_contigs            Maximum contig count allowed to pass QC
      --min_n50                Minimum N50 required to pass QC
      --min_gc                 Minimum GC percentage allowed to pass QC
      --max_gc                 Maximum GC percentage allowed to pass QC
      --max_ambiguous_bases    Maximum ambiguous/gap bases allowed to pass QC
      --checkm2_db             Optional CheckM2 database path
      --checkm2_auto_download_db Download CheckM2 database automatically when --checkm2_db is not provided [default: true]
      --checkm2_db_dir         Directory for automatic CheckM2 database download [default: <outdir>/databases/checkm2]
      --checkm2_db_download_retries Number of CheckM2 database download attempts [default: 3]
      --checkm2_threads        Threads for CheckM2 only. Defaults to min(--threads, 4) to reduce desktop RAM pressure.
      --min_completeness       Minimum CheckM2 completeness required to pass QC
      --max_contamination      Maximum CheckM2 contamination allowed to pass QC
      --checkm2_lowmem         Run CheckM2 in low-memory mode [default: true]
      --stop_after_qc          Stop after sequence QC and CheckM2 [default: false]
      --run_gtdbtk             Enable GTDB-Tk taxonomy QC [default: false]
      --gtdbtk_data_path       Optional GTDB-Tk reference data path
      --taxonomy_match_rank    Compare GTDB-Tk classification to metadata at genus or species [default: genus]
      --run_quast              Enable QUAST assembly-structure QC [default: false]
      --run_ani                Enable FastANI/skani pairwise ANI analysis [default: false]
      --ani_tool               ANI tool: fastani or skani [default: fastani]
      --ani_duplicate_threshold ANI threshold for near-duplicate clusters [default: 99.9]
      --ani_species_threshold  ANI warning threshold for species consistency [default: 95.0]
      --ani_large_run_strategy Large ANI strategy: auto, all, or skip. In auto mode, large-dataset runs above --ani_max_all_vs_all_genomes are skipped with an audit file [default: auto]
      --ani_max_all_vs_all_genomes Maximum genomes for automatic all-vs-all ANI in large-dataset mode [default: 200]
      --run_mash               Enable Mash sketch/distance pre-screen [default: false]
      --representative_only    Keep one representative per near-duplicate ANI cluster for PanR2 when --qc_filter true [default: false]
      --run_abricate           Run legacy ABRicate/PanR branch when PanR2 comprehensive mode is disabled [default: true]
      --export_panr2_inputs    Export standardized panr2_inputs handoff directory [default: true]
      --analysis_profile       Optional mode: custom, qc_only, amr_basic, amr_vp, amr_vp_mge, comprehensive [default: custom]
      --run_panr2_comprehensive Run comprehensive PanR2 analysis; PanResistome runs standard feature runners first when --panr2_native_feature_runners true [default: false]
      --panr2_run_defensefinder Add DefenseFinder to comprehensive PanR2 mode when a working installation is available [default: false]
      --defensefinder_dir      Existing DefenseFinder table directory to pass into PanR2
      --run_amrfinderplus      Run NCBI AMRFinderPlus and export standardized PanR2 feature tables [default: false]
      --amrfinderplus_dir      Existing AMRFinderPlus table directory to export under panr2_inputs/features
      --amrfinderplus_organism Optional AMRFinderPlus --organism value for mutation-aware calls
      --amrfinderplus_update_db Run amrfinder -u before AMRFinderPlus execution [default: true]
      --amrfinderplus_jobs     Parallel AMRFinderPlus sample jobs [default: --threads]
      --amrfinderplus_threads_per_sample Threads per AMRFinderPlus sample job [default: 1]
      --amrfinderplus_reuse_existing Reuse non-empty per-sample AMRFinderPlus TSVs on resume [default: true]
      --amrfinderplus_progress_every Print AMRFinderPlus progress every N completed samples [default: 10]
      --panr2_setup_abricate_db Run panr setup-db before comprehensive PanR2 analysis [default: true]
      --panr2_update_abricate_db Force-refresh requested ABRicate databases with abricate-get_db after setup when available [default: true]
      --panr2_abricate_dbs     ABRicate databases for PanR2 comprehensive mode [default: ncbi,vfdb,plasmidfinder; add isfinder only if installed]
      --panr2_min_identity     Minimum identity for PanR2 integrated feature analysis [default: 90]
      --panr2_plot_style       PanR2 plot style: publication, dashboard, compact [default: publication]
      --panr2_label_max_length Maximum feature-label length in PanR2 plots [default: 40]
      --large_dataset          Enable large-dataset report safeguards and compact defaults [default: false]
      --report_mode            PanR2 handoff report mode: compact, publication, exploratory [default: publication; compact when --large_dataset true]
      --max_features_heatmap   Maximum features retained in handoff presence/absence matrices [default: 300; 150 in large-dataset mode]
      --max_features_network   Maximum features used for cross-database co-occurrence/proximity summaries [default: --panr2_cross_database_max_features; 150 in large-dataset mode]
      --max_metadata_columns   Maximum metadata audit rows shown in handoff HTML pages [default: 80; 20 in large-dataset mode]
      --top_n_features_per_database
                              Number of top prevalent features per database summarized for report navigation [default: 25; 50 in large-dataset mode]
      --skip_heavy_interactive_plots
                              Mark heavy interactive plots as skipped/deprioritized in report controls [default: false; true in large-dataset mode]
      --panr2_native_feature_runners
                              Run ABRicate, IntegronFinder, MLST, and optional MobileElementFinder under PanResistome, then pass precomputed directories to PanR2 [default: true]
      --panr2_native_feature_runner_mode
                              Native feature-runner backend: serial or parallel. Keep serial as validated fallback; use parallel for ABRicate DBs plus per-assembly IntegronFinder/MLST [default: serial]
      --panr2_sample_map       Optional sample_id to Assembly Accession map for external PanR2 table inputs
      --panr2_run_mobileelementfinder
                              Run MobileElementFinder in the PanR2 feature-runner layer. Disabled by default because some valid assemblies trigger upstream parser failures.
      --panr2_mobileelementfinder_allow_failure
                              Keep MobileElementFinder failures nonfatal and write auditable header-only outputs [default: true]
      --run_isfinder           Run PanResistome ISfinder-compatible BLAST annotation [default: false]
      --isfinder_db_fasta      Authorized ISfinder nucleotide FASTA used to build a local BLAST database
      --isfinder_dir           Existing ISfinder-style result directory to pass into PanR2
      --isfinder_min_identity  Minimum ISfinder BLAST identity percentage [default: 90]
      --isfinder_min_coverage  Minimum ISfinder BLAST subject coverage percentage [default: 80]
      --run_mobsuite           Run MOB-suite plasmid reconstruction/typing and pass outputs to PanR2 [default: false]
      --mobsuite_dir           Existing MOB-suite table directory to pass into PanR2
      --mobsuite_db            Existing MOB-suite database directory for mob_recon --database_directory
      --mobsuite_db_dir        MOB-suite database cache/init directory when --mobsuite_db is not supplied [default: <outdir>/databases/mobsuite]
      --mobsuite_auto_init_db  Run mob_init for the MOB-suite cache directory when needed [default: true]
      --mobsuite_auto_init_taxa Initialize ETE taxa.sqlite for MOB-suite cache when needed [default: true]
      --mobsuite_jobs          Parallel MOB-suite sample jobs [default: min(--threads, 8)]
      --mobsuite_threads_per_sample Threads used by each MOB-suite sample job [default: 1]
      --mobsuite_reuse_existing Reuse existing non-empty per-sample MOB-suite outputs on resumed runs [default: true]
      --run_genomad            Run geNomad viral/prophage annotation and pass outputs to PanR2 [default: false]
      --prophage_dir           Existing prophage/viral-region table directory to pass into PanR2
      --genomad_db             geNomad database directory, required when --run_genomad true
      --genomad_db_dir         geNomad database download/cache directory when --genomad_db is not supplied [default: <outdir>/databases/genomad]
      --genomad_auto_download_db Download geNomad database into --genomad_db_dir when --run_genomad true and --genomad_db is not supplied [default: true]
      --genomad_use_host_env   Do not create the geNomad Conda env; use a prebuilt host/container geNomad executable [default: false]
      --genomad_splits         Split geNomad/MMseqs searches to reduce memory usage; higher is slower but safer [default: geNomad default]
      --genomad_sensitivity    geNomAD/MMseqs marker-search sensitivity; lower can reduce memory/time [default: geNomad default]
      --genomad_jobs           Parallel geNomad sample jobs inside GENOMAD_PROPHAGE [default: min(--threads, 2)]
      --genomad_threads_per_sample Threads used by each geNomad sample job [default: 1]
      --genomad_reuse_existing Reuse existing non-empty per-sample geNomad summaries on resumed/interrupted runs [default: true]
      --run_organism_specific_typing Run organism-specific typing helpers where applicable [default: false]
      --run_kleborate          Run Kleborate and pass outputs to PanR2 [default: false]
      --kleborate_dir          Existing Kleborate table directory to pass into PanR2
      --kleborate_jobs         Parallel Kleborate sample jobs [default: min(--threads, 8)]
      --kleborate_reuse_existing Reuse existing non-empty per-sample Kleborate outputs on resumed runs [default: true]
      --run_kaptive            Run Kaptive when --kaptive_db is provided [default: false]
      --kaptive_dir            Existing Kaptive table directory to pass into PanR2
      --kaptive_db             Kaptive database path
      --run_ectyper            Run ECTyper and pass outputs to PanR2 [default: false]
      --ectyper_dir            Existing ECTyper table directory to pass into PanR2
      --serotypefinder_dir     Existing SerotypeFinder table directory to pass into PanR2
      --sccmecfinder_dir       Existing SCCmecFinder table directory to pass into PanR2
      --container_image        Experimental container image used by docker/apptainer/singularity profiles
      --container_run_options  Extra runtime options passed to Docker/Apptainer/Singularity profiles
      --local_samples          Optional directory of prebuilt sample folders for offline tests
      --run_checkm2            Enable CheckM2 QC [default: true]

    🔧 Other options:
  --threads          Number of threads for CheckM2, GTDB-Tk, and abricate [default: 8]
                     CheckM2 is capped separately by --checkm2_threads unless explicitly raised.
      --db               Directory containing abricate databases [default: ./db]
      --help             Show this help message and exit

    Resource profiles can be combined with -profile conda,mamba:
      lowmem             threads=4, checkm2_threads=1, serial native runners
      desktop_parallel   threads=16, checkm2_threads=2, parallel native runners
      workstation        threads=16, checkm2_threads=4, parallel native runners
      docker             experimental; pass --container_image and validate first
      apptainer          experimental; pass --container_image and validate first
      singularity        experimental; pass --container_image and validate first

    Example:
       nextflow run main.nf --input test_small.tsv --outdir results_small -profile conda --threads 8 
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

// Validate database directory
if (!file(params.db).exists()) {
    log.error "Database directory does not exist: ${params.db}"
    exit 1
}

// Process 1: Capture versions from the FetchM2/seqkit environment
process FETCHM_ENV_VERSIONS {
    conda 'envs/fetchm.yaml'
    publishDir "${params.outdir}/pipeline_versions", mode: 'copy'

    output:
    path "fetchm_env_versions.txt", emit: fetchm_versions

    script:
    """
    {
        echo "[fetchm_env]"
        echo "pipeline_version=${params.pipeline_version}"
        date -u +"run_timestamp_utc=%Y-%m-%dT%H:%M:%SZ"
        python --version
        seqkit version || true
        fetchm2 --version || true
        fetchM --version || true
        python - <<'PY'
from importlib.metadata import PackageNotFoundError, version

for package in ("fetchm2", "fetchM"):
    try:
        print(f"{package}=={version(package)}")
    except PackageNotFoundError:
        print(f"{package}==NOT_FOUND")
PY
    } > fetchm_env_versions.txt
    """
}

// Process 2: Capture versions from the ABRicate environment
process ABRICATE_ENV_VERSIONS {
    conda 'envs/abricate.yaml'
    publishDir "${params.outdir}/pipeline_versions", mode: 'copy'

    output:
    path "abricate_env_versions.txt", emit: abricate_versions

    script:
    """
    {
        echo "[abricate_env]"
        echo "pipeline_version=${params.pipeline_version}"
        date -u +"run_timestamp_utc=%Y-%m-%dT%H:%M:%SZ"
        abricate --version || true
        perl -e 'print "perl==" . \$^V . "\\n"'
    } > abricate_env_versions.txt
    """
}

process AMRFINDERPLUS_ENV_VERSIONS {
    conda 'envs/amrfinderplus.yaml'
    publishDir "${params.outdir}/pipeline_versions", mode: 'copy'

    output:
    path "amrfinderplus_env_versions.txt", emit: amrfinderplus_versions

    script:
    """
    {
        echo "[amrfinderplus_env]"
        echo "pipeline_version=${params.pipeline_version}"
        date -u +"run_timestamp_utc=%Y-%m-%dT%H:%M:%SZ"
        amrfinder --version || true
        amrfinder -l || true
    } > amrfinderplus_env_versions.txt
    """
}

process PANR2_COMPREHENSIVE_ENV_VERSIONS {
    conda 'envs/panr2_comprehensive.yaml'
    publishDir "${params.outdir}/pipeline_versions", mode: 'copy'

    output:
    path "panr2_comprehensive_env_versions.txt", emit: panr2_comprehensive_versions

    script:
    """
    {
        echo "[panr2_comprehensive_env]"
        python --version
        panr --version || true
        panr doctor || true
        abricate --version || true
        abricate --list || true
        integron_finder --version || true
        mlst --version || true
        defense-finder --version || true
        mefinder --version || true
    } > panr2_comprehensive_env_versions.txt
    """
}

// Process 3: Capture versions from the CheckM2 environment
process CHECKM2_ENV_VERSIONS {
    conda 'envs/checkm2.yaml'
    publishDir "${params.outdir}/pipeline_versions", mode: 'copy'

    output:
    path "checkm2_env_versions.txt", emit: checkm2_versions

    script:
    """
    {
        echo "[checkm2_env]"
        echo "pipeline_version=${params.pipeline_version}"
        date -u +"run_timestamp_utc=%Y-%m-%dT%H:%M:%SZ"
        checkm2 --version || true
    } > checkm2_env_versions.txt
    """
}

// Process 4: Capture versions from the GTDB-Tk environment
process GTDBTK_ENV_VERSIONS {
    conda 'envs/gtdbtk.yaml'
    publishDir "${params.outdir}/pipeline_versions", mode: 'copy'

    output:
    path "gtdbtk_env_versions.txt", emit: gtdbtk_versions

    script:
    """
    {
        echo "[gtdbtk_env]"
        echo "pipeline_version=${params.pipeline_version}"
        date -u +"run_timestamp_utc=%Y-%m-%dT%H:%M:%SZ"
        gtdbtk --version || true
    } > gtdbtk_env_versions.txt
    """
}

process ANI_ENV_VERSIONS {
    conda 'envs/ani.yaml'
    publishDir "${params.outdir}/pipeline_versions", mode: 'copy'

    output:
    path "ani_env_versions.txt", emit: ani_versions

    script:
    """
    {
        echo "[ani_env]"
        fastANI --version || true
        skani --version || true
    } > ani_env_versions.txt
    """
}

process QUAST_ENV_VERSIONS {
    conda 'envs/quast.yaml'
    publishDir "${params.outdir}/pipeline_versions", mode: 'copy'

    output:
    path "quast_env_versions.txt", emit: quast_versions

    script:
    """
    {
        echo "[quast_env]"
        quast.py --version || true
    } > quast_env_versions.txt
    """
}

process MASH_ENV_VERSIONS {
    conda 'envs/mash.yaml'
    publishDir "${params.outdir}/pipeline_versions", mode: 'copy'

    output:
    path "mash_env_versions.txt", emit: mash_versions

    script:
    """
    {
        echo "[mash_env]"
        mash --version || true
    } > mash_env_versions.txt
    """
}

process MOBSUITE_ENV_VERSIONS {
    conda 'envs/mobsuite.yaml'
    publishDir "${params.outdir}/pipeline_versions", mode: 'copy'

    output:
    path "mobsuite_env_versions.txt", emit: mobsuite_versions

    script:
    """
    {
        echo "[mobsuite_env]"
        echo "pipeline_version=${params.pipeline_version}"
        date -u +"run_timestamp_utc=%Y-%m-%dT%H:%M:%SZ"
        mob_recon --version || true
    } > mobsuite_env_versions.txt
    """
}

process GENOMAD_ENV_VERSIONS {
    conda 'envs/genomad.yaml'
    publishDir "${params.outdir}/pipeline_versions", mode: 'copy'

    output:
    path "genomad_env_versions.txt", emit: genomad_versions

    script:
    """
    {
        echo "[genomad_env]"
        echo "pipeline_version=${params.pipeline_version}"
        date -u +"run_timestamp_utc=%Y-%m-%dT%H:%M:%SZ"
        genomad --version || true
    } > genomad_env_versions.txt
    """
}

process ORGANISM_TYPING_ENV_VERSIONS {
    conda 'envs/organism_typing.yaml'
    publishDir "${params.outdir}/pipeline_versions", mode: 'copy'

    output:
    path "organism_typing_env_versions.txt", emit: organism_typing_versions

    script:
    """
    {
        echo "[organism_typing_env]"
        echo "pipeline_version=${params.pipeline_version}"
        date -u +"run_timestamp_utc=%Y-%m-%dT%H:%M:%SZ"
        kleborate --version || true
        kaptive --version || true
        ectyper --version || true
    } > organism_typing_env_versions.txt
    """
}

// Process 5: Run FetchM2 by default, with a reversible legacy FetchM mode.
process FETCHM {
    conda 'envs/fetchm.yaml'
    
    input:
    path input_file
    
    output:
    path "fetchm_results", emit: fetchm_results
    
    script:
    def fetchm2Ani = paramList(params.ani) ?: ['all']
    def fetchm2Args = [
        "--input ${input_file}",
        "--outdir fetchm_results",
        "--ani ${fetchm2Ani.collect { shellQuote(it) }.join(' ')}",
        "--workers ${params.fetchm2_workers}",
        "--sleep ${params.sleep}",
    ]
    if (params.checkm) {
        fetchm2Args << "--checkm ${params.checkm}"
    }
    if (params.fetchm2_offline) {
        fetchm2Args << "--offline"
    }
    if (params.fetchm2_no_analysis) {
        fetchm2Args << "--no-analysis"
    }
    if (params.fetchm2_keep_assembly_duplicates) {
        fetchm2Args << "--keep-assembly-duplicates"
    }
    def explicitYearOptions = []
    if (params.year_from) {
        explicitYearOptions << "--year-from ${shellQuote(params.year_from)}"
    }
    if (params.year_to) {
        explicitYearOptions << "--year-to ${shellQuote(params.year_to)}"
    }
    def filterArgs = [
        joinedOption("--host", params.host),
        joinedOption("--country", params.country),
        joinedOption("--continent", params.cont),
        joinedOption("--subcontinent", params.subcont),
        joinedOption("--sample-type", params.sample_type),
        joinedOption("--isolation-source", params.isolation_source),
        joinedOption("--environment-medium", params.environment_medium),
        explicitYearOptions ? explicitYearOptions.join(' ') : yearRangeOptions(params.year),
    ].findAll { it }
    fetchm2Args.addAll(filterArgs)
    def nativeDownloadArgs = [
        "--input fetchm_results/metadata_output/fetchm2_clean.csv",
        "--outdir fetchm_results/sequence",
        "--download-workers ${params.fetchm2_download_workers}",
        "--retries ${params.fetchm2_retries}",
        "--retry-delay ${params.fetchm2_retry_delay}",
        params.fetchm2_max_genomes ? "--max-genomes ${params.fetchm2_max_genomes}" : "",
        params.fetchm2_keep_gz ? "--keep-gz" : "",
    ].findAll { it }
    nativeDownloadArgs.addAll(filterArgs)
    def panresistomeDownloadArgs = [
        "--input fetchm_results/metadata_output/fetchm2_clean.csv",
        "--outdir fetchm_results/sequence",
        "--workers ${params.fetchm2_download_workers}",
        "--retries ${params.fetchm2_retries}",
        "--retry-delay ${params.fetchm2_retry_delay}",
        params.fetchm2_max_genomes ? "--max-genomes ${params.fetchm2_max_genomes}" : "",
        params.fetchm2_keep_gz ? "--keep-gz" : "",
    ].findAll { it }
    panresistomeDownloadArgs.addAll(filterArgs)
    def legacyArgs = [
        "--input ${input_file}",
        "--outdir fetchm_results/",
        params.checkm ? "--checkm ${params.checkm}" : "",
        "--ani ${params.ani}",
        "--sleep ${params.sleep}",
        "--seq",
        joinedOption("--host", params.host),
        joinedOption("--year", params.year),
        joinedOption("--country", params.country),
        joinedOption("--cont", params.cont),
        joinedOption("--subcont", params.subcont),
    ].findAll { it }
    """
    export MPLCONFIGDIR="\$(pwd)/.matplotlib"
    mkdir -p "\${MPLCONFIGDIR}"

    if [ "${params.metadata_engine}" = "fetchm2" ]; then
        fetchm2 metadata ${fetchm2Args.join(' ')}
        if [ "${params.fetchm2_download}" = "true" ]; then
            if [ "${params.fetchm2_download_engine}" = "native" ]; then
                fetchm2 seq ${nativeDownloadArgs.join(' ')} || {
                    echo "Warning: native fetchm2 seq failed; falling back to PanResistome downloader" >&2
                    rm -rf fetchm_results/sequence
                    python ${baseDir}/scripts/download_fetchm2_sequences.py ${panresistomeDownloadArgs.join(' ')}
                }
            elif [ "${params.fetchm2_download_engine}" = "panresistome" ]; then
                python ${baseDir}/scripts/download_fetchm2_sequences.py ${panresistomeDownloadArgs.join(' ')}
            else
                echo "Unsupported fetchm2_download_engine: ${params.fetchm2_download_engine}" >&2
                exit 1
            fi
        fi
        python ${baseDir}/scripts/normalize_fetchm2_output.py --results-dir fetchm_results
    elif [ "${params.metadata_engine}" = "legacy_fetchm" ] || [ "${params.metadata_engine}" = "fetchm" ]; then
        fetchM ${legacyArgs.join(' ')}
    else
        echo "Unsupported metadata engine: ${params.metadata_engine}" >&2
        exit 1
    fi
    """
}

// Process 6: Generate assembly QC stats and enrich metadata
process SEQUENCE_QC {
    conda 'envs/fetchm.yaml'
    stageInMode { params.local_samples ? 'copy' : 'symlink' }

    input:
    path sample_dir

    output:
    path "${sample_dir}", emit: qc_results

    script:
    """
    mkdir -p ${sample_dir}/sequence_qc
    mkdir -p ${sample_dir}/metadata_output
    rm -rf ${sample_dir}/sequence_filtered
    mkdir -p ${sample_dir}/sequence_filtered

    if [ "${params.sequence_qc_engine}" = "python" ] && [ -d "${sample_dir}/sequence" ] && [ -n "\$(find ${sample_dir}/sequence -name "*.fna" -print -quit)" ]; then
        python - <<'PY'
from pathlib import Path

sample_dir = Path("${sample_dir}")
stats_path = sample_dir / "sequence_qc" / "assembly_stats.tsv"

def fasta_lengths(path):
    lengths = []
    current = 0
    gc = 0
    gaps = 0
    total = 0
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
                continue
            seq = line.upper()
            current += len(seq)
            total += len(seq)
            gc += seq.count("G") + seq.count("C")
            gaps += sum(1 for base in seq if base not in {"A", "C", "G", "T"})
    if current:
        lengths.append(current)
    return lengths, total, gc, gaps

def n50(lengths):
    if not lengths:
        return 0
    half = sum(lengths) / 2
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= half:
            return length
    return 0

header = [
    "file", "format", "type", "num_seqs", "sum_len", "min_len", "avg_len",
    "max_len", "Q1", "Q2", "Q3", "sum_gap", "N50", "Q20(%)", "Q30(%)", "GC(%)",
]
with stats_path.open("w") as out:
    out.write("\\t".join(header) + "\\n")
    for path in sorted((sample_dir / "sequence").glob("*.fna")):
        lengths, total, gc, gaps = fasta_lengths(path)
        row = {
            "file": str(path),
            "format": "FASTA",
            "type": "DNA",
            "num_seqs": len(lengths),
            "sum_len": total,
            "min_len": min(lengths) if lengths else 0,
            "avg_len": f"{(total / len(lengths)):.1f}" if lengths else "0",
            "max_len": max(lengths) if lengths else 0,
            "Q1": "",
            "Q2": "",
            "Q3": "",
            "sum_gap": gaps,
            "N50": n50(lengths),
            "Q20(%)": "",
            "Q30(%)": "",
            "GC(%)": f"{((gc / total) * 100):.2f}" if total else "0",
        }
        out.write("\\t".join(str(row[field]) for field in header) + "\\n")
PY
    elif [ -d "${sample_dir}/sequence" ] && [ -n "\$(find ${sample_dir}/sequence -name "*.fna" -print -quit)" ]; then
        seqkit stats -a -T ${sample_dir}/sequence/*.fna > ${sample_dir}/sequence_qc/assembly_stats.tsv
    else
        echo "Warning: No .fna files found in ${sample_dir}/sequence/ for sequence QC" >&2
        printf '%b\\n' "file\\tformat\\ttype\\tnum_seqs\\tsum_len\\tmin_len\\tavg_len\\tmax_len\\tQ1\\tQ2\\tQ3\\tsum_gap\\tN50\\tQ20(%)\\tQ30(%)\\tGC(%)" > ${sample_dir}/sequence_qc/assembly_stats.tsv
    fi

    python - <<'PY'
import csv
from pathlib import Path

sample_dir = Path("${sample_dir}")
metadata_dir = sample_dir / "metadata_output"
stats_path = sample_dir / "sequence_qc" / "assembly_stats.tsv"
clean_path = metadata_dir / "ncbi_clean.csv"
enriched_path = metadata_dir / "ncbi_enriched.csv"
pass_path = metadata_dir / "ncbi_clean_qc_pass.csv"
unfiltered_path = metadata_dir / "ncbi_clean_unfiltered.csv"
decision_path = sample_dir / "sequence_qc" / "qc_decisions.tsv"
filtered_sequence_dir = sample_dir / "sequence_filtered"

qc_filter = str("${params.qc_filter}").strip().lower() in {"true", "1", "yes", "y"}
thresholds = {
    "min_total_length": "${params.min_total_length}",
    "max_contigs": "${params.max_contigs}",
    "min_n50": "${params.min_n50}",
    "min_gc": "${params.min_gc}",
    "max_gc": "${params.max_gc}",
    "max_ambiguous_bases": "${params.max_ambiguous_bases}",
}

def parse_float(value):
    value = str(value or "").strip()
    if value in {"", "null", "None"}:
        return None
    return float(value.replace(",", ""))

thresholds = {key: parse_float(value) for key, value in thresholds.items()}

def normalize(value):
    return str(value or "").strip().lower().replace("-", "_").replace(".", "_")

def key_candidates(value):
    value = str(value or "").strip()
    if not value:
        return set()
    stem = Path(value).stem
    candidates = {value, stem}
    if stem.endswith("_genomic"):
        candidates.add(stem[:-8])
    parts = stem.split("_")
    if len(parts) >= 2 and parts[0] in {"GCF", "GCA"}:
        candidates.add("_".join(parts[:2]))
    return {normalize(candidate) for candidate in candidates if normalize(candidate)}

def numeric(value):
    try:
        return float(str(value or "").replace(",", ""))
    except ValueError:
        return None

def qc_decision(stats):
    reasons = []
    checks = [
        ("min_total_length", numeric(stats.get("sequence_total_length")), ">="),
        ("max_contigs", numeric(stats.get("sequence_num_contigs")), "<="),
        ("min_n50", numeric(stats.get("sequence_n50")), ">="),
        ("min_gc", numeric(stats.get("sequence_gc_percent")), ">="),
        ("max_gc", numeric(stats.get("sequence_gc_percent")), "<="),
        ("max_ambiguous_bases", numeric(stats.get("sequence_ambiguous_bases")), "<="),
    ]
    if not stats.get("sequence_num_contigs"):
        reasons.append("NO_SEQUENCE")
    for threshold_name, observed, operator in checks:
        threshold = thresholds[threshold_name]
        if threshold is None:
            continue
        if observed is None:
            reasons.append(f"{threshold_name}:missing")
        elif operator == ">=" and observed < threshold:
            reasons.append(f"{threshold_name}:{observed:g}<{threshold:g}")
        elif operator == "<=" and observed > threshold:
            reasons.append(f"{threshold_name}:{observed:g}>{threshold:g}")
    return ("PASS", "") if not reasons else ("FAIL", ";".join(reasons))

def read_stats(path):
    rows = []
    if not path.exists():
        return rows
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\\t")
        for row in reader:
            seq_path = Path(row.get("file", ""))
            accession = seq_path.stem
            stats = {
                "sequence_file": seq_path.name,
                "sequence_accession": accession,
                "sequence_num_contigs": row.get("num_seqs", ""),
                "sequence_total_length": row.get("sum_len", ""),
                "sequence_min_contig": row.get("min_len", ""),
                "sequence_avg_contig": row.get("avg_len", ""),
                "sequence_max_contig": row.get("max_len", ""),
                "sequence_n50": row.get("N50", ""),
                "sequence_gc_percent": row.get("GC(%)", ""),
                "sequence_ambiguous_bases": row.get("sum_gap", ""),
            }
            status, reason = qc_decision(stats)
            stats["sequence_qc_status"] = status
            stats["sequence_qc_fail_reasons"] = reason
            rows.append(stats)
    return rows

stats_rows = read_stats(stats_path)
stats_by_key = {}
for stats in stats_rows:
    keys = key_candidates(stats["sequence_accession"]) | key_candidates(stats["sequence_file"])
    for key in keys:
        if key:
            stats_by_key[key] = stats

metadata_keys = [
    "assembly_accession", "Assembly Accession", "assembly", "accession",
    "genome_accession", "Genome accession", "Assembly", "sample_id",
]

added_fields = [
    "sequence_file", "sequence_accession", "sequence_num_contigs",
    "sequence_total_length", "sequence_min_contig", "sequence_avg_contig",
    "sequence_max_contig", "sequence_n50", "sequence_gc_percent",
    "sequence_ambiguous_bases", "sequence_qc_status", "sequence_qc_fail_reasons",
]

if clean_path.exists():
    with clean_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        output_fields = fieldnames + [field for field in added_fields if field not in fieldnames]
        rows = []
        for row in reader:
            match = None
            for key in metadata_keys:
                if key in row:
                    for candidate in key_candidates(row[key]):
                        if candidate in stats_by_key:
                            match = stats_by_key[candidate]
                            break
                if match:
                    break
            enriched = dict(row)
            for field in added_fields:
                if match:
                    enriched[field] = match.get(field, "")
                elif field == "sequence_qc_status":
                    enriched[field] = "FAIL"
                elif field == "sequence_qc_fail_reasons":
                    enriched[field] = "NO_SEQUENCE"
                else:
                    enriched[field] = ""
            rows.append(enriched)
else:
    output_fields = added_fields
    rows = stats_rows

pass_rows = [row for row in rows if row.get("sequence_qc_status") == "PASS"]

with enriched_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=output_fields)
    writer.writeheader()
    writer.writerows(rows)

with pass_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=output_fields)
    writer.writeheader()
    writer.writerows(pass_rows)

with decision_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["sequence_file", "sequence_qc_status", "sequence_qc_fail_reasons"], delimiter="\\t")
    writer.writeheader()
    for row in stats_rows:
        writer.writerow({
            "sequence_file": row.get("sequence_file", ""),
            "sequence_qc_status": row.get("sequence_qc_status", ""),
            "sequence_qc_fail_reasons": row.get("sequence_qc_fail_reasons", ""),
        })

for row in stats_rows:
    if row.get("sequence_qc_status") != "PASS":
        continue
    source = sample_dir / "sequence" / row["sequence_file"]
    target = filtered_sequence_dir / row["sequence_file"]
    if source.exists():
        target.write_bytes(source.read_bytes())

if qc_filter and clean_path.exists():
    if not unfiltered_path.exists():
        unfiltered_path.write_bytes(clean_path.read_bytes())
    clean_path.write_bytes(pass_path.read_bytes())
PY
    """
}

// Process 7: Run CheckM2 and add genome quality metrics
process CHECKM2_QC {
    conda 'envs/checkm2.yaml'

    input:
    path sample_dir

    output:
    path "${sample_dir}", emit: checkm2_results

    script:
    def checkm2DbPath = params.checkm2_db ? launchPath(params.checkm2_db) : ""
    def checkm2_db_arg = checkm2DbPath ? "--database_path ${checkm2DbPath}" : ""
    def checkm2_lowmem_arg = params.checkm2_lowmem ? "--lowmem" : ""
    def checkm2Threads = params.checkm2_threads ? (params.checkm2_threads as int) : Math.min((params.threads as int), 4)
    def checkm2DownloadDir = params.checkm2_db_dir ? launchPath(params.checkm2_db_dir) : (params.outdir.toString().startsWith("/") ? "${params.outdir}/databases/checkm2" : "${launchDir}/${params.outdir}/databases/checkm2")
    """
    mkdir -p ${sample_dir}/checkm2
    checkm2_db_arg="${checkm2_db_arg}"
    echo "CheckM2 thread cap: ${checkm2Threads} thread(s). General --threads remains ${params.threads} for other stages."
    if [ "${checkm2Threads}" -ge 8 ]; then
        echo "Warning: CheckM2 is running with ${checkm2Threads} threads. On 16-32 GB desktops, use --checkm2_threads 2 or 4 to reduce DIAMOND memory pressure." >&2
    fi

    if [ -z "\${checkm2_db_arg}" ] && [ "${params.checkm2_auto_download_db}" = "true" ]; then
        mkdir -p "${checkm2DownloadDir}"
        existing_db=\$(find "${checkm2DownloadDir}" -name "*.dmnd" -type f -print -quit 2>/dev/null || true)
        if [ -z "\${existing_db}" ]; then
            echo "No CheckM2 database provided; downloading to ${checkm2DownloadDir}"
            attempt=1
            max_attempts="${params.checkm2_db_download_retries}"
            while true; do
                echo "CheckM2 database download attempt \${attempt}/\${max_attempts}"
                if checkm2 database --download --path "${checkm2DownloadDir}" --no_write_json_db; then
                    break
                fi
                if [ "\${attempt}" -ge "\${max_attempts}" ]; then
                    echo "CheckM2 database download failed after \${max_attempts} attempts" >&2
                    exit 1
                fi
                attempt=\$((attempt + 1))
                sleep \$((attempt * 30))
            done
            existing_db=\$(find "${checkm2DownloadDir}" -name "*.dmnd" -type f -print -quit 2>/dev/null || true)
        fi
        if [ -n "\${existing_db}" ]; then
            checkm2_db_arg="--database_path \${existing_db}"
            echo "Using CheckM2 database: \${existing_db}"
        else
            echo "CheckM2 database download did not produce a .dmnd file under ${checkm2DownloadDir}" >&2
            exit 1
        fi
    fi

    sequence_dir="${sample_dir}/sequence"
    if [ "${params.qc_filter}" = "true" ] && [ -d "${sample_dir}/sequence_filtered" ]; then
        sequence_dir="${sample_dir}/sequence_filtered"
    fi

    if [ -d "\${sequence_dir}" ] && [ -n "\$(find \${sequence_dir} -name "*.fna" -print -quit)" ]; then
        checkm2 predict \\
            --threads ${checkm2Threads} \\
            --input \${sequence_dir} \\
            --output-directory ${sample_dir}/checkm2 \\
            -x fna --force ${checkm2_lowmem_arg} \${checkm2_db_arg}
    else
        echo "Warning: No .fna files found in \${sequence_dir}/ for CheckM2" >&2
        printf "Name\\tCompleteness\\tContamination\\n" > ${sample_dir}/checkm2/quality_report.tsv
    fi

    python - <<'PY'
import csv
import shutil
from pathlib import Path

sample_dir = Path("${sample_dir}")
metadata_dir = sample_dir / "metadata_output"
checkm2_path = sample_dir / "checkm2" / "quality_report.tsv"
clean_path = metadata_dir / "ncbi_clean.csv"
enriched_path = metadata_dir / "ncbi_enriched.csv"
pass_path = metadata_dir / "ncbi_clean_qc_pass.csv"
unfiltered_path = metadata_dir / "ncbi_clean_unfiltered.csv"
decision_path = sample_dir / "sequence_qc" / "qc_decisions.tsv"
filtered_sequence_dir = sample_dir / "sequence_filtered"
source_sequence_dir = sample_dir / "sequence"

qc_filter = str("${params.qc_filter}").strip().lower() in {"true", "1", "yes", "y"}
thresholds = {
    "min_completeness": "${params.min_completeness}",
    "max_contamination": "${params.max_contamination}",
}

def parse_float(value):
    value = str(value or "").strip()
    if value in {"", "null", "None"}:
        return None
    return float(value.replace(",", ""))

thresholds = {key: parse_float(value) for key, value in thresholds.items()}

def normalize(value):
    return str(value or "").strip().lower().replace("-", "_").replace(".", "_")

def key_candidates(value):
    value = str(value or "").strip()
    if not value:
        return set()
    stem = Path(value).stem
    candidates = {value, stem}
    if stem.endswith("_genomic"):
        candidates.add(stem[:-8])
    parts = stem.split("_")
    if len(parts) >= 2 and parts[0] in {"GCF", "GCA"}:
        candidates.add("_".join(parts[:2]))
    return {normalize(candidate) for candidate in candidates if normalize(candidate)}

def numeric(value):
    try:
        return float(str(value or "").replace(",", ""))
    except ValueError:
        return None

def first(row, names):
    for name in names:
        if name in row:
            return row.get(name, "")
    return ""

def checkm2_decision(stats):
    reasons = []
    completeness = numeric(stats.get("checkm2_completeness"))
    contamination = numeric(stats.get("checkm2_contamination"))
    if not stats:
        reasons.append("NO_CHECKM2_RESULT")
    if thresholds["min_completeness"] is not None:
        if completeness is None:
            reasons.append("min_completeness:missing")
        elif completeness < thresholds["min_completeness"]:
            reasons.append(f"min_completeness:{completeness:g}<{thresholds['min_completeness']:g}")
    if thresholds["max_contamination"] is not None:
        if contamination is None:
            reasons.append("max_contamination:missing")
        elif contamination > thresholds["max_contamination"]:
            reasons.append(f"max_contamination:{contamination:g}>{thresholds['max_contamination']:g}")
    return ("PASS", "") if not reasons else ("FAIL", ";".join(reasons))

def read_checkm2(path):
    rows = {}
    if not path.exists():
        return rows
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\\t")
        for row in reader:
            name = first(row, ["Name", "Bin Id", "Bin_ID", "Genome", "genome"])
            stats = {
                "checkm2_name": name,
                "checkm2_completeness": first(row, ["Completeness", "completeness"]),
                "checkm2_contamination": first(row, ["Contamination", "contamination"]),
                "checkm2_model": first(row, ["Completeness_Model_Used", "Model", "model"]),
                "checkm2_coding_density": first(row, ["Coding_Density", "coding_density"]),
                "checkm2_genome_size": first(row, ["Genome_Size", "genome_size"]),
                "checkm2_gc_percent": first(row, ["GC_Content", "GC", "gc_content"]),
                "checkm2_notes": first(row, ["Additional_Notes", "Notes", "notes"]),
            }
            status, reason = checkm2_decision(stats)
            stats["checkm2_qc_status"] = status
            stats["checkm2_qc_fail_reasons"] = reason
            for key in key_candidates(name):
                rows[key] = stats
    return rows

metadata_keys = [
    "sequence_accession", "sequence_file", "assembly_accession", "Assembly Accession",
    "assembly", "accession", "genome_accession", "Genome accession", "Assembly", "sample_id",
]

checkm2_fields = [
    "checkm2_name", "checkm2_completeness", "checkm2_contamination",
    "checkm2_model", "checkm2_coding_density", "checkm2_genome_size",
    "checkm2_gc_percent", "checkm2_notes", "checkm2_qc_status",
    "checkm2_qc_fail_reasons", "combined_qc_status", "combined_qc_fail_reasons",
]

checkm2_by_key = read_checkm2(checkm2_path)

input_path = enriched_path if enriched_path.exists() else clean_path
if input_path.exists():
    with input_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        output_fields = fieldnames + [field for field in checkm2_fields if field not in fieldnames]
        rows = []
        for row in reader:
            match = None
            for key in metadata_keys:
                if key in row:
                    for candidate in key_candidates(row[key]):
                        if candidate in checkm2_by_key:
                            match = checkm2_by_key[candidate]
                            break
                if match:
                    break
            enriched = dict(row)
            for field in checkm2_fields:
                if field in {"combined_qc_status", "combined_qc_fail_reasons"}:
                    continue
                if match:
                    enriched[field] = match.get(field, "")
                elif field == "checkm2_qc_status":
                    enriched[field] = "FAIL"
                elif field == "checkm2_qc_fail_reasons":
                    enriched[field] = "NO_CHECKM2_RESULT"
                else:
                    enriched[field] = ""

            fail_reasons = []
            if enriched.get("sequence_qc_status") == "FAIL":
                fail_reasons.append(enriched.get("sequence_qc_fail_reasons") or "SEQUENCE_QC_FAIL")
            if enriched.get("checkm2_qc_status") == "FAIL":
                fail_reasons.append(enriched.get("checkm2_qc_fail_reasons") or "CHECKM2_QC_FAIL")
            enriched["combined_qc_status"] = "FAIL" if fail_reasons else "PASS"
            enriched["combined_qc_fail_reasons"] = ";".join(reason for reason in fail_reasons if reason)
            rows.append(enriched)
else:
    output_fields = checkm2_fields
    rows = []
    for stats in checkm2_by_key.values():
        row = dict(stats)
        row["combined_qc_status"] = stats["checkm2_qc_status"]
        row["combined_qc_fail_reasons"] = stats["checkm2_qc_fail_reasons"]
        rows.append(row)

pass_rows = [row for row in rows if row.get("combined_qc_status") == "PASS"]

with enriched_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=output_fields)
    writer.writeheader()
    writer.writerows(rows)

with pass_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=output_fields)
    writer.writeheader()
    writer.writerows(pass_rows)

with decision_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "sequence_file", "sequence_qc_status", "sequence_qc_fail_reasons",
            "checkm2_qc_status", "checkm2_qc_fail_reasons",
            "combined_qc_status", "combined_qc_fail_reasons",
        ],
        delimiter="\\t",
    )
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row.get(field, "") for field in writer.fieldnames})

if qc_filter:
    if clean_path.exists() and not unfiltered_path.exists():
        unfiltered_path.write_bytes(clean_path.read_bytes())
    if clean_path.exists():
        clean_path.write_bytes(pass_path.read_bytes())

    shutil.rmtree(filtered_sequence_dir, ignore_errors=True)
    filtered_sequence_dir.mkdir(parents=True, exist_ok=True)
    for row in pass_rows:
        sequence_file = row.get("sequence_file")
        if not sequence_file:
            continue
        source = source_sequence_dir / sequence_file
        target = filtered_sequence_dir / sequence_file
        if source.exists():
            target.write_bytes(source.read_bytes())
PY
    """
}

// Process 8: Run GTDB-Tk and add taxonomy match QC
process GTDBTK_QC {
    conda 'envs/gtdbtk.yaml'

    input:
    path sample_dir

    output:
    path "${sample_dir}", emit: gtdbtk_results

    script:
    def data_export = params.gtdbtk_data_path ? "export GTDBTK_DATA_PATH=${params.gtdbtk_data_path}" : "true"
    """
    mkdir -p ${sample_dir}/gtdbtk

    ${data_export}

    sequence_dir="${sample_dir}/sequence"
    if [ "${params.qc_filter}" = "true" ] && [ -d "${sample_dir}/sequence_filtered" ]; then
        sequence_dir="${sample_dir}/sequence_filtered"
    fi

    if [ -d "\${sequence_dir}" ] && [ -n "\$(find \${sequence_dir} -name "*.fna" -print -quit)" ]; then
        gtdbtk classify_wf \\
            --genome_dir \${sequence_dir} \\
            --out_dir ${sample_dir}/gtdbtk \\
            --extension fna \\
            --cpus ${params.threads}
    else
        echo "Warning: No .fna files found in \${sequence_dir}/ for GTDB-Tk" >&2
        printf "user_genome\\tclassification\\n" > ${sample_dir}/gtdbtk/gtdbtk.empty.summary.tsv
    fi

    python - <<'PY'
import csv
import shutil
from pathlib import Path

sample_dir = Path("${sample_dir}")
metadata_dir = sample_dir / "metadata_output"
gtdbtk_dir = sample_dir / "gtdbtk"
clean_path = metadata_dir / "ncbi_clean.csv"
enriched_path = metadata_dir / "ncbi_enriched.csv"
pass_path = metadata_dir / "ncbi_clean_qc_pass.csv"
unfiltered_path = metadata_dir / "ncbi_clean_unfiltered.csv"
decision_path = sample_dir / "sequence_qc" / "qc_decisions.tsv"
filtered_sequence_dir = sample_dir / "sequence_filtered"
source_sequence_dir = sample_dir / "sequence"

qc_filter = str("${params.qc_filter}").strip().lower() in {"true", "1", "yes", "y"}
match_rank = str("${params.taxonomy_match_rank}" or "genus").strip().lower()
if match_rank not in {"genus", "species"}:
    match_rank = "genus"

def normalize_key(value):
    return str(value or "").strip().lower().replace("-", "_").replace(".", "_")

def key_candidates(value):
    value = str(value or "").strip()
    if not value:
        return set()
    stem = Path(value).stem
    candidates = {value, stem}
    if stem.endswith("_genomic"):
        candidates.add(stem[:-8])
    parts = stem.split("_")
    if len(parts) >= 2 and parts[0] in {"GCF", "GCA"}:
        candidates.add("_".join(parts[:2]))
    return {normalize_key(candidate) for candidate in candidates if normalize_key(candidate)}

def normalize_taxon(value):
    value = str(value or "").strip()
    value = value.replace("[", "").replace("]", "")
    value = value.replace("_", " ")
    return " ".join(value.split()).lower()

def first(row, names):
    for name in names:
        if name in row:
            return row.get(name, "")
    return ""

def parse_gtdb_classification(classification):
    ranks = {}
    for item in str(classification or "").split(";"):
        item = item.strip()
        if "__" not in item:
            continue
        prefix, value = item.split("__", 1)
        ranks[prefix] = value.strip()
    return {
        "domain": ranks.get("d", ""),
        "phylum": ranks.get("p", ""),
        "class": ranks.get("c", ""),
        "order": ranks.get("o", ""),
        "family": ranks.get("f", ""),
        "genus": ranks.get("g", ""),
        "species": ranks.get("s", ""),
    }

def expected_from_metadata(row):
    expected = first(row, [
        "organism_name", "Organism Name", "organism", "Organism",
        "species", "Species", "scientific_name", "Scientific Name",
    ])
    tokens = normalize_taxon(expected).split()
    genus = tokens[0] if tokens else ""
    species = " ".join(tokens[:2]) if len(tokens) >= 2 else ""
    return expected, genus, species

def compare_taxonomy(row, match):
    expected_raw, expected_genus, expected_species = expected_from_metadata(row)
    observed_genus = normalize_taxon(match.get("gtdbtk_genus", ""))
    observed_species = normalize_taxon(match.get("gtdbtk_species", ""))

    if match_rank == "species":
        if not expected_species:
            return "FAIL", "NO_EXPECTED_SPECIES"
        if not observed_species:
            return "FAIL", "NO_GTDBTK_SPECIES"
        if observed_species == expected_species or observed_species.startswith(expected_species + " "):
            return "PASS", ""
        return "FAIL", f"SPECIES_MISMATCH:expected={expected_species};observed={observed_species}"

    if not expected_genus:
        return "FAIL", "NO_EXPECTED_GENUS"
    if not observed_genus:
        return "FAIL", "NO_GTDBTK_GENUS"
    if observed_genus == expected_genus:
        return "PASS", ""
    return "FAIL", f"GENUS_MISMATCH:expected={expected_genus};observed={observed_genus}"

def read_gtdbtk(path):
    rows = {}
    for summary_path in path.glob("*.summary.tsv"):
        with summary_path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\\t")
            for row in reader:
                genome = first(row, ["user_genome", "Name", "genome"])
                classification = first(row, ["classification", "Classification"])
                ranks = parse_gtdb_classification(classification)
                stats = {
                    "gtdbtk_user_genome": genome,
                    "gtdbtk_classification": classification,
                    "gtdbtk_domain": ranks["domain"],
                    "gtdbtk_phylum": ranks["phylum"],
                    "gtdbtk_class": ranks["class"],
                    "gtdbtk_order": ranks["order"],
                    "gtdbtk_family": ranks["family"],
                    "gtdbtk_genus": ranks["genus"],
                    "gtdbtk_species": ranks["species"],
                    "gtdbtk_classification_method": first(row, ["classification_method", "Classification_Method"]),
                    "gtdbtk_fastani_ani": first(row, ["fastani_ani", "FastANI_ANI"]),
                    "gtdbtk_fastani_reference": first(row, ["fastani_reference", "FastANI_Reference"]),
                    "gtdbtk_warnings": first(row, ["warnings", "Warnings"]),
                }
                for key in key_candidates(genome):
                    rows[key] = stats
    return rows

metadata_keys = [
    "sequence_accession", "sequence_file", "assembly_accession", "Assembly Accession",
    "assembly", "accession", "genome_accession", "Genome accession", "Assembly", "sample_id",
]

gtdbtk_fields = [
    "gtdbtk_user_genome", "gtdbtk_classification", "gtdbtk_domain",
    "gtdbtk_phylum", "gtdbtk_class", "gtdbtk_order", "gtdbtk_family",
    "gtdbtk_genus", "gtdbtk_species", "gtdbtk_classification_method",
    "gtdbtk_fastani_ani", "gtdbtk_fastani_reference", "gtdbtk_warnings",
    "gtdbtk_match_rank", "gtdbtk_qc_status", "gtdbtk_qc_fail_reasons",
    "combined_qc_status", "combined_qc_fail_reasons",
]

gtdbtk_by_key = read_gtdbtk(gtdbtk_dir)

input_path = enriched_path if enriched_path.exists() else clean_path
if input_path.exists():
    with input_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        output_fields = fieldnames + [field for field in gtdbtk_fields if field not in fieldnames]
        rows = []
        for row in reader:
            match = None
            for key in metadata_keys:
                if key in row:
                    for candidate in key_candidates(row[key]):
                        if candidate in gtdbtk_by_key:
                            match = gtdbtk_by_key[candidate]
                            break
                if match:
                    break
            enriched = dict(row)
            for field in gtdbtk_fields:
                if field in {"combined_qc_status", "combined_qc_fail_reasons"}:
                    continue
                if match:
                    enriched[field] = match.get(field, "")
                elif field == "gtdbtk_match_rank":
                    enriched[field] = match_rank
                elif field == "gtdbtk_qc_status":
                    enriched[field] = "FAIL"
                elif field == "gtdbtk_qc_fail_reasons":
                    enriched[field] = "NO_GTDBTK_RESULT"
                else:
                    enriched[field] = ""

            if match:
                status, reason = compare_taxonomy(enriched, match)
                enriched["gtdbtk_match_rank"] = match_rank
                enriched["gtdbtk_qc_status"] = status
                enriched["gtdbtk_qc_fail_reasons"] = reason

            fail_reasons = []
            if enriched.get("sequence_qc_status") == "FAIL":
                fail_reasons.append(enriched.get("sequence_qc_fail_reasons") or "SEQUENCE_QC_FAIL")
            if enriched.get("checkm2_qc_status") == "FAIL":
                fail_reasons.append(enriched.get("checkm2_qc_fail_reasons") or "CHECKM2_QC_FAIL")
            if enriched.get("gtdbtk_qc_status") == "FAIL":
                fail_reasons.append(enriched.get("gtdbtk_qc_fail_reasons") or "GTDBTK_QC_FAIL")
            enriched["combined_qc_status"] = "FAIL" if fail_reasons else "PASS"
            enriched["combined_qc_fail_reasons"] = ";".join(reason for reason in fail_reasons if reason)
            rows.append(enriched)
else:
    output_fields = gtdbtk_fields
    rows = []
    for stats in gtdbtk_by_key.values():
        row = dict(stats)
        row["gtdbtk_match_rank"] = match_rank
        row["gtdbtk_qc_status"] = "FAIL"
        row["gtdbtk_qc_fail_reasons"] = "NO_METADATA"
        row["combined_qc_status"] = "FAIL"
        row["combined_qc_fail_reasons"] = "NO_METADATA"
        rows.append(row)

pass_rows = [row for row in rows if row.get("combined_qc_status") == "PASS"]

with enriched_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=output_fields)
    writer.writeheader()
    writer.writerows(rows)

with pass_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=output_fields)
    writer.writeheader()
    writer.writerows(pass_rows)

with decision_path.open("w", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "sequence_file", "sequence_qc_status", "sequence_qc_fail_reasons",
            "checkm2_qc_status", "checkm2_qc_fail_reasons",
            "gtdbtk_match_rank", "gtdbtk_qc_status", "gtdbtk_qc_fail_reasons",
            "combined_qc_status", "combined_qc_fail_reasons",
        ],
        delimiter="\\t",
    )
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row.get(field, "") for field in writer.fieldnames})

if qc_filter:
    if clean_path.exists() and not unfiltered_path.exists():
        unfiltered_path.write_bytes(clean_path.read_bytes())
    if clean_path.exists():
        clean_path.write_bytes(pass_path.read_bytes())

    shutil.rmtree(filtered_sequence_dir, ignore_errors=True)
    filtered_sequence_dir.mkdir(parents=True, exist_ok=True)
    for row in pass_rows:
        sequence_file = row.get("sequence_file")
        if not sequence_file:
            continue
        source = source_sequence_dir / sequence_file
        target = filtered_sequence_dir / sequence_file
        if source.exists():
            target.write_bytes(source.read_bytes())
PY
    """
}

// Process 9: Optional QUAST assembly-structure QC
process QUAST_QC {
    conda 'envs/quast.yaml'

    input:
    path sample_dir

    output:
    path "${sample_dir}", emit: quast_results

    script:
    """
    mkdir -p ${sample_dir}/quast
    sequence_dir="${sample_dir}/sequence"
    if [ "${params.qc_filter}" = "true" ] && [ -d "${sample_dir}/sequence_filtered" ]; then
        sequence_dir="${sample_dir}/sequence_filtered"
    fi

    if [ -d "\${sequence_dir}" ] && [ -n "\$(find \${sequence_dir} -name "*.fna" -print -quit)" ]; then
        quast.py -t ${params.threads} -o ${sample_dir}/quast \${sequence_dir}/*.fna || true
    else
        echo "Warning: No .fna files found in \${sequence_dir}/ for QUAST" >&2
    fi
    python ${baseDir}/scripts/quast_summary.py --sample-dir ${sample_dir}
    """
}

// Process 10: Optional FastANI/skani comparative ANI
process ANI_ANALYSIS {
    conda 'envs/ani.yaml'

    input:
    path sample_dir

    output:
    path "${sample_dir}", emit: ani_results

    script:
    def aniStrategy = params.ani_large_run_strategy ?: "auto"
    def aniMaxAllVsAll = params.ani_max_all_vs_all_genomes ? params.ani_max_all_vs_all_genomes as int : 200
    """
    mkdir -p ${sample_dir}/ani/analysis
    sequence_dir="${sample_dir}/sequence"
    if [ "${params.qc_filter}" = "true" ] && [ -d "${sample_dir}/sequence_filtered" ]; then
        sequence_dir="${sample_dir}/sequence_filtered"
    fi
    find \${sequence_dir} -name "*.fna" | sort > ${sample_dir}/ani/genomes.list || true
    genome_count=\$(wc -l < ${sample_dir}/ani/genomes.list || echo 0)
    estimated_comparisons=\$((genome_count * genome_count))
    ani_decision="run_all_vs_all"
    ani_status="PASS"
    ani_message="Running all-vs-all ANI."
    if [ "\${genome_count}" -ge 2 ]; then
        if [ "${aniStrategy}" = "skip" ]; then
            ani_decision="skip_requested"
            ani_status="SKIPPED"
            ani_message="ANI skipped because --ani_large_run_strategy skip was requested."
        elif [ "${aniStrategy}" = "auto" ] && [ "${params.large_dataset}" = "true" ] && [ "\${genome_count}" -gt "${aniMaxAllVsAll}" ]; then
            ani_decision="skip_large_dataset"
            ani_status="SKIPPED_LARGE_DATASET"
            ani_message="ANI skipped automatically: \${genome_count} genomes exceed --ani_max_all_vs_all_genomes ${aniMaxAllVsAll} in large-dataset mode."
        fi

        printf "tool\\tgenome_count\\testimated_comparisons\\tstrategy\\tmax_all_vs_all_genomes\\tlarge_dataset\\tdecision\\tstatus\\tmessage\\n" > ${sample_dir}/ani/analysis/ani_run_status.tsv
        printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" "${params.ani_tool}" "\${genome_count}" "\${estimated_comparisons}" "${aniStrategy}" "${aniMaxAllVsAll}" "${params.large_dataset}" "\${ani_decision}" "\${ani_status}" "\${ani_message}" >> ${sample_dir}/ani/analysis/ani_run_status.tsv

        if [ "\${ani_status}" = "SKIPPED" ] || [ "\${ani_status}" = "SKIPPED_LARGE_DATASET" ]; then
            printf "query\\treference\\tani\\tfragments_mapped\\tfragments_total\\n" > ${sample_dir}/ani/fastani_pairs.tsv
            python ${baseDir}/scripts/ani_summary.py \\
                --sample-dir ${sample_dir} \\
                --pairs ${sample_dir}/ani/fastani_pairs.tsv \\
                --genomes-list ${sample_dir}/ani/genomes.list \\
                --tool ${params.ani_tool} \\
                --duplicate-threshold ${params.ani_duplicate_threshold} \\
                --species-threshold ${params.ani_species_threshold}
        else
            if [ "${params.ani_tool}" = "skani" ]; then
                skani triangle -t ${params.threads} \$(cat ${sample_dir}/ani/genomes.list) > ${sample_dir}/ani/skani_pairs.tsv || true
                pair_file="${sample_dir}/ani/skani_pairs.tsv"
            else
                fastANI --ql ${sample_dir}/ani/genomes.list --rl ${sample_dir}/ani/genomes.list -t ${params.threads} -o ${sample_dir}/ani/fastani_pairs.tsv || true
                pair_file="${sample_dir}/ani/fastani_pairs.tsv"
            fi
            python ${baseDir}/scripts/ani_summary.py \\
                --sample-dir ${sample_dir} \\
                --pairs \${pair_file} \\
                --genomes-list ${sample_dir}/ani/genomes.list \\
                --tool ${params.ani_tool} \\
                --duplicate-threshold ${params.ani_duplicate_threshold} \\
                --species-threshold ${params.ani_species_threshold}
        fi
    else
        printf "tool\\tgenome_count\\testimated_comparisons\\tstrategy\\tmax_all_vs_all_genomes\\tlarge_dataset\\tdecision\\tstatus\\tmessage\\n" > ${sample_dir}/ani/analysis/ani_run_status.tsv
        printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" "${params.ani_tool}" "\${genome_count}" "\${estimated_comparisons}" "${aniStrategy}" "${aniMaxAllVsAll}" "${params.large_dataset}" "insufficient_genomes" "SKIPPED_INAPPLICABLE" "ANI requires at least two genomes; pairwise ANI was skipped." >> ${sample_dir}/ani/analysis/ani_run_status.tsv
        printf "query\\treference\\tani\\tfragments_mapped\\tfragments_total\\n" > ${sample_dir}/ani/fastani_pairs.tsv
        python ${baseDir}/scripts/ani_summary.py --sample-dir ${sample_dir} --pairs ${sample_dir}/ani/fastani_pairs.tsv --genomes-list ${sample_dir}/ani/genomes.list --tool ${params.ani_tool}
    fi
    """
}

// Process 11: Optional Mash sketching pre-screen
process MASH_PRESCREEN {
    conda 'envs/mash.yaml'

    input:
    path sample_dir

    output:
    path "${sample_dir}", emit: mash_results

    script:
    """
    mkdir -p ${sample_dir}/mash/analysis
    sequence_dir="${sample_dir}/sequence"
    if [ "${params.qc_filter}" = "true" ] && [ -d "${sample_dir}/sequence_filtered" ]; then
        sequence_dir="${sample_dir}/sequence_filtered"
    fi
    find \${sequence_dir} -name "*.fna" | sort > ${sample_dir}/mash/genomes.list || true
    genome_count=\$(wc -l < ${sample_dir}/mash/genomes.list || echo 0)
    if [ "\${genome_count}" -ge 2 ]; then
        mash sketch -o ${sample_dir}/mash/genomes \$(cat ${sample_dir}/mash/genomes.list)
        mash dist ${sample_dir}/mash/genomes.msh \$(cat ${sample_dir}/mash/genomes.list) > ${sample_dir}/mash/mash_dist.tsv || true
    else
        touch ${sample_dir}/mash/mash_dist.tsv
    fi
    python ${baseDir}/scripts/mash_summary.py \\
        --sample-dir ${sample_dir} \\
        --dist ${sample_dir}/mash/mash_dist.tsv \\
        --genomes-list ${sample_dir}/mash/genomes.list
    """
}

process MOBSUITE_ANALYSIS {
    conda 'envs/mobsuite.yaml'

    input:
    path sample_dir

    output:
    path "${sample_dir}", emit: mobsuite_results

    script:
    def mobsuiteDbPath = effectiveMobsuiteDbPath()
    def mobsuiteAutoInit = params.mobsuite_db ? "false" : params.mobsuite_auto_init_db
    def jobs = params.mobsuite_jobs ? params.mobsuite_jobs as int : Math.min(params.threads as int, 8)
    def threadsPerSample = params.mobsuite_threads_per_sample ? params.mobsuite_threads_per_sample as int : 1
    """
    sequence_dir="${sample_dir}/sequence"
    if [ "${params.qc_filter}" = "true" ] && [ -d "${sample_dir}/sequence_filtered" ]; then
        sequence_dir="${sample_dir}/sequence_filtered"
    fi

    mkdir -p ${sample_dir}/mobsuite/raw ${sample_dir}/mobsuite/tables
    export HOME="\${PWD}/mobsuite_home"
    mkdir -p "\${HOME}/.etetoolkit"
    status_file=${sample_dir}/mobsuite/module_status.tsv
    printf "module\\tenabled\\tstarted\\tcompleted\\tstatus\\tsamples_input\\tsamples_processed\\tsamples_failed\\traw_tables_created\\tfeature_rows_created\\tunique_features_created\\toutput_dir\\tmessage\\n" > "\${status_file}"
    started=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")
    python ${baseDir}/scripts/setup_mobsuite_database.py \\
        --db-dir ${shellQuote(mobsuiteDbPath)} \\
        --out ${sample_dir}/mobsuite/mobsuite_database_setup_status.tsv \\
        --auto-init ${mobsuiteAutoInit} \\
        --auto-init-taxa ${params.mobsuite_auto_init_taxa}
    db_setup_status=\$(awk -F '\\t' 'NR == 2 {print \$8}' ${sample_dir}/mobsuite/mobsuite_database_setup_status.tsv 2>/dev/null || echo FAIL)
    samples_input=0
    samples_processed=0
    samples_failed=0
    if [ "\${db_setup_status}" = "FAIL" ]; then
        echo "MOB-suite database setup failed; inspect ${sample_dir}/mobsuite/mobsuite_database_setup_status.tsv." > ${sample_dir}/mobsuite/tables/mobsuite_warning.txt
    elif [ -d "\${sequence_dir}" ]; then
        echo "Running MOB-suite with ${jobs} parallel sample job(s) and ${threadsPerSample} thread(s) per sample."
        python ${baseDir}/scripts/run_mobsuite_parallel.py \\
            --sequence-dir "\${sequence_dir}" \\
            --raw-dir ${sample_dir}/mobsuite/raw \\
            --status-file ${sample_dir}/mobsuite/tables/mobsuite_sample_status.tsv \\
            --database-directory ${shellQuote(mobsuiteDbPath)} \\
            --jobs ${jobs} \\
            --threads-per-sample ${threadsPerSample} \\
            --reuse-existing ${params.mobsuite_reuse_existing} || true
        samples_input=\$(awk 'NR > 1 {count++} END {print count + 0}' ${sample_dir}/mobsuite/tables/mobsuite_sample_status.tsv 2>/dev/null || echo 0)
        samples_failed=\$(awk -F '\\t' 'NR > 1 && \$4 == "FAIL" {count++} END {print count + 0}' ${sample_dir}/mobsuite/tables/mobsuite_sample_status.tsv 2>/dev/null || echo 0)
        samples_processed=\$(awk -F '\\t' 'NR > 1 && (\$4 == "PASS" || \$4 == "REUSED") {count++} END {print count + 0}' ${sample_dir}/mobsuite/tables/mobsuite_sample_status.tsv 2>/dev/null || echo 0)
    else
        echo "No sequence directory was found for MOB-suite." > ${sample_dir}/mobsuite/tables/mobsuite_warning.txt
    fi
    python ${baseDir}/scripts/collect_optional_tool_tables.py \\
        --raw-dir ${sample_dir}/mobsuite/raw \\
        --out ${sample_dir}/mobsuite/tables/mobsuite.tsv \\
        --tool mobsuite
    feature_rows=\$(awk 'NR > 1 {count++} END {print count + 0}' ${sample_dir}/mobsuite/tables/mobsuite.tsv 2>/dev/null || echo 0)
    unique_features=\$(awk -F '\\t' 'NR > 1 {for (i=1; i<=NF; i++) if (\$i != "") seen[\$i]=1} END {print length(seen)}' ${sample_dir}/mobsuite/tables/mobsuite.tsv 2>/dev/null || echo 0)
    status=PASS
    message="MOB-suite completed and parseable tables were collected."
    if [ "\${db_setup_status}" = "FAIL" ]; then
        status=FAIL
        message="MOB-suite database setup failed; inspect mobsuite/mobsuite_database_setup_status.tsv."
    fi
    if [ "\${status}" = "FAIL" ]; then
        :
    elif [ "\${samples_input}" -eq 0 ]; then
        status=WARNING_EMPTY
        message="No FASTA files found for MOB-suite."
    elif [ "\${samples_failed}" -gt 0 ]; then
        status=WARNING_PARTIAL
        message="MOB-suite had sample failures; inspect mobsuite/raw/*/mob_recon.stderr."
    fi
    if [ "\${feature_rows}" -eq 0 ]; then
        status=WARNING_EMPTY
        message="\${message} No MOB-suite feature rows were collected."
    fi
    printf "mobsuite\\ttrue\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
        "\${started}" "\$(date -u +"%Y-%m-%dT%H:%M:%SZ")" "\${status}" "\${samples_input}" "\${samples_processed}" "\${samples_failed}" "\${samples_processed}" "\${feature_rows}" "\${unique_features}" "${sample_dir}/mobsuite" "\${message}" >> "\${status_file}"
    """
}

process GENOMAD_PROPHAGE {
    conda 'envs/genomad.yaml'

    input:
    path sample_dir

    output:
    path "${sample_dir}", emit: genomad_results

    script:
    def genomadDbPath = effectiveGenomadDbPath()
    def genomadAutoDownload = params.genomad_db ? "false" : params.genomad_auto_download_db
    def genomadExtraArgs = []
    if (params.genomad_splits != null && params.genomad_splits.toString().trim()) {
        genomadExtraArgs << "--splits ${params.genomad_splits}"
    }
    if (params.genomad_sensitivity != null && params.genomad_sensitivity.toString().trim()) {
        genomadExtraArgs << "--sensitivity ${params.genomad_sensitivity}"
    }
    def genomadExtraArgString = genomadExtraArgs.join(' ')
    def genomadJobs = params.genomad_jobs ? params.genomad_jobs as int : Math.min(params.threads as int, 2)
    def genomadThreadsPerSample = params.genomad_threads_per_sample ? params.genomad_threads_per_sample as int : 1
    def genomadReuseExisting = params.genomad_reuse_existing
    """
    sequence_dir="${sample_dir}/sequence"
    if [ "${params.qc_filter}" = "true" ] && [ -d "${sample_dir}/sequence_filtered" ]; then
        sequence_dir="${sample_dir}/sequence_filtered"
    fi

    mkdir -p ${sample_dir}/prophage/raw ${sample_dir}/prophage/tables
    status_file=${sample_dir}/prophage/module_status.tsv
    sample_status_file=${sample_dir}/prophage/genomad_sample_status.tsv
    printf "module\\tenabled\\tstarted\\tcompleted\\tstatus\\tsamples_input\\tsamples_processed\\tsamples_failed\\traw_tables_created\\tfeature_rows_created\\tunique_features_created\\toutput_dir\\tmessage\\n" > "\${status_file}"
    printf "sample_id\\tstatus\\tseconds\\toutput_dir\\tmessage\\n" > "\${sample_status_file}"
    started=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")
    python ${baseDir}/scripts/setup_genomad_database.py \\
        --db-dir ${shellQuote(genomadDbPath)} \\
        --out ${sample_dir}/prophage/genomad_database_setup_status.tsv \\
        --auto-download ${genomadAutoDownload}
    genomad_db_path=\$(awk -F '\\t' 'NR == 2 {print \$2}' ${sample_dir}/prophage/genomad_database_setup_status.tsv 2>/dev/null || true)
    db_setup_status=\$(awk -F '\\t' 'NR == 2 {print \$6}' ${sample_dir}/prophage/genomad_database_setup_status.tsv 2>/dev/null || echo FAIL)
    samples_input=0
    samples_processed=0
    samples_failed=0
    if [ "\${db_setup_status}" = "FAIL" ]; then
        echo "geNomad database setup failed; inspect ${sample_dir}/prophage/genomad_database_setup_status.tsv." > ${sample_dir}/prophage/tables/genomad_warning.txt
    elif command -v genomad >/dev/null 2>&1 && [ -d "\${sequence_dir}" ]; then
        genomad_jobs=${genomadJobs}
        if [ "\${genomad_jobs}" -lt 1 ]; then
            genomad_jobs=1
        fi
        run_genomad_one() {
            fasta="\$1"
            prefix=\$(basename "\${fasta}" .fna)
            outdir="${sample_dir}/prophage/raw/\${prefix}"
            start_epoch=\$(date +%s)
            mkdir -p "\${outdir}"
            if [ "${genomadReuseExisting}" = "true" ] && find "\${outdir}" -name "*_summary.tsv" -size +0 -print -quit | grep -q .; then
                seconds=\$((\$(date +%s) - start_epoch))
                printf "%s\\tPASS_REUSED\\t%s\\t%s\\t%s\\n" "\${prefix}" "\${seconds}" "\${outdir}" "existing geNomad summary reused" >> "\${sample_status_file}"
            elif genomad end-to-end "\${fasta}" "\${outdir}" "\${genomad_db_path}" --threads ${genomadThreadsPerSample} ${genomadExtraArgString}; then
                seconds=\$((\$(date +%s) - start_epoch))
                printf "%s\\tPASS\\t%s\\t%s\\t%s\\n" "\${prefix}" "\${seconds}" "\${outdir}" "geNomad completed" >> "\${sample_status_file}"
            else
                seconds=\$((\$(date +%s) - start_epoch))
                printf "%s\\tFAIL\\t%s\\t%s\\t%s\\n" "\${prefix}" "\${seconds}" "\${outdir}" "geNomad failed; inspect raw output/logs" >> "\${sample_status_file}"
            fi
        }
        active_jobs=0
        for fasta in \$(find "\${sequence_dir}" -name "*.fna" | sort); do
            samples_input=\$((samples_input + 1))
            run_genomad_one "\${fasta}" &
            active_jobs=\$((active_jobs + 1))
            if [ "\${active_jobs}" -ge "\${genomad_jobs}" ]; then
                wait -n || true
                active_jobs=\$((active_jobs - 1))
            fi
        done
        while [ "\${active_jobs}" -gt 0 ]; do
            wait -n || true
            active_jobs=\$((active_jobs - 1))
        done
        samples_processed=\$(awk -F '\\t' 'NR > 1 && (\$2 == "PASS" || \$2 == "PASS_REUSED") {count++} END {print count + 0}' "\${sample_status_file}" 2>/dev/null || echo 0)
        samples_failed=\$(awk -F '\\t' 'NR > 1 && \$2 == "FAIL" {count++} END {print count + 0}' "\${sample_status_file}" 2>/dev/null || echo 0)
    else
        echo "geNomad executable was not available or no sequence directory was found." > ${sample_dir}/prophage/tables/genomad_warning.txt
    fi
    python ${baseDir}/scripts/collect_optional_tool_tables.py \\
        --raw-dir ${sample_dir}/prophage/raw \\
        --out ${sample_dir}/prophage/tables/prophage.tsv \\
        --tool prophage
    feature_rows=\$(awk 'NR > 1 {count++} END {print count + 0}' ${sample_dir}/prophage/tables/prophage.tsv 2>/dev/null || echo 0)
    unique_features=\$(awk -F '\\t' 'NR > 1 {for (i=1; i<=NF; i++) if (\$i != "") seen[\$i]=1} END {print length(seen)}' ${sample_dir}/prophage/tables/prophage.tsv 2>/dev/null || echo 0)
    status=PASS
    message="geNomad completed and parseable tables were collected."
    if [ "\${db_setup_status}" = "FAIL" ]; then
        status=FAIL
        message="geNomad database setup failed; inspect prophage/genomad_database_setup_status.tsv."
    elif [ "\${samples_input}" -eq 0 ]; then
        status=WARNING_EMPTY
        message="No FASTA files found for geNomad."
    elif [ "\${samples_failed}" -gt 0 ]; then
        status=WARNING_PARTIAL
        message="geNomad had sample failures; inspect prophage/raw/* outputs."
    fi
    if [ "\${status}" != "FAIL" ] && [ "\${feature_rows}" -eq 0 ]; then
        status=WARNING_EMPTY
        message="\${message} No geNomad feature rows were collected."
    fi
    printf "genomad\\ttrue\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
        "\${started}" "\$(date -u +"%Y-%m-%dT%H:%M:%SZ")" "\${status}" "\${samples_input}" "\${samples_processed}" "\${samples_failed}" "\${samples_processed}" "\${feature_rows}" "\${unique_features}" "${sample_dir}/prophage" "\${message}" >> "\${status_file}"
    """
}

process ISFINDER_BLAST {
    conda 'envs/panr2_comprehensive.yaml'

    input:
    path sample_dir

    output:
    path "${sample_dir}", emit: isfinder_results

    script:
    def isfinderDbFasta = params.isfinder_db_fasta ? launchPath(params.isfinder_db_fasta) : ""
    def dbVersion = params.isfinder_database_version ?: ""
    """
    sequence_dir="${sample_dir}/sequence"
    if [ "${params.qc_filter}" = "true" ] && [ -d "${sample_dir}/sequence_filtered" ]; then
        sequence_dir="${sample_dir}/sequence_filtered"
    fi

    mkdir -p ${sample_dir}/isfinder/raw ${sample_dir}/isfinder/tables ${sample_dir}/isfinder/db
    status_file=${sample_dir}/isfinder/module_status.tsv
    printf "module\\tenabled\\tstarted\\tcompleted\\tstatus\\tsamples_input\\tsamples_processed\\tsamples_failed\\traw_tables_created\\tfeature_rows_created\\tunique_features_created\\toutput_dir\\tmessage\\n" > "\${status_file}"
    started=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")

    if [ -z "${isfinderDbFasta}" ] || [ ! -s "${isfinderDbFasta}" ]; then
        printf "isfinder\\ttrue\\t%s\\t%s\\tFAIL\\t0\\t0\\t0\\t0\\t0\\t0\\t%s\\t%s\\n" "\${started}" "\$(date -u +"%Y-%m-%dT%H:%M:%SZ")" "${sample_dir}/isfinder" "Missing --isfinder_db_fasta. ISfinder database download/redistribution requires written authorization; provide an authorized local FASTA." >> "\${status_file}"
        echo "Missing --isfinder_db_fasta. ISfinder database download/redistribution requires written authorization; provide an authorized local FASTA." >&2
        exit 1
    fi

    makeblastdb -in "${isfinderDbFasta}" -dbtype nucl -out ${sample_dir}/isfinder/db/isfinder >/dev/null

    samples_input=0
    samples_processed=0
    samples_failed=0
    raw_tables_created=0
    for fasta in \$(find "\${sequence_dir}" -name "*.fna" | sort 2>/dev/null || true); do
        samples_input=\$((samples_input + 1))
        prefix=\$(basename "\${fasta}" .fna)
        raw=${sample_dir}/isfinder/raw/\${prefix}.blast.tsv
        results=${sample_dir}/isfinder/tables/\${prefix}_results.tab
        summary=${sample_dir}/isfinder/tables/\${prefix}_summary.tab
        if blastn \\
            -query "\${fasta}" \\
            -db ${sample_dir}/isfinder/db/isfinder \\
            -out "\${raw}" \\
            -outfmt "6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore" \\
            -num_threads ${params.threads}; then
            samples_processed=\$((samples_processed + 1))
            raw_tables_created=\$((raw_tables_created + 1))
        else
            samples_failed=\$((samples_failed + 1))
            : > "\${raw}"
        fi
        python ${baseDir}/scripts/isfinder_blast_to_abricate.py \\
            --blast "\${raw}" \\
            --sample-id "\${prefix}" \\
            --out-results "\${results}" \\
            --out-summary "\${summary}" \\
            --min-identity ${params.isfinder_min_identity} \\
            --min-coverage ${params.isfinder_min_coverage} \\
            --database-version ${shellQuote(dbVersion)}
    done

    feature_rows=\$(awk 'FNR > 1 {count++} END {print count + 0}' ${sample_dir}/isfinder/tables/*_results.tab 2>/dev/null || echo 0)
    unique_features=\$(awk -F '\\t' 'FNR > 1 && \$5 != "" {seen[\$5]=1} END {print length(seen)}' ${sample_dir}/isfinder/tables/*_results.tab 2>/dev/null || echo 0)
    status=PASS
    message="ISfinder-compatible BLAST completed with authorized local database FASTA."
    if [ "\${samples_input}" -eq 0 ]; then
        status=WARNING_EMPTY
        message="No FASTA files found for ISfinder-compatible BLAST."
    fi
    printf "isfinder\\ttrue\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
        "\${started}" "\$(date -u +"%Y-%m-%dT%H:%M:%SZ")" "\${status}" "\${samples_input}" "\${samples_processed}" "\${samples_failed}" "\${raw_tables_created}" "\${feature_rows}" "\${unique_features}" "${sample_dir}/isfinder" "\${message}" >> "\${status_file}"
    """
}

process ORGANISM_SPECIFIC_TYPING {
    conda 'envs/organism_typing.yaml'

    input:
    path sample_dir

    output:
    path "${sample_dir}", emit: organism_typing_results

    script:
    def runKleborate = (params.run_organism_specific_typing || params.run_kleborate) ? "true" : "false"
    def runKaptive = (params.run_organism_specific_typing || params.run_kaptive) ? "true" : "false"
    def runEctyper = (params.run_organism_specific_typing || params.run_ectyper) ? "true" : "false"
    def kaptiveDbCheck = params.kaptive_db ? "true" : "false"
    def kleborateJobs = params.kleborate_jobs ? params.kleborate_jobs as int : Math.min(params.threads as int, 8)
    """
    sequence_dir="${sample_dir}/sequence"
    if [ "${params.qc_filter}" = "true" ] && [ -d "${sample_dir}/sequence_filtered" ]; then
        sequence_dir="${sample_dir}/sequence_filtered"
    fi
    fasta_files=\$(find "\${sequence_dir}" -name "*.fna" | sort 2>/dev/null || true)

    if [ "${runKleborate}" = "true" ]; then
        mkdir -p ${sample_dir}/kleborate/raw ${sample_dir}/kleborate/tables
        kleborate_status_file=${sample_dir}/kleborate/module_status.tsv
        printf "module\\tenabled\\tstarted\\tcompleted\\tstatus\\tsamples_input\\tsamples_processed\\tsamples_failed\\traw_tables_created\\tfeature_rows_created\\tunique_features_created\\toutput_dir\\tmessage\\n" > "\${kleborate_status_file}"
        kleborate_started=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")
        if command -v kleborate >/dev/null 2>&1 && [ -n "\${fasta_files}" ]; then
            echo "Running Kleborate with ${kleborateJobs} parallel sample job(s)."
            python ${baseDir}/scripts/run_kleborate_parallel.py \\
                --sequence-dir "\${sequence_dir}" \\
                --raw-dir ${sample_dir}/kleborate/raw \\
                --status-file ${sample_dir}/kleborate/tables/kleborate_sample_status.tsv \\
                --jobs ${kleborateJobs} \\
                --reuse-existing ${params.kleborate_reuse_existing} || true
            python ${baseDir}/scripts/collect_optional_tool_tables.py --raw-dir ${sample_dir}/kleborate/raw --out ${sample_dir}/kleborate/tables/kleborate.tsv --tool kleborate
            kleborate_samples_input=\$(awk 'NR > 1 {count++} END {print count + 0}' ${sample_dir}/kleborate/tables/kleborate_sample_status.tsv 2>/dev/null || echo 0)
            kleborate_samples_failed=\$(awk -F '\\t' 'NR > 1 && \$4 == "FAIL" {count++} END {print count + 0}' ${sample_dir}/kleborate/tables/kleborate_sample_status.tsv 2>/dev/null || echo 0)
            kleborate_samples_processed=\$(awk -F '\\t' 'NR > 1 && (\$4 == "PASS" || \$4 == "REUSED") {count++} END {print count + 0}' ${sample_dir}/kleborate/tables/kleborate_sample_status.tsv 2>/dev/null || echo 0)
        else
            printf "sample_id\\tstatus\\n" > ${sample_dir}/kleborate/tables/kleborate.tsv
            kleborate_samples_input=0
            kleborate_samples_processed=0
            kleborate_samples_failed=0
        fi
        kleborate_feature_rows=\$(awk 'NR > 1 {count++} END {print count + 0}' ${sample_dir}/kleborate/tables/kleborate.tsv 2>/dev/null || echo 0)
        kleborate_unique_features=\$(awk -F '\\t' 'NR > 1 {for (i=1; i<=NF; i++) if (\$i != "") seen[\$i]=1} END {print length(seen)}' ${sample_dir}/kleborate/tables/kleborate.tsv 2>/dev/null || echo 0)
        kleborate_status=PASS
        kleborate_message="Kleborate completed and parseable tables were collected."
        if [ "\${kleborate_samples_input}" -eq 0 ]; then
            kleborate_status=WARNING_EMPTY
            kleborate_message="Kleborate executable was unavailable or no FASTA files were found."
        elif [ "\${kleborate_samples_failed}" -gt 0 ]; then
            kleborate_status=WARNING_PARTIAL
            kleborate_message="Kleborate had sample failures; inspect kleborate/raw/*/kleborate.stderr."
        fi
        if [ "\${kleborate_feature_rows}" -eq 0 ]; then
            kleborate_status=WARNING_EMPTY
            kleborate_message="\${kleborate_message} No Kleborate feature rows were collected."
        fi
        printf "kleborate\\ttrue\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \\
            "\${kleborate_started}" "\$(date -u +"%Y-%m-%dT%H:%M:%SZ")" "\${kleborate_status}" "\${kleborate_samples_input}" "\${kleborate_samples_processed}" "\${kleborate_samples_failed}" "\${kleborate_samples_processed}" "\${kleborate_feature_rows}" "\${kleborate_unique_features}" "${sample_dir}/kleborate" "\${kleborate_message}" >> "\${kleborate_status_file}"
    fi

    if [ "${runKaptive}" = "true" ]; then
        mkdir -p ${sample_dir}/kaptive/raw ${sample_dir}/kaptive/tables
        if [ "${kaptiveDbCheck}" = "true" ] && command -v kaptive >/dev/null 2>&1 && [ -n "\${fasta_files}" ]; then
            for fasta in \${fasta_files}; do
                prefix=\$(basename "\${fasta}" .fna)
                kaptive assembly "${params.kaptive_db}" "\${fasta}" > ${sample_dir}/kaptive/raw/\${prefix}.tsv || true
            done
            python ${baseDir}/scripts/collect_optional_tool_tables.py --raw-dir ${sample_dir}/kaptive/raw --out ${sample_dir}/kaptive/tables/kaptive.tsv --tool kaptive
        else
            printf "sample_id\\tstatus\\n" > ${sample_dir}/kaptive/tables/kaptive.tsv
        fi
    fi

    if [ "${runEctyper}" = "true" ]; then
        mkdir -p ${sample_dir}/ectyper/raw ${sample_dir}/ectyper/tables
        if command -v ectyper >/dev/null 2>&1 && [ -n "\${fasta_files}" ]; then
            ectyper -i \${fasta_files} -o ${sample_dir}/ectyper/raw --cores ${params.threads} || true
            python ${baseDir}/scripts/collect_optional_tool_tables.py --raw-dir ${sample_dir}/ectyper/raw --out ${sample_dir}/ectyper/tables/ectyper.tsv --tool ectyper
        else
            printf "sample_id\\tstatus\\n" > ${sample_dir}/ectyper/tables/ectyper.tsv
        fi
    fi
    """
}

process AMRFINDERPLUS_ANALYSIS {
    conda 'envs/amrfinderplus.yaml'

    input:
    path sample_dir

    output:
    path "${sample_dir}", emit: amrfinderplus_results

    script:
    def organismArg = params.amrfinderplus_organism ? "--organism ${shellQuote(params.amrfinderplus_organism)}" : ""
    def jobs = params.amrfinderplus_jobs ? params.amrfinderplus_jobs as int : params.threads as int
    def threadsPerSample = params.amrfinderplus_threads_per_sample ? params.amrfinderplus_threads_per_sample as int : 1
    """
    sequence_dir="${sample_dir}/sequence"
    if [ "${params.qc_filter}" = "true" ] && [ -d "${sample_dir}/sequence_filtered" ]; then
        sequence_dir="${sample_dir}/sequence_filtered"
    fi

    mkdir -p ${sample_dir}/amrfinderplus/raw ${sample_dir}/amrfinderplus/tables
    status_file="${sample_dir}/amrfinderplus/tables/amrfinderplus_sample_status.tsv"
    echo "Running AMRFinderPlus with ${jobs} parallel sample job(s) and ${threadsPerSample} thread(s) per sample."
    python ${baseDir}/scripts/run_amrfinderplus_parallel.py \\
        --sequence-dir "\${sequence_dir}" \\
        --raw-dir ${sample_dir}/amrfinderplus/raw \\
        --status-file "\${status_file}" \\
        --jobs ${jobs} \\
        --threads-per-sample ${threadsPerSample} \\
        --update-db ${params.amrfinderplus_update_db} \\
        --reuse-existing ${params.amrfinderplus_reuse_existing} \\
        --progress-every ${params.amrfinderplus_progress_every} \\
        ${organismArg}
    python ${baseDir}/scripts/collect_optional_tool_tables.py \\
        --raw-dir ${sample_dir}/amrfinderplus/raw \\
        --out ${sample_dir}/amrfinderplus/tables/amrfinderplus.tsv \\
        --tool amrfinderplus
    """
}

// Process 12: Combined QC decision engine and optional filtering
process COMBINED_QC {
    conda 'envs/fetchm.yaml'

    input:
    path sample_dir

    output:
    path "${sample_dir}", emit: combined_qc_results

    script:
    def max_contigs_arg = params.max_contigs ? "--max-contigs ${params.max_contigs}" : ""
    def min_n50_arg = params.min_n50 ? "--min-n50 ${params.min_n50}" : ""
    def min_completeness_arg = params.min_completeness ? "--min-completeness ${params.min_completeness}" : ""
    def max_contamination_arg = params.max_contamination ? "--max-contamination ${params.max_contamination}" : ""
    """
    python ${baseDir}/scripts/qc_master.py \\
        --sample-dir ${sample_dir} \\
        --qc-filter ${params.qc_filter} \\
        --representative-only ${params.representative_only} \\
        --ani-species-threshold ${params.ani_species_threshold} \\
        ${max_contigs_arg} ${min_n50_arg} ${min_completeness_arg} ${max_contamination_arg}
    """
}

// Process 13: Run abricate
process ABRICATE {
    conda 'envs/abricate.yaml'
    
    input:
    path sample_dir
    
    output:
    path "${sample_dir}", emit: abricate_results
    
    script:
    def sample_name = sample_dir.name
    """
    mkdir -p ${sample_dir}/abricate
    
    # Check if sequence directory exists and has .fna files
    sequence_dir="${sample_dir}/sequence"
    if [ "${params.qc_filter}" = "true" ] && [ -d "${sample_dir}/sequence_filtered" ]; then
        sequence_dir="${sample_dir}/sequence_filtered"
    fi

    if [ -d "\${sequence_dir}" ] && [ -n "\$(find \${sequence_dir} -name "*.fna" -print -quit)" ]; then
        echo "Processing ${sample_name} with \$(find \${sequence_dir} -name "*.fna" | wc -l) .fna files"
        abricate --threads ${params.threads} --datadir ${params.db} \${sequence_dir}/*.fna > ${sample_dir}/abricate/ncbi_results.tab
        
        # Only create summary if results file is not empty
        if [ -s "${sample_dir}/abricate/ncbi_results.tab" ]; then
            abricate --summary ${sample_dir}/abricate/ncbi_results.tab > ${sample_dir}/abricate/ncbi_summary.tab
        else
            printf '#FILE\\tNUM_FOUND\\n' > ${sample_dir}/abricate/ncbi_summary.tab
            printf '#FILE\\tSEQUENCE\\tSTART\\tEND\\tGENE\\tCOVERAGE\\tCOVERAGE_MAP\\tGAPS\\t%%COVERAGE\\t%%IDENTITY\\tDATABASE\\tACCESSION\\tPRODUCT\\tRESISTANCE\\n' > ${sample_dir}/abricate/ncbi_results.tab
        fi
    else
        echo "Warning: No .fna files found in \${sequence_dir}/" >&2
        printf '#FILE\\tNUM_FOUND\\n' > ${sample_dir}/abricate/ncbi_summary.tab
        printf '#FILE\\tSEQUENCE\\tSTART\\tEND\\tGENE\\tCOVERAGE\\tCOVERAGE_MAP\\tGAPS\\t%%COVERAGE\\t%%IDENTITY\\tDATABASE\\tACCESSION\\tPRODUCT\\tRESISTANCE\\n' > ${sample_dir}/abricate/ncbi_results.tab
    fi
    """
}

// Process 14: Export PanR2-ready handoff directory
process EXPORT_PANR2_INPUTS {
    conda 'envs/fetchm.yaml'

    input:
    path sample_dir

    output:
    path "${sample_dir}", emit: panr2_inputs_results

    script:
    def versionReportsDir = params.outdir.toString().startsWith("/") ? "${params.outdir}/pipeline_versions" : "${launchDir}/${params.outdir}/pipeline_versions"
    def externalAmrfinderDir = params.amrfinderplus_dir ? launchPath(params.amrfinderplus_dir) : ""
    def largeDatasetFlag = truthyParam(params.large_dataset) ? "--large-dataset" : ""
    def skipHeavyPlotsFlag = effectiveSkipHeavyInteractivePlots() ? "--skip-heavy-interactive-plots" : ""
    def profileStack = workflow.profile ? workflow.profile.toString() : ""
    def commandLine = workflow.commandLine ? workflow.commandLine.toString() : "unknown"
    def sessionId = workflow.sessionId ? workflow.sessionId.toString() : "unknown"
    def runName = workflow.runName ? workflow.runName.toString() : "unknown"
    def containerImage = params.container_image ? params.container_image.toString() : ""
    """
    if [ "${params.export_panr2_inputs}" = "true" ]; then
        if [ -n "${externalAmrfinderDir}" ] && [ -d "${externalAmrfinderDir}" ]; then
            mkdir -p ${sample_dir}/amrfinderplus/tables
            cp -r "${externalAmrfinderDir}"/* ${sample_dir}/amrfinderplus/tables/ || true
        fi
        python ${baseDir}/scripts/export_panr2_inputs.py \\
            --sample-dir ${sample_dir} \\
            --versions-dir ${shellQuote(versionReportsDir)} \\
            --report-mode ${effectiveReportMode()} \\
            --max-features-heatmap ${effectiveMaxFeaturesHeatmap()} \\
            --max-features-network ${effectiveMaxFeaturesNetwork()} \\
            --max-metadata-columns ${effectiveMaxMetadataColumns()} \\
            --top-n-features-per-database ${effectiveTopNFeaturesPerDatabase()} \\
            --core-feature-threshold ${params.core_feature_threshold} \\
            --rare-feature-threshold ${params.rare_feature_threshold} \\
            --pipeline-version ${shellQuote(params.pipeline_version)} \\
            --repo-dir ${shellQuote(baseDir.toString())} \\
            --run-command ${shellQuote(commandLine)} \\
            --profile-stack ${shellQuote(profileStack)} \\
            --launch-dir ${shellQuote(launchDir.toString())} \\
            --pipeline-outdir ${shellQuote(params.outdir.toString())} \\
            --container-image ${shellQuote(containerImage)} \\
            --nextflow-session-id ${shellQuote(sessionId)} \\
            --nextflow-run-name ${shellQuote(runName)} \\
            ${largeDatasetFlag} \\
            ${skipHeavyPlotsFlag}
    fi
    """
}

// Process 15: Run panR2
process PANR {
    conda 'envs/panr2_comprehensive.yaml'
    
    input:
    path sample_dir
    
    output:
    path "${sample_dir}", emit: panr_results
    
    script:
    def sample_name = sample_dir.name
    """
    # Check if required directories exist
    if [ -d "${sample_dir}/metadata_output" ] && [ -d "${sample_dir}/abricate" ]; then
        if python ${baseDir}/scripts/should_run_panr.py --sample-dir ${sample_dir} --reason-file ${sample_dir}/panr_output/panr2_input_status.txt; then
            echo "Running panR2 for ${sample_name}"
            panr --ncbi-dir ${sample_dir}/metadata_output/ --abricate-dir ${sample_dir}/abricate/ --output-dir ${sample_dir}/ --format ${params.format}
        else
            echo "Skipping panR2 for ${sample_name}: \$(cat ${sample_dir}/panr_output/panr2_input_status.txt)" >&2
        fi
    else
        echo "Warning: Required directories not found for ${sample_name}" >&2
        if [ ! -d "${sample_dir}/metadata_output" ]; then
            echo "Missing: ${sample_dir}/metadata_output" >&2
        fi
        if [ ! -d "${sample_dir}/abricate" ]; then
            echo "Missing: ${sample_dir}/abricate" >&2
        fi
        # Create a placeholder file to indicate processing was attempted
        mkdir -p ${sample_dir}/panr_output
        echo "Processing failed: missing required directories" > ${sample_dir}/panr_output/error.log
    fi
    """
}

process PANR2_FEATURE_RUNNERS {
    conda 'envs/panr2_comprehensive.yaml'

    input:
    path sample_dir

    output:
    path "${sample_dir}", emit: panr2_feature_runner_results

    script:
    def panr2Dbs = effectivePanr2Dbs()
    def setupCmd = abricateSetupCommand(panr2Dbs, sample_dir)
    def checkm2DbPath = params.checkm2_db ? launchPath(params.checkm2_db) : ""
    def checkm2DownloadDir = params.checkm2_db_dir ? launchPath(params.checkm2_db_dir) : (params.outdir.toString().startsWith("/") ? "${params.outdir}/databases/checkm2" : "${launchDir}/${params.outdir}/databases/checkm2")
    def gtdbtkDataPath = params.gtdbtk_data_path ? launchPath(params.gtdbtk_data_path) : ""
    def isfinderDbFasta = params.isfinder_db_fasta ? launchPath(params.isfinder_db_fasta) : ""
    def mobsuiteDbPath = effectiveMobsuiteDbPath()
    def genomadDbPath = effectiveGenomadDbPath()
    def kaptiveDbPath = params.kaptive_db ? launchPath(params.kaptive_db) : ""
    """
    sequence_dir="${sample_dir}/sequence"
    if [ "${params.qc_filter}" = "true" ] && [ -d "${sample_dir}/sequence_filtered" ]; then
        sequence_dir="${sample_dir}/sequence_filtered"
    fi

    mkdir -p ${sample_dir}/panr_output ${sample_dir}/panr2_inputs/manifest
    if [ ! -d "${sample_dir}/metadata_output" ]; then
        echo "Processing failed: missing ${sample_dir}/metadata_output" > ${sample_dir}/panr_output/error.log
        exit 1
    fi
    if [ ! -d "\${sequence_dir}" ] || [ -z "\$(find "\${sequence_dir}" -name "*.fna" -print -quit)" ]; then
        echo "Processing failed: no .fna files found in \${sequence_dir}" > ${sample_dir}/panr_output/error.log
        exit 1
    fi

    echo "Preparing ABRicate databases for PanResistome-native PanR2 feature runners: ${panr2Dbs}"
    ${setupCmd}

    python ${baseDir}/scripts/database_setup_status.py \\
        --sample-dir ${sample_dir} \\
        --out ${sample_dir}/panr2_inputs/manifest/database_setup_status.tsv \\
        --analysis-profile ${analysisProfile()} \\
        --panr2-dbs ${panr2Dbs} \\
        --qc-filter ${params.qc_filter} \\
        --run-checkm2 ${params.run_checkm2} \\
        --checkm2-db ${shellQuote(checkm2DbPath)} \\
        --checkm2-db-dir ${shellQuote(checkm2DownloadDir)} \\
        --checkm2-auto-download-db ${params.checkm2_auto_download_db} \\
        --run-gtdbtk ${params.run_gtdbtk} \\
        --gtdbtk-data-path ${shellQuote(gtdbtkDataPath)} \\
        --run-quast ${params.run_quast} \\
        --run-ani ${params.run_ani} \\
        --run-mash ${params.run_mash} \\
        --run-panr2-comprehensive ${effectiveRunPanr2Comprehensive()} \\
        --panr2-update-abricate-db ${params.panr2_update_abricate_db} \\
        --run-integronfinder ${effectiveRunIntegronFinder()} \\
        --run-mlst ${effectiveRunMlst()} \\
        --run-mobileelementfinder ${effectiveRunMobileElementFinder()} \\
        --run-defensefinder ${params.panr2_run_defensefinder} \\
        --run-isfinder ${effectiveRunIsfinder()} \\
        --isfinder-db-fasta ${shellQuote(isfinderDbFasta)} \\
        --run-amrfinderplus ${params.run_amrfinderplus} \\
        --amrfinderplus-update-db ${params.amrfinderplus_update_db} \\
        --run-mobsuite ${params.run_mobsuite} \\
        --mobsuite-db ${shellQuote(mobsuiteDbPath)} \\
        --mobsuite-auto-init-db ${params.mobsuite_auto_init_db} \\
        --mobsuite-auto-init-taxa ${params.mobsuite_auto_init_taxa} \\
        --run-genomad ${params.run_genomad} \\
        --genomad-db ${shellQuote(genomadDbPath)} \\
        --genomad-auto-download-db ${params.genomad_auto_download_db} \\
        --run-kaptive ${params.run_kaptive} \\
        --kaptive-db ${shellQuote(kaptiveDbPath)} \\
        --strict

    python ${baseDir}/scripts/run_panr2_native_features.py \\
        --sample-dir ${sample_dir} \\
        --sequence-dir "\${sequence_dir}" \\
        --abricate-dbs ${panr2Dbs} \\
        --threads ${params.threads} \\
        --mode ${params.panr2_native_feature_runner_mode} \\
        --force ${params.panr2_force_tool_run} \\
        --run-integronfinder ${effectiveRunIntegronFinder()} \\
        --run-mlst ${effectiveRunMlst()} \\
        --run-mobileelementfinder ${effectiveRunMobileElementFinder()} \\
        --mobileelementfinder-allow-failure ${params.panr2_mobileelementfinder_allow_failure}
    """
}

process PANR2_COMPREHENSIVE {
    conda 'envs/panr2_comprehensive.yaml'

    input:
    path sample_dir

    output:
    path "${sample_dir}", emit: panr2_comprehensive_results

    script:
    def sample_name = sample_dir.name
    def optionalArgs = []
    if (params.panr2_label_max_length) {
        optionalArgs << "--label-max-length ${params.panr2_label_max_length}"
    }
    if (params.panr2_force_tool_run) {
        optionalArgs << "--force-tool-run"
    }
    if (params.panr2_run_defensefinder) {
        optionalArgs << "--run-defensefinder"
    }
    def optionalArgText = optionalArgs.join(' ')
    def externalFeatureArgs = []
    if (params.defensefinder_dir) {
        externalFeatureArgs << "--defensefinder-dir ${launchPath(params.defensefinder_dir)}"
    }
    if (params.isfinder_dir) {
        externalFeatureArgs << "--isfinder-dir ${launchPath(params.isfinder_dir)}"
    }
    if (params.mobsuite_dir) {
        externalFeatureArgs << "--mobsuite-dir ${launchPath(params.mobsuite_dir)}"
    }
    if (params.prophage_dir) {
        externalFeatureArgs << "--prophage-dir ${launchPath(params.prophage_dir)}"
    }
    if (params.kleborate_dir) {
        externalFeatureArgs << "--kleborate-dir ${launchPath(params.kleborate_dir)}"
    }
    if (params.kaptive_dir) {
        externalFeatureArgs << "--kaptive-dir ${launchPath(params.kaptive_dir)}"
    }
    if (params.ectyper_dir) {
        externalFeatureArgs << "--ectyper-dir ${launchPath(params.ectyper_dir)}"
    }
    if (params.serotypefinder_dir) {
        externalFeatureArgs << "--serotypefinder-dir ${launchPath(params.serotypefinder_dir)}"
    }
    if (params.sccmecfinder_dir) {
        externalFeatureArgs << "--sccmecfinder-dir ${launchPath(params.sccmecfinder_dir)}"
    }
    def externalFeatureArgText = externalFeatureArgs.join(' ')
    def panr2Dbs = effectivePanr2Dbs()
    def setupCmd = abricateSetupCommand(panr2Dbs, sample_dir)
    def configuredSampleMap = params.panr2_sample_map ? launchPath(params.panr2_sample_map) : ""
    def abricateRunFlag = params.panr2_native_feature_runners ? "" : "--run-abricate"
    def mobileElementFinderFlag = (!params.panr2_native_feature_runners && effectiveRunMobileElementFinder()) ? "--run-mobileelementfinder" : ""
    def integronFinderFlag = (!params.panr2_native_feature_runners && effectiveRunIntegronFinder()) ? "--run-integronfinder" : ""
    def mlstFlag = (!params.panr2_native_feature_runners && effectiveRunMlst()) ? "--run-mlst" : ""
    def checkm2DbPath = params.checkm2_db ? launchPath(params.checkm2_db) : ""
    def checkm2DownloadDir = params.checkm2_db_dir ? launchPath(params.checkm2_db_dir) : (params.outdir.toString().startsWith("/") ? "${params.outdir}/databases/checkm2" : "${launchDir}/${params.outdir}/databases/checkm2")
    def gtdbtkDataPath = params.gtdbtk_data_path ? launchPath(params.gtdbtk_data_path) : ""
    def isfinderDbFasta = params.isfinder_db_fasta ? launchPath(params.isfinder_db_fasta) : ""
    def mobsuiteDbPath = effectiveMobsuiteDbPath()
    def genomadDbPath = effectiveGenomadDbPath()
    def kaptiveDbPath = params.kaptive_db ? launchPath(params.kaptive_db) : ""
    """
    sequence_dir="${sample_dir}/sequence"
    if [ "${params.qc_filter}" = "true" ] && [ -d "${sample_dir}/sequence_filtered" ]; then
        sequence_dir="${sample_dir}/sequence_filtered"
    fi

    mkdir -p ${sample_dir}/panr_output
    if [ ! -d "${sample_dir}/metadata_output" ]; then
        echo "Processing failed: missing ${sample_dir}/metadata_output" > ${sample_dir}/panr_output/error.log
        exit 1
    fi
    if [ ! -d "\${sequence_dir}" ] || [ -z "\$(find "\${sequence_dir}" -name "*.fna" -print -quit)" ]; then
        echo "Processing failed: no .fna files found in \${sequence_dir}" > ${sample_dir}/panr_output/error.log
        exit 1
    fi

    if [ "${params.panr2_native_feature_runners}" != "true" ]; then
        echo "Preparing ABRicate databases for comprehensive PanR2 integrated analysis: ${panr2Dbs}"
        ${setupCmd}

        python ${baseDir}/scripts/database_setup_status.py \\
            --sample-dir ${sample_dir} \\
            --out ${sample_dir}/panr2_inputs/manifest/database_setup_status.tsv \\
            --analysis-profile ${analysisProfile()} \\
            --panr2-dbs ${panr2Dbs} \\
            --qc-filter ${params.qc_filter} \\
            --run-checkm2 ${params.run_checkm2} \\
            --checkm2-db ${shellQuote(checkm2DbPath)} \\
            --checkm2-db-dir ${shellQuote(checkm2DownloadDir)} \\
            --checkm2-auto-download-db ${params.checkm2_auto_download_db} \\
            --run-gtdbtk ${params.run_gtdbtk} \\
            --gtdbtk-data-path ${shellQuote(gtdbtkDataPath)} \\
            --run-quast ${params.run_quast} \\
            --run-ani ${params.run_ani} \\
            --run-mash ${params.run_mash} \\
            --run-panr2-comprehensive ${effectiveRunPanr2Comprehensive()} \\
            --panr2-update-abricate-db ${params.panr2_update_abricate_db} \\
            --run-integronfinder ${effectiveRunIntegronFinder()} \\
            --run-mlst ${effectiveRunMlst()} \\
            --run-mobileelementfinder ${effectiveRunMobileElementFinder()} \\
            --run-defensefinder ${params.panr2_run_defensefinder} \\
            --run-isfinder ${effectiveRunIsfinder()} \\
            --isfinder-db-fasta ${shellQuote(isfinderDbFasta)} \\
            --run-amrfinderplus ${params.run_amrfinderplus} \\
            --amrfinderplus-update-db ${params.amrfinderplus_update_db} \\
            --run-mobsuite ${params.run_mobsuite} \\
            --mobsuite-db ${shellQuote(mobsuiteDbPath)} \\
            --mobsuite-auto-init-db ${params.mobsuite_auto_init_db} \\
            --mobsuite-auto-init-taxa ${params.mobsuite_auto_init_taxa} \\
            --run-genomad ${params.run_genomad} \\
            --genomad-db ${shellQuote(genomadDbPath)} \\
            --genomad-auto-download-db ${params.genomad_auto_download_db} \\
            --run-kaptive ${params.run_kaptive} \\
            --kaptive-db ${shellQuote(kaptiveDbPath)} \\
            --strict
    fi

    sample_map_arg=""
    if [ -n "${configuredSampleMap}" ]; then
        sample_map_arg="--sample-map ${configuredSampleMap}"
    elif [ -f "${sample_dir}/metadata_output/sample_map.csv" ]; then
        sample_map_arg="--sample-map ${sample_dir}/metadata_output/sample_map.csv"
    fi

    extra_feature_args="${externalFeatureArgText}"
    if [ "${params.panr2_native_feature_runners}" = "true" ]; then
        for spec in \\
            "abricate --abricate-dir ${sample_dir}/tool_results/abricate/ncbi table" \\
            "vfdb --vfdb-dir ${sample_dir}/tool_results/abricate/vfdb table" \\
            "plasmidfinder --plasmidfinder-dir ${sample_dir}/tool_results/abricate/plasmidfinder table" \\
            "integronfinder --integronfinder-dir ${sample_dir}/tool_results/integronfinder/panr2_inputs table" \\
            "mobileelementfinder --mobileelementfinder-dir ${sample_dir}/tool_results/mobileelementfinder/panr2_inputs table" \\
            "mlst --mlst-dir ${sample_dir}/tool_results/mlst/raw mlst"; do
            name=\$(echo "\${spec}" | cut -d' ' -f1)
            option=\$(echo "\${spec}" | cut -d' ' -f2)
            directory=\$(echo "\${spec}" | cut -d' ' -f3)
            mode=\$(echo "\${spec}" | cut -d' ' -f4)
            if [ "\${mode}" = "mlst" ]; then
                if [ -d "\${directory}" ] && find "\${directory}" -type f \\( -name "*.tsv" -o -name "*.tab" -o -name "*.csv" \\) -print -quit | grep -q .; then
                    extra_feature_args="\${extra_feature_args} \${option} \${directory}"
                    echo "Passing PanResistome-native \${name} tables to PanR2: \${directory}"
                fi
            elif [ -d "\${directory}" ] && python ${baseDir}/scripts/has_feature_table_rows.py "\${directory}"; then
                extra_feature_args="\${extra_feature_args} \${option} \${directory}"
                echo "Passing PanResistome-native \${name} tables to PanR2: \${directory}"
            fi
        done
    fi
    for spec in \\
        "isfinder --isfinder-dir ${sample_dir}/isfinder/tables" \\
        "mobsuite --mobsuite-dir ${sample_dir}/mobsuite/tables" \\
        "kleborate --kleborate-dir ${sample_dir}/kleborate/tables" \\
        "kaptive --kaptive-dir ${sample_dir}/kaptive/tables" \\
        "ectyper --ectyper-dir ${sample_dir}/ectyper/tables" \\
        "prophage --prophage-dir ${sample_dir}/prophage/tables"; do
        name=\$(echo "\${spec}" | cut -d' ' -f1)
        option=\$(echo "\${spec}" | cut -d' ' -f2)
        directory=\$(echo "\${spec}" | cut -d' ' -f3)
        if [ -d "\${directory}" ] && python ${baseDir}/scripts/has_feature_table_rows.py "\${directory}"; then
            extra_feature_args="\${extra_feature_args} \${option} \${directory}"
            echo "Passing \${name} tables to PanR2: \${directory}"
        fi
    done
    echo "Running comprehensive PanR2 for ${sample_name}"
    panr \\
        --ncbi-dir ${sample_dir}/metadata_output/ \\
        --sequence-dir "\${sequence_dir}" \\
        --output-dir ${sample_dir}/ \\
        ${abricateRunFlag} \\
        ${mobileElementFinderFlag} \\
        ${integronFinderFlag} \\
        ${mlstFlag} \\
        --format ${params.format} \\
        --abricate-dbs ${panr2Dbs} \\
        --min-identity ${params.panr2_min_identity} \\
        --plot-style ${effectivePanr2PlotStyle()} \\
        --cross-database-max-features ${effectiveMaxFeaturesNetwork()} \\
        --mobileelementfinder-threads ${params.threads} \\
        --integronfinder-threads ${params.threads} \\
        ${optionalArgText} \\
        \${sample_map_arg} \\
        \${extra_feature_args}
    """
}

// Process 11: Collect final results
process COLLECT_RESULTS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path sample_dir

    output:
    path "${sample_dir.name}", emit: final_results

    script:
    """
    mkdir -p ${params.outdir}/${sample_dir.name}
    cp -r ${sample_dir}/* ${params.outdir}/${sample_dir.name}/
    """
}

// Workflow
workflow {
    // Create input channel
    input_ch = Channel.fromPath(params.input, checkIfExists: true)

    // Capture reproducibility metadata for each Conda environment
    if (params.capture_versions) {
        FETCHM_ENV_VERSIONS()
        ABRICATE_ENV_VERSIONS()
        version_reports = FETCHM_ENV_VERSIONS.out.fetchm_versions
            .mix(ABRICATE_ENV_VERSIONS.out.abricate_versions)
        if (effectiveRunPanr2Comprehensive()) {
            PANR2_COMPREHENSIVE_ENV_VERSIONS()
            version_reports = version_reports.mix(PANR2_COMPREHENSIVE_ENV_VERSIONS.out.panr2_comprehensive_versions)
        }
        if (params.run_amrfinderplus) {
            AMRFINDERPLUS_ENV_VERSIONS()
            version_reports = version_reports.mix(AMRFINDERPLUS_ENV_VERSIONS.out.amrfinderplus_versions)
        }
        if (params.run_checkm2) {
            CHECKM2_ENV_VERSIONS()
            version_reports = version_reports.mix(CHECKM2_ENV_VERSIONS.out.checkm2_versions)
        }
        if (params.run_gtdbtk) {
            GTDBTK_ENV_VERSIONS()
            version_reports = version_reports.mix(GTDBTK_ENV_VERSIONS.out.gtdbtk_versions)
        }
        if (params.run_ani) {
            ANI_ENV_VERSIONS()
            version_reports = version_reports.mix(ANI_ENV_VERSIONS.out.ani_versions)
        }
        if (params.run_quast) {
            QUAST_ENV_VERSIONS()
            version_reports = version_reports.mix(QUAST_ENV_VERSIONS.out.quast_versions)
        }
        if (params.run_mash) {
            MASH_ENV_VERSIONS()
            version_reports = version_reports.mix(MASH_ENV_VERSIONS.out.mash_versions)
        }
        if (params.run_mobsuite) {
            MOBSUITE_ENV_VERSIONS()
            version_reports = version_reports.mix(MOBSUITE_ENV_VERSIONS.out.mobsuite_versions)
        }
        if (params.run_genomad) {
            GENOMAD_ENV_VERSIONS()
            version_reports = version_reports.mix(GENOMAD_ENV_VERSIONS.out.genomad_versions)
        }
        if (params.run_organism_specific_typing || params.run_kleborate || params.run_kaptive || params.run_ectyper) {
            ORGANISM_TYPING_ENV_VERSIONS()
            version_reports = version_reports.mix(ORGANISM_TYPING_ENV_VERSIONS.out.organism_typing_versions)
        }
        version_reports.view { version_file -> "Version report saved: ${version_file}" }
    }
    
    if (params.local_samples) {
        sample_dirs = Channel.fromPath("${params.local_samples}/*", type: 'dir', checkIfExists: true)
            .filter { sample_dir -> !(sample_dir.name in ['pipeline_versions', 'work', 'report', 'trace']) }
    } else {
        // Run FetchM2 metadata/download adapter
        FETCHM(input_ch)

        // Create channel for sample directories
        sample_dirs = FETCHM.out.fetchm_results
            .map { results_dir ->
                results_dir.listFiles().findAll { it.isDirectory() }
            }
            .flatten()
    }
    
    // Generate sequence QC stats and enriched metadata
    SEQUENCE_QC(sample_dirs)

    qc_ready_ch = SEQUENCE_QC.out.qc_results

    if (params.run_checkm2) {
        // Add CheckM2 completeness/contamination QC
        CHECKM2_QC(SEQUENCE_QC.out.qc_results)
        qc_ready_ch = CHECKM2_QC.out.checkm2_results
    }

    if (params.run_gtdbtk) {
        // Add GTDB-Tk taxonomy match QC
        GTDBTK_QC(qc_ready_ch)
        qc_ready_ch = GTDBTK_QC.out.gtdbtk_results
    }

    if (params.run_quast) {
        QUAST_QC(qc_ready_ch)
        qc_ready_ch = QUAST_QC.out.quast_results
    }

    if (params.run_ani) {
        ANI_ANALYSIS(qc_ready_ch)
        qc_ready_ch = ANI_ANALYSIS.out.ani_results
    }

    if (params.run_mash) {
        MASH_PRESCREEN(qc_ready_ch)
        qc_ready_ch = MASH_PRESCREEN.out.mash_results
    }

    COMBINED_QC(qc_ready_ch)
    qc_ready_ch = COMBINED_QC.out.combined_qc_results

    if (effectiveStopAfterQc()) {
        // Collect QC-only results to output directory
        COLLECT_RESULTS(qc_ready_ch)
    } else {
        if (params.run_amrfinderplus) {
            AMRFINDERPLUS_ANALYSIS(qc_ready_ch)
            qc_ready_ch = AMRFINDERPLUS_ANALYSIS.out.amrfinderplus_results
        }

        if (effectiveRunIsfinder()) {
            ISFINDER_BLAST(qc_ready_ch)
            qc_ready_ch = ISFINDER_BLAST.out.isfinder_results
        }

        if (params.run_mobsuite) {
            MOBSUITE_ANALYSIS(qc_ready_ch)
            qc_ready_ch = MOBSUITE_ANALYSIS.out.mobsuite_results
        }

        if (params.run_genomad) {
            GENOMAD_PROPHAGE(qc_ready_ch)
            qc_ready_ch = GENOMAD_PROPHAGE.out.genomad_results
        }

        if (params.run_organism_specific_typing || params.run_kleborate || params.run_kaptive || params.run_ectyper) {
            ORGANISM_SPECIFIC_TYPING(qc_ready_ch)
            qc_ready_ch = ORGANISM_SPECIFIC_TYPING.out.organism_typing_results
        }

        if (effectiveRunPanr2Comprehensive()) {
            if (params.panr2_native_feature_runners) {
                PANR2_FEATURE_RUNNERS(qc_ready_ch)
                qc_ready_ch = PANR2_FEATURE_RUNNERS.out.panr2_feature_runner_results
            }
            PANR2_COMPREHENSIVE(qc_ready_ch)
            EXPORT_PANR2_INPUTS(PANR2_COMPREHENSIVE.out.panr2_comprehensive_results)
            COLLECT_RESULTS(EXPORT_PANR2_INPUTS.out.panr2_inputs_results)
        } else {
            if (params.run_abricate) {
                // Run abricate on each sample directory
                ABRICATE(qc_ready_ch)

                EXPORT_PANR2_INPUTS(ABRICATE.out.abricate_results)

                // Run panR on each sample directory after abricate
                PANR(EXPORT_PANR2_INPUTS.out.panr2_inputs_results)

                // Collect final results to output directory
                COLLECT_RESULTS(PANR.out.panr_results)
            } else {
                // Optional-runner/table-input validation path: export current standardized outputs without legacy ABRicate/PanR.
                EXPORT_PANR2_INPUTS(qc_ready_ch)
                COLLECT_RESULTS(EXPORT_PANR2_INPUTS.out.panr2_inputs_results)
            }
        }
    }
    
    // Display completion message
    COLLECT_RESULTS.out.final_results.collect().view { "Pipeline completed. Results saved to: ${params.outdir}" }
}

workflow.onComplete {
    def tracePath = "${launchDir}/pipeline_trace.txt"
    def traceFile = file(tracePath)
    if (!traceFile.exists() || traceFile.lastModified() < pipelineStartMillis) {
        log.info "Runtime/resource summary skipped because no current-run Nextflow trace was found."
        return
    }
    def summaryPath = "${params.outdir}/pipeline_runtime_summary.tsv"
    def taskPath = "${params.outdir}/pipeline_runtime_tasks.tsv"
    def pythonExe = "python3"
    try {
        def pythonCheck = [pythonExe, "--version"].execute()
        pythonCheck.waitFor()
        if (pythonCheck.exitValue() != 0) {
            pythonExe = "python"
        }
    } catch (Exception ignored) {
        pythonExe = "python"
    }
    def command = [
        pythonExe,
        "${baseDir}/scripts/summarize_nextflow_trace.py",
        "--trace",
        tracePath,
        "--out",
        summaryPath,
        "--tasks-out",
        taskPath,
    ]
    try {
        def process = command.execute()
        process.waitFor()
        if (process.exitValue() == 0) {
            log.info "Runtime/resource summary saved to: ${summaryPath}"
        } else {
            log.warn "Runtime/resource summary failed: ${process.err.text}"
        }
    } catch (Exception error) {
        log.warn "Runtime/resource summary could not be generated: ${error.message}"
    }
}
