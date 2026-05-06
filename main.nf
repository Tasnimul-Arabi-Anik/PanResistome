#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params.pipeline_version = '0.2.1-dev'
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
params.run_mash = false
params.representative_only = false
params.export_panr2_inputs = true


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
      --run_mash               Enable Mash sketch/distance pre-screen [default: false]
      --representative_only    Keep one representative per near-duplicate ANI cluster for PanR2 when --qc_filter true [default: false]
      --export_panr2_inputs    Export standardized panr2_inputs handoff directory [default: true]
      --local_samples          Optional directory of prebuilt sample folders for offline tests
      --run_checkm2            Enable CheckM2 QC [default: true]

    🔧 Other options:
  --threads          Number of threads for CheckM2, GTDB-Tk, and abricate [default: 8]
      --db               Directory containing abricate databases [default: ./db]
      --help             Show this help message and exit

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

// Process 1: Capture versions from the FetchM2/PanR2/seqkit environment
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

for package in ("fetchm2", "fetchM", "PanR2"):
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
    def checkm2_db_arg = params.checkm2_db ? "--database_path ${params.checkm2_db}" : ""
    def checkm2_lowmem_arg = params.checkm2_lowmem ? "--lowmem" : ""
    """
    mkdir -p ${sample_dir}/checkm2

    sequence_dir="${sample_dir}/sequence"
    if [ "${params.qc_filter}" = "true" ] && [ -d "${sample_dir}/sequence_filtered" ]; then
        sequence_dir="${sample_dir}/sequence_filtered"
    fi

    if [ -d "\${sequence_dir}" ] && [ -n "\$(find \${sequence_dir} -name "*.fna" -print -quit)" ]; then
        checkm2 predict \\
            --threads ${params.threads} \\
            --input \${sequence_dir} \\
            --output-directory ${sample_dir}/checkm2 \\
            -x fna --force ${checkm2_lowmem_arg} ${checkm2_db_arg}
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
    """
    mkdir -p ${sample_dir}/ani/analysis
    sequence_dir="${sample_dir}/sequence"
    if [ "${params.qc_filter}" = "true" ] && [ -d "${sample_dir}/sequence_filtered" ]; then
        sequence_dir="${sample_dir}/sequence_filtered"
    fi
    find \${sequence_dir} -name "*.fna" | sort > ${sample_dir}/ani/genomes.list || true
    genome_count=\$(wc -l < ${sample_dir}/ani/genomes.list || echo 0)
    if [ "\${genome_count}" -ge 2 ]; then
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
            --tool ${params.ani_tool} \\
            --duplicate-threshold ${params.ani_duplicate_threshold} \\
            --species-threshold ${params.ani_species_threshold}
    else
        printf "query\\treference\\tani\\tfragments_mapped\\tfragments_total\\n" > ${sample_dir}/ani/fastani_pairs.tsv
        python ${baseDir}/scripts/ani_summary.py --sample-dir ${sample_dir} --pairs ${sample_dir}/ani/fastani_pairs.tsv --tool ${params.ani_tool}
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
    python ${baseDir}/scripts/mash_summary.py --sample-dir ${sample_dir} --dist ${sample_dir}/mash/mash_dist.tsv
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
    """
    if [ "${params.export_panr2_inputs}" = "true" ]; then
        python ${baseDir}/scripts/export_panr2_inputs.py --sample-dir ${sample_dir} --versions-dir ${shellQuote(versionReportsDir)}
    fi
    """
}

// Process 15: Run panR2
process PANR {
    conda 'envs/fetchm.yaml'
    
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

    if (params.stop_after_qc) {
        // Collect QC-only results to output directory
        COLLECT_RESULTS(qc_ready_ch)
    } else {
        // Run abricate on each sample directory
        ABRICATE(qc_ready_ch)

        EXPORT_PANR2_INPUTS(ABRICATE.out.abricate_results)
        
        // Run panR on each sample directory after abricate
        PANR(EXPORT_PANR2_INPUTS.out.panr2_inputs_results)
        
        // Collect final results to output directory
        COLLECT_RESULTS(PANR.out.panr_results)
    }
    
    // Display completion message
    COLLECT_RESULTS.out.final_results.collect().view { "Pipeline completed. Results saved to: ${params.outdir}" }
}
