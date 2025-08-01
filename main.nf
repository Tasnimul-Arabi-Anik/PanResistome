#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Parameters
params.checkm = null
params.ani = 'all'
params.sleep = 0.5
params.host = []
params.year = []
params.country = []
params.cont = []
params.subcont = []

params.input = "test.tsv"
params.outdir = "results"
params.threads = 24
params.db = "$baseDir/db"
params.help = false
params.format = "png"
params.genep   = null
params.nseq    = null


// Help message
def helpMessage() {
    log.info """
    📦 PanResistome Pipeline
    ------------------------
    A Nextflow pipeline for downloading genomic data, identifying resistance genes, and visualizing pangenome resistome profiles.

    ▶️ Usage:
      nextflow run main.nf --input <input.tsv> --outdir <output_dir> [options]

    ✅ Required arguments:
      --input            Input TSV file listing genome accessions
      --outdir           Output directory for results

    ⚙️ Optional arguments for fetchM:
      --checkm           Minimum CheckM completeness threshold (e.g. 90. Default: null)
      --ani              ANI filter status (Choices: OK, Inconclusive, Failed, all. Default: all)
      --sleep            Time to wait between fetch requests (default: 0.5s)

        🧬 Instead of global resistance analysis, you may do specific analysis by providing: 
      --host             Host species (e.g. "Homo sapiens" "Bos taurus")
      --year             Filter by year or year range (e.g. "2015" or "2015-2023")
      --country          Country filter (e.g. "Bangladesh" "USA")
      --cont             Continent filter (e.g. "Asia", "Africa")
      --subcont          Subcontinent filter (e.g. "Southern Asia")

    🧬 Optional arguments for PanR2:
      --genep            Minimum % gene presence to include in heatmap (float)
      --nseq             Minimum number of sequences per group in heatmaps (int)
      --format           Output format for plots (tiff, svg, png, pdf) [default: png]

    🔧 Other options:
      --threads          Number of threads for abricate [default: 24]
      --db               Directory containing abricate databases [default: ./db]
      --help             Show this help message and exit

    Example:
       nextflow run main1.nf --input test.tsv --outdir results -profile conda --threads 24 
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

// Process 1: Run fetchM
process FETCHM {
    conda 'envs/fetchm.yaml'
    
    input:
    path input_file
    
    output:
    path "results", emit: fetchm_results
    
    script:
    """
    fetchM \\
        --input ${input_file} \\
        --outdir ${params.outdir}/ \\
        ${params.checkm ? "--checkm ${params.checkm}" : ""} \\
        --ani ${params.ani} \\
        --sleep ${params.sleep} \\
        --seq
        ${params.host ? "--host ${params.host.join(' ')}" : ""} \\
        ${params.year ? "--year ${params.year.join(' ')}" : ""} \\
        ${params.country ? "--country ${params.country.join(' ')}" : ""} \\
        ${params.cont ? "--cont ${params.cont.join(' ')}" : ""} \\
        ${params.subcont ? "--subcont ${params.subcont.join(' ')}" : ""}
    """
}

// Process 2: Run abricate
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
    if [ -d "${sample_dir}/sequence" ] && [ -n "\$(find ${sample_dir}/sequence -name "*.fna" -print -quit)" ]; then
        echo "Processing ${sample_name} with \$(find ${sample_dir}/sequence -name "*.fna" | wc -l) .fna files"
        abricate --threads ${params.threads} --datadir ${params.db} ${sample_dir}/sequence/*.fna > ${sample_dir}/abricate/ncbi_results.tab
        
        # Only create summary if results file is not empty
        if [ -s "${sample_dir}/abricate/ncbi_results.tab" ]; then
            abricate --summary ${sample_dir}/abricate/ncbi_results.tab > ${sample_dir}/abricate/ncbi_summary.tab
        else
            echo "No results found for ${sample_name}" > ${sample_dir}/abricate/ncbi_summary.tab
        fi
    else
        echo "Warning: No .fna files found in ${sample_dir}/sequence/" >&2
        touch ${sample_dir}/abricate/ncbi_results.tab
        echo "No .fna files found" > ${sample_dir}/abricate/ncbi_summary.tab
    fi
    """
}

// Process 3: Run panR2
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
        echo "Running panR2 for ${sample_name}"
        panr --ncbi-dir ${sample_dir}/metadata_output/ --abricate-dir ${sample_dir}/abricate/ --output-dir ${sample_dir}/ --format ${params.format}
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

// Process 4: Collect final results
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
    
    // Run fetchM
    FETCHM(input_ch)
    
    // Create channel for sample directories
    sample_dirs = FETCHM.out.fetchm_results
        .map { results_dir -> 
            results_dir.listFiles().findAll { it.isDirectory() }
        }
        .flatten()
    
    // Run abricate on each sample directory
    ABRICATE(sample_dirs)
    
    // Run panR on each sample directory after abricate
    PANR(ABRICATE.out.abricate_results)
    
    // Collect final results to output directory
    COLLECT_RESULTS(PANR.out.panr_results)
    
    // Display completion message
    COLLECT_RESULTS.out.final_results.collect().view { "Pipeline completed. Results saved to: ${params.outdir}" }
}
