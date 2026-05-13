# GHCR Docker 100-record validation results

Date: 2026-05-13

Validation target: confirm that the public GHCR image can run the 100-record
`Klebsiella pneumoniae` large-mode biological workflow through the Docker
profile, not only the locally built image.

Image:

```text
ghcr.io/tasnimul-arabi-anik/panresistome:experimental
```

Command:

```bash
nextflow run main.nf \
  --input validation/klebsiella_pneumoniae_100/ncbi_dataset.tsv \
  --outdir validation_runs/klebsiella_100_ghcr_docker_large \
  -profile docker,large \
  --container_image ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_checkm2 false \
  --run_quast false \
  --run_ani false \
  --run_mash true \
  --run_amrfinderplus false \
  --panr2_native_feature_runners true \
  --panr2_native_feature_runner_mode parallel \
  --threads 8 \
  --fetchm2_download_workers 2
```

The command was launched from a shell where `anik` was run with active group
`docker`, matching the non-sudo Docker permission model expected after a fresh
login session. The normal interactive shell still requires a new login before
its supplementary group list includes `docker`.

## Result

```text
Status: PASS
Nextflow processes: 12/12 succeeded
Duration: 21m35s
CPU hours: 4.7
Output: validation_runs/klebsiella_100_ghcr_docker_large
```

Input and QC:

```text
Input records: 100
Downloaded/analyzed genomes: 99
Failed download: GCF_055382775.1
Sequence QC PASS: 99
QC master PASS: 99
QC master FAIL: 1
```

The single failed accession is the same public-record download failure observed
in earlier Klebsiella validations. The pipeline recorded it and continued with
the remaining 99 genomes.

Feature-contract validation:

```text
feature_files_checked=5
feature_rows=11488
databases_seen=amr,integronfinder,mlst,plasmidfinder,vfdb
metadata_accessions=100
unmatched_feature_rows=0
invalid_feature_rows=0
duplicate_feature_rows=0
```

Feature rows by database:

```text
amr.features.tsv:             1234
vfdb.features.tsv:            7821
plasmidfinder.features.tsv:    393
integronfinder.features.tsv:   418
mlst.features.tsv:            1622
all_features.tsv:            11488
```

Database and runner audits:

```text
ABRicate setup/update: PASS for ncbi, vfdb, plasmidfinder
ABRicate raw tables:   297 expected, 297 observed
ABRicate samples:      99 processed, 0 failed
IntegronFinder:        99 expected, 99 observed; 99 processed, 0 failed
MLST:                  99 processed, 0 failed
MobileElementFinder:   SKIPPED, not requested
```

Large-mode report controls:

```text
large_dataset=true
report_mode=compact
max_features_heatmap=150
max_features_network=150
max_metadata_columns=20
top_n_features_per_database=50
skip_heavy_interactive_plots=true
```

Runtime/resource summary:

```text
FETCHM:                 7.27m
PANR2_FEATURE_RUNNERS: 11.57m
PANR2_COMPREHENSIVE:    1.93m
max peak RSS:           2.2 GiB
max peak vmem:         33.3 GiB
```

## Interpretation

This closes the previous deployment validation gap:

```text
Docker local image, 100-record biological validation: PASS
Docker GHCR image, 2-genome geNomad validation: PASS
Docker GHCR image, 100-record large-mode validation: PASS
Singularity GHCR image, 100-record large-mode validation: PASS
Apptainer GHCR image validation: pending
```

The run also confirms that default ABRicate setup/update automation works in
the GHCR Docker route for the default comprehensive databases.

