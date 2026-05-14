# v0.6.0 ANI-enabled Singularity/GHCR lineage validation

Date: 2026-05-14

This validation exercises the v0.5.0 lineage-aware interpretation layer with
real pairwise ANI enabled, using the public GHCR image through Singularity CE.
It also confirms that native IntegronFinder handoff rows are preserved in the
strict PanR2 feature contract after the v0.6.0 native-handoff export fix.

## Command

The successful patched validation used the cached GHCR image and a longer
Singularity pull timeout:

```bash
NXF_DISABLE_CHECK_LATEST=true \
SINGULARITY_CACHEDIR=/tmp/panresistome_singularity_cache \
SINGULARITY_TMPDIR=/tmp/panresistome_singularity_tmp \
NXF_SINGULARITY_CACHEDIR=/tmp/panresistome_singularity_cache \
nextflow -c /tmp/panresistome_singularity_long_pull.config run main.nf \
  --input validation/klebsiella_pneumoniae_10/ncbi_dataset.tsv \
  --outdir validation_runs/klebsiella_10_ani_lineage_singularity_fixed \
  -profile singularity,large \
  --container_image docker://ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_checkm2 false \
  --run_gtdbtk false \
  --run_quast false \
  --run_ani true \
  --ani_tool skani \
  --ani_large_run_strategy all \
  --run_mash true \
  --run_amrfinderplus false \
  --panr2_native_feature_runners true \
  --panr2_native_feature_runner_mode parallel \
  --threads 8 \
  --fetchm2_download_workers 2 \
  -w /tmp/panresistome_v06_ani_fixed_work
```

The temporary config contained:

```nextflow
singularity.pullTimeout = '2h'
```

The first uncached GHCR-to-SIF conversion exceeded the default Nextflow
20-minute Singularity pull timeout. Retrying with a 2-hour timeout and persistent
cache completed successfully.

## Result

- Processes: 14/14 succeeded.
- Runtime after image cache was available: 3m44s.
- CPU hours: 0.7.
- Input records: 10.
- Genomes analyzed: 10.
- CheckM2: disabled for speed.
- GTDB-Tk: disabled.
- QUAST: disabled.
- ANI: enabled with skani.
- Mash: enabled.
- AMRFinderPlus: disabled.
- Native PanR2 feature runners: enabled.
- Native runner mode: parallel.

## ANI and lineage evidence

`ani/analysis/ani_run_status.tsv`:

```text
tool	genome_count	estimated_comparisons	strategy	max_all_vs_all_genomes	large_dataset	decision	status	message
skani	10	100	all	200	true	run_all_vs_all	PASS	Running all-vs-all ANI.
```

Generated ANI outputs included:

- `ani/analysis/ani_matrix.csv`
- `ani/analysis/pairwise_ani_long.csv`
- `ani/analysis/duplicate_clusters.csv`
- `ani/analysis/panr2_ani_summary.csv`

Generated lineage-aware PanR2 outputs included:

- `panr2_inputs/metadata_feature_analysis/lineage_summary.tsv`
- `panr2_inputs/metadata_feature_analysis/lineage_feature_burden.tsv`
- `panr2_inputs/metadata_feature_analysis/lineage_metadata_confounding.tsv`
- `panr2_inputs/metadata_feature_analysis/lineage_adjusted_warnings.tsv`
- `panr2_inputs/report/lineage_context.html`

`lineage_summary.tsv` contains ANI clusters for all 10 assemblies.

## Feature contract

`panr2_inputs/manifest/schema_validation_summary.txt`:

```text
feature_files_checked=5
feature_rows=1248
databases_seen=amr,integronfinder,mlst,plasmidfinder,vfdb
samples_seen=20
metadata_accessions=10
unmatched_feature_rows=0
invalid_feature_rows=0
duplicate_feature_rows=0
```

`samples_seen` is larger than `metadata_accessions` because the validator sees
both assembly accessions and normalized sample IDs in feature rows. The release
gate is that unmatched feature rows remained zero.

Feature rows by database:

| Database | Rows |
| --- | ---: |
| AMR | 132 |
| IntegronFinder | 44 |
| MLST | 170 |
| PlasmidFinder | 42 |
| VFDB | 860 |

## Native IntegronFinder handoff fix

The first ANI-enabled run exposed a v0.6.0 bug: native IntegronFinder generated
44 real handoff rows under `tool_results/integronfinder/panr2_inputs/`, and the
native runner audit reported `PASS`, but the strict exported
`integronfinder.features.tsv` was header-only because feature-row discovery
skipped all paths inside `panr2_inputs` directories.

The fix adds a targeted exception for native runner handoff result tables under:

- `tool_results/integronfinder/panr2_inputs/*_results.*`
- `tool_results/mobileelementfinder/panr2_inputs/*_results.*`

After the fix, the clean end-to-end validation produced 44 standardized
IntegronFinder rows and `feature_completeness_audit.tsv` reported
IntegronFinder `PASS`.

## Report outputs

The validation generated:

- `panr2_inputs/report/panr2_handoff_index.html`
- `panr2_inputs/report/lineage_context.html`
- `panr2_inputs/report/diversity_summary.html`
- `panr2_inputs/report/statistical_summary.html`

## Interpretation

This validation confirms that compact Singularity/GHCR runs can exercise
ANI-enabled lineage context and native positive IntegronFinder output while
preserving the PanR2 feature contract. It does not revalidate CheckM2,
AMRFinderPlus, QUAST, GTDB-Tk, geNomAD, or MOB-suite because those stages were
intentionally disabled for speed.
