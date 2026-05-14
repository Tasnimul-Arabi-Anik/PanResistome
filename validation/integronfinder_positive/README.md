# IntegronFinder-positive Klebsiella validation subset

This directory contains a compact 10-record `Klebsiella pneumoniae` input
selected from the validated 100-record GHCR/Docker large-mode run. The accessions
were chosen because they produced biological IntegronFinder feature rows in:

`validation_runs/klebsiella_100_ghcr_docker_large/Klebsiella_pneumoniae/panr2_inputs/features/integronfinder.features.tsv`

The subset is intended for fast v0.6.0 positive-call validation of the
PanResistome-native IntegronFinder runner and PanR2 standardized feature export.

## Selection basis

The source run produced:

- 99 analyzed genomes.
- 99/99 QC PASS.
- 418 IntegronFinder feature rows.
- 99 observed IntegronFinder raw tables in the native-runner merge audit.
- Zero unmatched, invalid, or duplicate feature rows in the PanR2 feature
  contract validation.

The selected accessions are among the samples with the highest IntegronFinder
feature counts in the source validation.

## Intended fast validation profile

Use this input when the goal is to exercise positive IntegronFinder output while
keeping unrelated heavy stages disabled:

```bash
nextflow run main.nf \
  --input validation/integronfinder_positive/ncbi_dataset.tsv \
  --outdir validation_runs/integronfinder_positive_fast \
  -profile docker,large \
  --container_image ghcr.io/tasnimul-arabi-anik/panresistome:experimental \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_checkm2 false \
  --run_gtdbtk false \
  --run_quast false \
  --run_ani false \
  --run_mash false \
  --run_amrfinderplus false \
  --panr2_native_feature_runners true \
  --panr2_native_feature_runner_mode parallel \
  --threads 8 \
  --fetchm2_download_workers 2
```

Docker requires either non-root Docker access or an operator-run `sudo docker`
command. Singularity can be used instead when Docker group access is not
available.
