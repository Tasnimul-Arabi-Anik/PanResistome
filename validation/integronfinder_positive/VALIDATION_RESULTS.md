# IntegronFinder-positive biological validation evidence

Date documented: 2026-05-14

This validation note records real biological IntegronFinder-positive evidence
from the previously completed 100-record `Klebsiella pneumoniae` GHCR/Docker
large-mode run. It is used as v0.6.0 planning evidence because the compact
10-record subset in this directory was selected directly from these positive
samples.

## Source run

- Output directory: `validation_runs/klebsiella_100_ghcr_docker_large`
- Organism: `Klebsiella pneumoniae`
- Input records: 100
- Genomes analyzed: 99
- QC PASS: 99
- Container path: GHCR/Docker
- Large mode: enabled
- CheckM2: disabled
- GTDB-Tk: disabled
- ANI: disabled
- AMRFinderPlus: disabled
- PanResistome-native runners: enabled
- Native runner mode: parallel

## IntegronFinder evidence

From
`Klebsiella_pneumoniae/panr2_inputs/manifest/native_runner_merge_audit.tsv`:

| Module | Expected raw tables | Observed raw tables | Samples processed | Samples failed | Feature rows | Unique features | Status |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| integronfinder | 99 | 99 | 99 | 0 | 418 | 3 | PASS |

The standardized feature table
`Klebsiella_pneumoniae/panr2_inputs/features/integronfinder.features.tsv`
contains 418 biological feature rows plus a header.

## PanR2 feature-contract validation

From
`Klebsiella_pneumoniae/panr2_inputs/manifest/schema_validation_summary.txt`:

```text
feature_files_checked=5
feature_rows=11488
databases_seen=amr,integronfinder,mlst,plasmidfinder,vfdb
unmatched_feature_rows=0
invalid_feature_rows=0
duplicate_feature_rows=0
```

The `samples_seen` value in this report is larger than the metadata accession
count because PanR2 validates feature sample identifiers across both assembly
accessions and normalized sample IDs. The important release gate is that
unmatched feature rows remained zero.

## Interpretation

This evidence confirms that the PanResistome-native IntegronFinder path can
produce real biological positive calls, merge raw per-sample outputs, export
standardized PanR2 feature rows, and pass strict feature-contract validation.

It does not revalidate CheckM2, GTDB-Tk, ANI, or AMRFinderPlus because those
stages were intentionally disabled for this large-mode run.
