# PanResistome Validation Matrix

This matrix summarizes the major validation evidence used for release hardening. Generated run directories, work directories, downloaded FASTA files, and downloaded databases are intentionally not committed.

| Validation | Organism | Input records | Genomes downloaded/analyzed | QC PASS | GTDB-Tk | CheckM2 | QUAST | ANI | Mash | AMRFinderPlus | Native runners | Runner mode | Threads | CheckM2 threads | Feature rows | Unmatched | Invalid | Duplicate | Dashboard | Status | Report |
| --- | --- | ---: | ---: | ---: | --- | --- | --- | --- | --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- |
| Delftia fresh-clone v0.2.2 | `Delftia tsuruhatensis` | 45 | 45 | 45 | off | auto DB | on | on | on | on | no | PanR2 internal | 4 | 4 default cap | 291 | 0 | 0 | 0 | yes | PASS | `validation/delftia_tsuruhatensis_current/FRESH_CLONE_VALIDATION_RESULTS.md` |
| Delftia native-runner serial | `Delftia tsuruhatensis` | 45 | 45 | 45 | off | cached DB | on | on | on | on | yes | serial | 4 | 4 default cap | 201 | 0 | 0 | 0 | yes | PASS | `validation/delftia_tsuruhatensis_current/NATIVE_FEATURE_RUNNER_VALIDATION_RESULTS.md` |
| Delftia native-runner parallel | `Delftia tsuruhatensis` | 45 local sample | 45 | 45 | off | not rerun | not rerun | not rerun | not rerun | cached output | yes | parallel | 16 | not run | 201 | 0 | 0 | 0 | yes | PASS | `validation/delftia_tsuruhatensis_current/NATIVE_PARALLEL_COMPARISON.md` |
| Klebsiella pneumoniae 100 parallel | `Klebsiella pneumoniae` | 100 | 99 | 99 | off | cached auto DB | on | on | on | on | yes | parallel | 16 | 2 | 12838 | 0 | 0 | 0 | yes | PASS | `validation/klebsiella_pneumoniae_100/VALIDATION_RESULTS.md` |
| Klebsiella pneumoniae 300 large mode | `Klebsiella pneumoniae` | 300 | 299 | 299 | off | off | on | off | on | off | yes | parallel | 16 | not run | 36638 | 0 | 0 | 0 | yes | PASS | `validation/klebsiella_pneumoniae_300/LARGE_MODE_CHECKM2_OFF_VALIDATION_RESULTS.md` |
| Klebsiella pneumoniae 5 geNomAD Docker/GHCR | `Klebsiella pneumoniae` | 5 | 5 | 5 | off | off | off | off | off | off | no | geNomAD jobs=2 | 2 | not run | 601 | 0 | 0 | 0 | yes | PASS | `validation/optional_runner_biological/GENOMAD_POSITIVE_CALL_VALIDATION_RESULTS.md` |
| Klebsiella pneumoniae 10 v0.5.0 Docker/GHCR | `Klebsiella pneumoniae` | 10 | 10 | 10 | off | off | off | off | on | on | yes | parallel | 8 | not run | 1348 | 0 | 0 | 0 | yes | PASS | `validation/klebsiella_pneumoniae_10/V0_5_0_DOCKER_VALIDATION_RESULTS.md` |
| Klebsiella pneumoniae IntegronFinder-positive GHCR evidence | `Klebsiella pneumoniae` | 100 | 99 | 99 | off | off | off | off | on | off | yes | parallel | 8 | not run | 11488 | 0 | 0 | 0 | yes | PASS | `validation/integronfinder_positive/VALIDATION_RESULTS.md` |
| Klebsiella pneumoniae 10 ANI-lineage Singularity/GHCR | `Klebsiella pneumoniae` | 10 | 10 | 10 | off | off | off | on | on | off | yes | parallel | 8 | not run | 1248 | 0 | 0 | 0 | yes | PASS | `validation/klebsiella_pneumoniae_10/ANI_LINEAGE_SINGULARITY_VALIDATION_RESULTS.md` |
| Acinetobacter pittii 5 Docker/GHCR comprehensive | `Acinetobacter pittii` | 5 | 5 | 5 | off | auto DB | on | on | on | on | yes | parallel | 8 | 4 | 630 | 0 | 0 | 0 | yes | PASS | `validation/deployment/ACINETOBACTER_PITTII_5_DOCKER_STABILITY.md` |

## Current Interpretation

- The fresh-clone Delftia run historically proved the public command could install environments, fetch/build legal databases, and complete without a user-provided CheckM2 database path. A 2026-05-16 fixture reproduced a stale CheckM2 TensorFlow/Keras package failure; the updated CheckM2 build-1/CPU-TensorFlow route now loads the model and writes real prediction rows in both Conda and Docker/GHCR fixture smokes. The 5-genome Acinetobacter Docker/GHCR comprehensive run confirms the current container route works with CheckM2 auto-download and no manual database path.
- The native-runner Delftia runs prove ABRicate, IntegronFinder, and MLST can be owned by PanResistome and passed to PanR2 as precomputed results.
- The Klebsiella run proves the parallel native-runner path handles a biologically richer 100-record input with clean PanR2 feature-contract validation.
- The 300-record large-mode run proves report safeguards and native feature runners scale on a desktop-safe profile when heavyweight optional stages are intentionally disabled.
- The 5-genome Docker/GHCR geNomAD run proves bounded optional geNomAD execution can produce positive biological prophage/plasmid-like region calls with clean PanR2 handoff output.
- The 10-genome Docker/GHCR run proves the v0.5.0 reproducibility manifest,
  feature-contract manifest, lineage-aware summaries, pan-feature diversity
  summaries, statistical summary, and new HTML pages work on real
  feature-rich Klebsiella data with clean PanR2 feature-contract validation.
- The IntegronFinder-positive evidence confirms that PanResistome-native
  IntegronFinder can produce real biological feature rows at scale, merge 99/99
  raw tables, and pass PanR2 feature-contract validation with zero unmatched,
  invalid, or duplicate rows.
- The 10-genome ANI-lineage Singularity/GHCR run proves the public image can
  exercise skani all-vs-all ANI, populate ANI-cluster lineage context, preserve
  positive native IntegronFinder handoff rows in strict PanR2 features, and pass
  feature-contract validation after the v0.6.0 native-handoff export fix.
- The 5-genome Acinetobacter Docker/GHCR run proves the current recommended
  CheckM2-on comprehensive route can auto-download CheckM2 and geNomAD
  databases, run QUAST/ANI/Mash/AMRFinderPlus/native ABRicate/IntegronFinder/MLST,
  and pass PanR2 schema validation with zero unmatched, invalid, or duplicate
  feature rows.
- GTDB-Tk remains intentionally excluded from these validations because its reference database is large and should remain opt-in.
