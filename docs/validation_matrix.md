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

## Current Interpretation

- The fresh-clone Delftia run proves the public command can install environments, fetch/build legal databases, and complete without a user-provided CheckM2 database path.
- The native-runner Delftia runs prove ABRicate, IntegronFinder, and MLST can be owned by PanResistome and passed to PanR2 as precomputed results.
- The Klebsiella run proves the parallel native-runner path handles a biologically richer 100-record input with clean PanR2 feature-contract validation.
- The 300-record large-mode run proves report safeguards and native feature runners scale on a desktop-safe profile when heavyweight optional stages are intentionally disabled.
- The 5-genome Docker/GHCR geNomAD run proves bounded optional geNomAD execution can produce positive biological prophage/plasmid-like region calls with clean PanR2 handoff output.
- GTDB-Tk remains intentionally excluded from these validations because its reference database is large and should remain opt-in.
