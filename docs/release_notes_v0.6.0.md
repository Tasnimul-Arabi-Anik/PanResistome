# PanResistome v0.6.0 Release Notes

## Highlights

- ANI-enabled Singularity/GHCR validation on a compact 10-genome
  `Klebsiella pneumoniae` dataset.
- Positive IntegronFinder biological validation evidence.
- Native handoff export fix so IntegronFinder and MobileElementFinder result
  rows under `tool_results/*/panr2_inputs/` are preserved in strict PanR2
  feature tables.
- Documentation for slow first-time Singularity GHCR-to-SIF conversion using a
  longer `singularity.pullTimeout`.

## Key Validation Result

The v0.6.0 ANI-lineage validation completed 14/14 Nextflow processes through
Singularity/GHCR with:

- 10 input records.
- 10 genomes analyzed.
- skani all-vs-all ANI enabled.
- Mash enabled.
- PanResistome-native ABRicate, IntegronFinder, and MLST runners enabled.
- 1,248 standardized feature rows.
- 44 standardized IntegronFinder rows.
- 0 unmatched feature rows.
- 0 invalid feature rows.
- 0 duplicate feature rows.
- ANI-cluster lineage context and PanR2 HTML handoff pages generated.

## Recommended Announcement Wording

PanResistome v0.6.0 is a deployment and validation release. It adds
ANI-enabled Singularity/GHCR validation, documents positive IntegronFinder
biological output, and fixes native handoff export so positive IntegronFinder
and MobileElementFinder rows cannot disappear from PanR2 feature tables.

## Scope

This release does not add new default databases. GTDB-Tk and ISfinder remain
explicit path-based opt-in modules, and heavyweight optional modules remain
opt-in.
