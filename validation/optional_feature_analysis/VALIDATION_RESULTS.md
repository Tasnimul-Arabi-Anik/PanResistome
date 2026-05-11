# Optional Feature Analysis Validation

Date: 2026-05-11

Purpose: verify that opt-in module outputs can be converted into PanR2-compatible feature tables and analyzed through the same standardized layer used for AMR and VFDB.

This validation uses small local fixture tables, not large external databases. It validates the PanR2 handoff/analysis behavior for optional outputs and complements the documented biological Kleborate/MOB-suite runner validation.

## Result

Status: PASS

Feature rows by database are in `feature_counts.tsv`.

Schema validation summary:

```text
feature_files_checked=13
feature_rows=23
databases_seen=amr,defensefinder,ectyper,isfinder,kaptive,kleborate,mobileelementfinder,mobsuite,plasmidfinder,prophage,sccmecfinder,serotypefinder,vfdb
samples_seen=6
metadata_accessions=6
unmatched_feature_rows=0
invalid_feature_rows=0
duplicate_feature_rows=0
```

Generated feature tables:

- `amr.features.tsv`
- `defensefinder.features.tsv`
- `ectyper.features.tsv`
- `isfinder.features.tsv`
- `kaptive.features.tsv`
- `kleborate.features.tsv`
- `mobileelementfinder.features.tsv`
- `mobsuite.features.tsv`
- `plasmidfinder.features.tsv`
- `prophage.features.tsv`
- `sccmecfinder.features.tsv`
- `serotypefinder.features.tsv`
- `vfdb.features.tsv`

Generated analysis/report checks are in `analysis_outputs.tsv`.

## Interpretation

This confirms that MobileElementFinder, ISfinder-style BLAST, MOB-suite, geNomad/prophage, DefenseFinder, Kleborate, Kaptive, ECTyper, SerotypeFinder, and SCCmecFinder tables can produce standardized feature rows, `all_features.tsv`, feature matrices, co-occurrence/proximity outputs, metadata usability outputs, top-feature summaries, and HTML report pages.

This does not claim that every external runner/database can be installed and run biologically on a fresh desktop. Runner-mode status remains separated in `docs/optional_module_validation_matrix.md`.
