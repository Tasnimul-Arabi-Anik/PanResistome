# PanResistome Release Reliability Checklist

This checklist turns the public-user reliability goals into release gates. A PanResistome release should not be marked stable until the standard comprehensive path passes these checks from a fresh clone.

## Release Gates

1. Fresh clone starts cleanly
   - Evidence: `git clone`, `cd PanResistome`, and `nextflow run main.nf -profile test` complete without local-only files.

2. One canonical comprehensive command exists
   - Evidence: the command in `docs/remote_user_validation.md` runs with GTDB-Tk disabled and without user-supplied CheckM2, AMRFinderPlus, or ABRicate setup steps.

3. FetchM2 is the default metadata engine
   - Evidence: each sample directory contains FetchM2 outputs plus legacy-compatible `metadata_output/ncbi_clean.csv`.

4. Assembly download is automatic
   - Evidence: FetchM2 native downloader writes downloaded FASTA files and failed-download reports under the run directory.

5. Sequence QC always runs before annotation
   - Evidence: `qc/sequence_qc.csv` exists and downstream filtered FASTA files are derived from the combined QC decision.

6. CheckM2 works without a pre-existing local database
   - Evidence: when `--checkm2_db` is omitted, the run downloads/caches the database under `<outdir>/databases/checkm2/` and records the path in logs/manifests.

7. QUAST, ANI/skani, and Mash are optional but validated in the standard comprehensive command
   - Evidence: their output tables are present when `--run_quast true --run_ani true --run_mash true` are used.

8. Combined QC decision is explicit
   - Evidence: `qc/qc_master_report.csv`, pass/fail/warning sample lists, and `qc/excluded_for_panr2.csv` exist.

9. QC filtering cannot silently leave zero samples
   - Evidence: `--qc_filter true` fails at combined QC with a clear message if no FASTA files remain.

10. Required database/tool setup is audited
    - Evidence: `panr2_inputs/manifest/database_setup_status.tsv` exists and required rows have no `FAIL` status.

11. ABRicate NCBI/VFDB/PlasmidFinder are available after setup
    - Evidence: comprehensive mode confirms `ncbi`, `vfdb`, and `plasmidfinder` through `abricate --list` after `panr setup-db`.

12. AMRFinderPlus setup is automatic when enabled
    - Evidence: `--run_amrfinderplus true` records the AMRFinderPlus executable version, database version, setup action, and per-sample status.

13. Fragile or restricted modules are honest opt-ins
    - Evidence: GTDB-Tk, MobileElementFinder, ISfinder FASTA, MOB-suite, geNomad, DefenseFinder, and organism-specific typing are documented as optional or table-input modules unless their required databases are supplied.

14. Every enabled feature-like module exports the PanR2 contract
    - Evidence: `panr2_inputs/features/<database>.features.tsv` exists for enabled modules that create feature rows.

15. A merged feature table is produced
    - Evidence: `panr2_inputs/features/all_features.tsv` exists and includes all successful feature-family exports.

16. Feature schema validation is strict
    - Evidence: `schema_validation_summary.txt`, `unmatched_features.csv`, `duplicate_features.csv`, and `invalid_feature_rows.csv` exist; release validation should have zero unmatched, duplicate, and invalid rows unless documented.

17. Module completeness is audited
    - Evidence: `feature_completeness_audit.tsv` and `module_status_summary.tsv` list every expected module as pass, warning, failure, or skipped with a reason.

18. PanR2 produces a combined report
    - Evidence: `<sample>/report/index.html` exists and links QC, metadata, database-specific outputs, cross-database results, citations, and software versions.

19. Validation evidence is summarized
    - Evidence: `scripts/summarize_validation_run.py --run-dir <run> --out-dir <run>` writes `validation_summary.csv` and `validation_summary.md`.

20. Release documentation separates stable, optional, and experimental paths
    - Evidence: README, changelog, and validation documents distinguish the stable public command from optional heavy modules and table-input-only modules.

## Current Standard Validation Command

```bash
nextflow run main.nf \
  --input validation/delftia_tsuruhatensis_current/ncbi_dataset.tsv \
  --outdir validation_runs/delftia_fresh \
  -profile conda,mamba \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_quast true \
  --run_ani true \
  --run_mash true \
  --run_amrfinderplus true \
  --threads 4 \
  --fetchm2_download_workers 2
```

## Known Release Caveats

- GTDB-Tk is intentionally excluded from the public comprehensive validation because its reference database is large.
- ISfinder is not auto-downloaded or redistributed. Use `--run_isfinder true --isfinder_db_fasta <authorized.fasta>` when the user has an authorized local database.
- MobileElementFinder remains opt-in because upstream output parsing has failed on otherwise valid assemblies during real validation.
- Standard ABRicate, IntegronFinder, MLST, and opt-in MobileElementFinder execution is now owned by the PanResistome native feature-runner stage when `--panr2_native_feature_runners true`. Large fragmented assemblies can still make IntegronFinder slow inside that stage; future optimization should split per-assembly execution into finer-grained Nextflow channels.
