# AGENTS.md - PanResistome Agent Run Rules

## Mission

PanResistome is the execution and reporting layer for bacterial genome feature
analysis. It downloads or accepts genome assemblies, runs QC and annotation
tools, exports strict PanR2-compatible feature tables, and builds user-facing
`basic/`, `important/`, or complete result bundles.

PanResistome should remain practical:

- heavy tool execution, database setup, QC, and provenance belong here;
- PanR2-compatible tables and reports are the downstream contract;
- large raw outputs and databases should not be committed;
- report polish must not weaken complete TSV preservation.

## Start Here

Before changing behavior, read:

```text
README.md
docs/future_agent_runbook.md
docs/database_automation_matrix.md
docs/optional_module_validation_matrix.md
docs/feature_contract_spec.md
docs/troubleshooting.md
```

For report/UI work, also inspect:

```text
scripts/panr2_contract.py
scripts/check_important_report_outputs.py
scripts/check_important_report_visual_layout.py
```

## Default Engineering Rules

- Do not add new external databases by default.
- Do not remove existing `panr2_inputs/` outputs.
- Do not weaken the feature contract.
- Do not use clinical-risk wording for notable-genome prioritization.
- Keep GTDB-Tk, ISfinder, geNomAD, MOB-suite, MobileElementFinder,
  DefenseFinder, and organism-specific typing opt-in unless validation says
  otherwise.
- For large datasets, preserve complete TSVs and cap only report-facing static
  figures or interactive summaries.
- Prefer reusable parameters, scripts, and docs over one-off run commands.

## Shared Database Policy

Large reusable databases should live outside result directories. On shared
workstations, use a stable database root such as:

```text
/mnt/storage/db
```

Use explicit path parameters when available:

```text
--checkm2_db /path/to/uniref100.KO.1.dmnd
--checkm2_db_dir /path/to/checkm2/cache
--gtdbtk_data_path /path/to/gtdbtk/release
--genomad_db /path/to/genomad_db
--genomad_db_dir /path/to/genomad/cache
--mobsuite_db /path/to/mobsuite/db
--mobsuite_db_dir /path/to/mobsuite/cache
--db /path/to/abricate/db
```

When a database is already present, prefer reuse or tool-supported update over
duplicating it into `<outdir>/databases/`. If a module does not yet expose a
stable external database path, document the limitation instead of pretending the
shared path is wired.

## Validation Expectations

For code changes, run the smallest meaningful checks:

```bash
pytest -q
python -m py_compile scripts/*.py
bash -n scripts/bootstrap.sh
nextflow config -profile docker,large -flat
git diff --check
```

For important-report changes, validate a real or existing output bundle:

```bash
python scripts/check_important_report_outputs.py <sample-dir>
python scripts/check_important_report_visual_layout.py <sample-dir> --browser skip
```

Use browser screenshots for visual changes when a local or remote report bundle
is available.

## Done Means

A task is not done until:

- generated commands point to stable paths;
- important warnings/cautions are documented;
- complete TSVs remain available;
- tests or targeted validation pass, or failures are explicitly reported;
- future agents can reproduce the next command without guessing.
