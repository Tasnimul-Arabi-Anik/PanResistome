# Database Automation Matrix

PanResistome should reduce setup friction, but it should not hide large,
licensed, fragile, or unvalidated database downloads. The release rule is:

```text
Automate database setup for default modules only when the path is validated in
PanResistome and has an audit trail.
For opt-in modules, an audited helper can be available when the upstream tool
officially supports download/init, but the module should remain opt-in until
fresh-user validation is complete.
Require an explicit path when the database is restricted, very large, or not
safe to fetch automatically.
```

## Current Behavior

| Module/database | Automated by default? | When enabled? | User path needed? | Notes |
| --- | --- | --- | --- | --- |
| FetchM2 metadata/download | Yes | Default metadata workflow | No | Fetches metadata and assemblies from user input. |
| CheckM2 database | Yes | `--run_checkm2 true` default | No, unless offline/restricted | Downloads under `<outdir>/databases/checkm2`; pass `--checkm2_db` for cached/reproducible runs. |
| ABRicate `ncbi` | Yes | Comprehensive profiles | No | `panr setup-db` plus default forced refresh with `abricate-get_db --force`; audited in `abricate_database_setup_status.tsv`. |
| ABRicate `vfdb` | Yes | Comprehensive profiles | No | Same as above. |
| ABRicate `plasmidfinder` | Yes | Comprehensive profiles | No | Same as above. |
| AMRFinderPlus database | Yes | `--run_amrfinderplus true` | No | `--amrfinderplus_update_db true` by default. |
| IntegronFinder | Yes | Comprehensive profiles | No | Tool/environment driven; no separate user database path in the supported route. |
| MLST schemes | Yes | Comprehensive profiles | No | Uses schemes bundled by the `mlst` package; unsupported organisms produce no-call feature tables. |
| QUAST | Yes | `--run_quast true` | No | Reference-free assembly metrics by default. |
| FastANI/skani | Yes | `--run_ani true` | No | Tool only; large mode guards expensive all-vs-all ANI. |
| Mash | Yes | `--run_mash true` | No | Tool only. |
| MOB-suite | Yes, when explicitly enabled | `--run_mobsuite true` | No, unless offline/restricted | Uses `--mobsuite_db_dir` under `<outdir>/databases/mobsuite`, runs `mob_init` when needed, and attempts ETE `taxa.sqlite` initialization. Still opt-in. |
| Kleborate | Yes | `--run_kleborate true` | No | Small biological validation passed on Klebsiella; still opt-in because organism-specific. |
| MobileElementFinder | Yes, when explicitly enabled | `--panr2_run_mobileelementfinder true` | No | Small 5-genome Klebsiella biological validation passed and produced clean PanR2 feature exports. Still opt-in because broader multi-species/runtime validation is pending. |
| Kaptive | No | `--run_kaptive true` | Yes | Requires `--kaptive_db`; keep explicit until database setup is validated. |
| geNomad | Automated through Docker/cached DB route; Conda first-run still heavy | `--run_genomad true` | No in principle, but cached/container setup is recommended | Uses `--genomad_db_dir` under `<outdir>/databases/genomad` and runs `genomad download-database` when no `--genomad_db` is supplied. Docker validation downloaded geNomad DB v1.9 to a mounted path, and later Docker/GHCR positive-call validations produced 6 region rows across 2 genomes with clean feature-contract validation. On memory-constrained desktops use `--threads 1 --genomad_splits 8 --genomad_sensitivity 3.0`. The earlier Conda route was bottlenecked by first-run environment creation. Use Docker or `-profile genomad_host` with `--genomad_use_host_env true` and a prebuilt/cached geNomad executable/database for easier use. |
| ISfinder-compatible BLAST | No | `--run_isfinder true` | Yes | Requires authorized `--isfinder_db_fasta`; PanResistome must not download or redistribute ISfinder. |
| GTDB-Tk | No | `--run_gtdbtk true` | Yes | Requires large GTDB-Tk reference data; remains off by default. |
| DefenseFinder | No | `--panr2_run_defensefinder true` | Depends on local installation/database | Kept opt-in until dependency/database setup is stable. |
| SerotypeFinder/SCCmecFinder | No runner | Table input only | Yes, for runner workflows outside PanResistome | PanR2-compatible table inputs are supported. |

## What Still Needs Validation?

The helper paths exist, but these should be validated before changing defaults:

1. **MOB-suite cached DB initialization validation** under `<outdir>/databases/mobsuite`, confirming `mob_init` and `taxa.sqlite` creation on a fresh machine.
2. **geNomad scale validation**, expanding the successful two-genome positive-call run to 5-10 genomes and then larger datasets with `--genomad_splits` tuned for available memory.
3. **MobileElementFinder broader validation** on another organism and/or a larger subset before considering default inclusion.
4. **Kaptive database helper** for Klebsiella-focused workflows, but not as a default module because it is organism-specific.

Do not automate by default:

1. **ISfinder**, because users need authorization for the database.
2. **GTDB-Tk**, because the reference database is very large and belongs in an explicit user/HPC setup path.
3. **DefenseFinder runner mode** until its dependency stack and database setup are stable in fresh installs.

## User-Facing Recommendation

For low-hassle comprehensive runs, use the default comprehensive profile and let
PanResistome set up the supported databases. For publication-grade reproducible
reruns, cache the databases and record the manifest files:

```text
panr2_inputs/manifest/database_setup_status.tsv
panr2_inputs/manifest/abricate_database_setup_status.tsv
panr2_inputs/manifest/mobsuite_database_setup_status.tsv
panr2_inputs/manifest/genomad_database_setup_status.tsv
pipeline_versions/
```
