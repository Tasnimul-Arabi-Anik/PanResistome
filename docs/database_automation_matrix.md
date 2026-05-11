# Database Automation Matrix

PanResistome should reduce setup friction, but it should not hide large,
licensed, fragile, or unvalidated database downloads. The release rule is:

```text
Automate database setup when it is legally redistributable or tool-supported,
validated in PanResistome, and has an audit trail.
Require an explicit path when the database is restricted, very large, or not
yet validated as a fresh-user download.
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
| MOB-suite | Partly | `--run_mobsuite true` | Recommended for reliable/offline runs | Runner can use MOB-suite runtime initialization, but biological validation used a preinitialized DB with `taxa.sqlite`. Keep opt-in. |
| Kleborate | Yes | `--run_kleborate true` | No | Small biological validation passed on Klebsiella; still opt-in because organism-specific. |
| Kaptive | No | `--run_kaptive true` | Yes | Requires `--kaptive_db`; keep explicit until database setup is validated. |
| geNomad | No | `--run_genomad true` | Yes | Requires `--genomad_db`; database is large, so keep explicit until fresh-user download is validated. |
| ISfinder-compatible BLAST | No | `--run_isfinder true` | Yes | Requires authorized `--isfinder_db_fasta`; PanResistome must not download or redistribute ISfinder. |
| GTDB-Tk | No | `--run_gtdbtk true` | Yes | Requires large GTDB-Tk reference data; remains off by default. |
| DefenseFinder | No | `--panr2_run_defensefinder true` | Depends on local installation/database | Kept opt-in until dependency/database setup is stable. |
| SerotypeFinder/SCCmecFinder | No runner | Table input only | Yes, for runner workflows outside PanResistome | PanR2-compatible table inputs are supported. |

## Should More Modules Be Automated?

Yes, but only module-by-module after validation.

Good candidates for future automation:

1. **MOB-suite cached DB initialization** under `<outdir>/databases/mobsuite`, because the runner is now biologically validated on a small Klebsiella set. This still needs a fresh-machine test that confirms `taxa.sqlite` is created reliably.
2. **Kaptive database helper** for Klebsiella-focused workflows, but not as a default module because it is organism-specific.
3. **geNomad database helper** only after testing disk size, download stability, and resume behavior.

Do not automate by default:

1. **ISfinder**, because users need authorization for the database.
2. **GTDB-Tk**, because the reference database is very large and belongs in an explicit user/HPC setup path.
3. **Fragile optional runners** such as MobileElementFinder or DefenseFinder until their dependency stacks are stable in fresh installs.

## User-Facing Recommendation

For low-hassle comprehensive runs, use the default comprehensive profile and let
PanResistome set up the supported databases. For publication-grade reproducible
reruns, cache the databases and record the manifest files:

```text
panr2_inputs/manifest/database_setup_status.tsv
panr2_inputs/manifest/abricate_database_setup_status.tsv
pipeline_versions/
```
