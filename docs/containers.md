# Container execution notes

PanResistome currently ships validated Conda/Mamba workflows. Docker and Apptainer/Singularity execution are v0.4.0 deployment targets and should be treated as experimental until a container smoke test and real-data validation are documented.

## Design goal

The container strategy should preserve the current architecture:

```text
PanResistome = tool execution, QC, database setup, feature export
PanR2        = standardized feature analysis and reporting
```

Each heavy tool may need its own image or carefully separated environment. A single monolithic image can be convenient, but it can also become difficult to rebuild and debug.

## Database paths

Containerized runs must mount large external databases into the container:

- CheckM2 database
- GTDB-Tk data, if enabled
- geNomad database, if enabled
- Kaptive database, if enabled
- authorized ISfinder FASTA/BLAST database, if enabled

Keep these databases outside `work/` and outside generated result folders.

## Practical first container targets

1. Keep the existing `-profile test` workflow working without containers.
2. Add a small container smoke test for local fixtures.
3. Validate the standard comprehensive command without GTDB-Tk.
4. Document database mounts for CheckM2 and optional heavy databases.
5. Only then advertise Docker/Apptainer profiles as supported.

## Example future pattern

This is a planning example, not a validated command:

```bash
nextflow run main.nf \
  --input validation/delftia_tsuruhatensis_current/ncbi_dataset.tsv \
  --outdir validation_runs/delftia_container \
  -profile apptainer,large \
  --analysis_profile comprehensive \
  --qc_filter true \
  --run_gtdbtk false \
  --run_amrfinderplus true \
  --checkm2_db /databases/checkm2/uniref100.KO.1.dmnd
```

Do not assume host database paths are visible inside a container. Mounts must expose the same paths that the Nextflow command passes to PanResistome.

