# PanR2 Feature Contract Specification

Current contract version: `1.0`

This specification defines the standardized feature tables that PanResistome
exports for PanR2 and downstream comparative analysis. Raw tool output remains
available for traceability, but `panr2_inputs/features/*.features.tsv` is the
strict downstream layer.

## Required columns

Every feature-like module table must contain these columns:

```text
sample_id
assembly_accession
database
feature_id
feature_category
presence
identity
coverage
contig
start
end
tool
tool_version
database_version
```

`presence` must be `1` for present feature calls or `0` for explicit absence
records. Current PanResistome feature exports are presence-oriented and normally
write `1`.

## Optional standardized columns

Downstream tools should preserve these columns when present:

```text
feature_name
feature_description
feature_subcategory
mechanism
drug_class
product
sequence_id
strand
source_table
source_file
source_database
raw_feature_id
raw_category
raw_method
evidence_type
confidence
notes
```

## Known database names

Current known values include:

```text
amr
amrfinderplus
vfdb
plasmidfinder
isfinder
mobileelementfinder
integronfinder
mlst
mobsuite
defensefinder
prophage
genomad
iceberg
kleborate
kaptive
ectyper
serotypefinder
sccmecfinder
ani
assembly_qc
quast
mash
custom
```

New module names may be added in minor releases when the module exports raw
output, status/audit rows, a standardized feature table if feature-like, and an
entry in `feature_completeness_audit.tsv`.

## Evidence and confidence values

Recommended values:

```text
evidence_type: sequence_match, tool_call, typing_call, assembly_metric,
               cooccurrence, proximity, unknown
confidence:    high, medium, low, unknown
```

Cross-database proximity outputs use:

```text
evidence_level: same_genome, same_contig, within_10kb, overlapping, adjacent,
                unknown
```

These evidence levels describe context strength only. They do not prove gene
transfer, phenotype, expression, or plasmid localization.

## Machine-readable manifest

Each PanR2 handoff export writes:

```text
panr2_inputs/manifest/feature_contract.json
```

This JSON records `contract_version`, `schema_version`, required columns,
optional columns, known database names, recommended controlled values, and the
backward-compatibility statement.

## Compatibility policy

Feature tables produced by v0.3.x and v0.4.0 remain valid under contract `1.0`
when they contain all required columns. Optional columns may be absent. Parsers
must preserve optional columns they understand and ignore unknown extra columns
unless a specific analysis requires them.

Breaking changes to required columns or value semantics require a new minor
release note and an updated `contract_version`.
