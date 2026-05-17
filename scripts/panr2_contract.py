#!/usr/bin/env python3
"""Create and validate PanR2 contract feature tables from PanResistome outputs."""

from __future__ import annotations

import csv
import html
import itertools
import json
import math
import re
import struct
import zipfile
import zlib
from collections import Counter, defaultdict
from pathlib import Path


CONTRACT_VERSION = "1.0"
SCHEMA_VERSION = "1.0"

CONTRACT_COLUMNS = [
    "sample_id",
    "assembly_accession",
    "database",
    "feature_id",
    "feature_category",
    "presence",
    "identity",
    "coverage",
    "contig",
    "start",
    "end",
    "tool",
    "tool_version",
    "database_version",
]

OPTIONAL_CONTRACT_COLUMNS = [
    "feature_name",
    "feature_description",
    "feature_subcategory",
    "mechanism",
    "drug_class",
    "product",
    "sequence_id",
    "strand",
    "source_table",
    "source_file",
    "source_database",
    "raw_feature_id",
    "raw_category",
    "raw_method",
    "evidence_type",
    "confidence",
    "notes",
]

FEATURE_COLUMNS = CONTRACT_COLUMNS + OPTIONAL_CONTRACT_COLUMNS

FEATURE_CONTRACT_ALLOWED_VALUES = {
    "presence": ["0", "1"],
    "database": sorted([
        "amr",
        "amrfinderplus",
        "vfdb",
        "plasmidfinder",
        "isfinder",
        "mobileelementfinder",
        "integronfinder",
        "mlst",
        "mobsuite",
        "defensefinder",
        "prophage",
        "genomad",
        "iceberg",
        "kleborate",
        "kaptive",
        "ectyper",
        "serotypefinder",
        "sccmecfinder",
        "ani",
        "assembly_qc",
        "quast",
        "mash",
        "custom",
    ]),
    "evidence_type": [
        "sequence_match",
        "tool_call",
        "typing_call",
        "sequence_type_call",
        "allele_call",
        "kleborate_call",
        "kleborate_amr_marker",
        "assembly_metric",
        "cooccurrence",
        "proximity",
        "unknown",
    ],
    "confidence": ["high", "medium", "low", "unknown"],
    "evidence_level": ["same_genome", "same_contig", "within_10kb", "overlapping", "adjacent", "unknown"],
    "module_status": [
        "PASS",
        "WARNING_EMPTY",
        "SKIPPED_NOT_ENABLED",
        "SKIPPED_INAPPLICABLE",
        "FAIL",
        "FAIL_MISSING_FEATURE_TABLE",
        "ERROR",
    ],
}

KNOWN_DATABASES = {
    "amr",
    "amrfinderplus",
    "vfdb",
    "plasmidfinder",
    "isfinder",
    "mobileelementfinder",
    "integronfinder",
    "mlst",
    "mobsuite",
    "defensefinder",
    "prophage",
    "genomad",
    "iceberg",
    "kleborate",
    "kaptive",
    "ectyper",
    "serotypefinder",
    "sccmecfinder",
    "ani",
    "assembly_qc",
    "quast",
    "mash",
    "custom",
}

METADATA_ALIASES = {
    "country": ["Country", "Geographic Location", "country", "Country_SD", "geo_loc_name_country"],
    "continent": ["Continent", "continent"],
    "subcontinent": ["Subcontinent", "subcontinent"],
    "host": ["Host", "Host_SD", "host", "host_scientific_name"],
    "host_group": ["Host_Group", "Host_Group_SD", "Host_Rank", "host_group"],
    "isolation_source": ["Isolation Source", "Isolation_Source", "Isolation_Source_SD", "isolation_source"],
    "sample_type": ["Sample_Type", "Sample_Type_SD", "sample_type"],
    "environment_medium": ["Environment_Medium", "Environment_Medium_SD", "environmental_medium"],
    "collection_year": ["Collection_Year", "Collection Date", "year", "collection_year"],
    "organism": ["Organism Name", "Organism", "organism"],
    "species": ["Species", "species"],
    "genus": ["Genus", "genus"],
    "bioproject": ["Assembly BioProject Accession", "BioProject", "BioProject Accession"],
    "biosample": ["Assembly BioSample Accession", "BioSample", "BioSample Accession"],
    "assembly_level": ["Assembly Level", "assembly_level"],
    "submitter": ["Assembly Submitter", "Submitter", "submitter"],
}

ACCESSION_RE = re.compile(r"(GC[AF]_\d+\.\d+)")
FASTA_SUFFIXES = [
    ".fasta.gz",
    ".fna.gz",
    ".fa.gz",
    ".fasta",
    ".fna",
    ".fa",
    ".ffn",
    ".gz",
]


ABRICATE_DATABASE_MAP = {
    "ncbi": "amr",
    "amr": "amr",
    "vfdb": "vfdb",
    "plasmidfinder": "plasmidfinder",
    "mobileelementfinder": "mobileelementfinder",
    "isfinder": "isfinder",
    "integronfinder": "integronfinder",
    "iceberg": "iceberg",
    "defensefinder": "defensefinder",
    "prophage": "prophage",
    "mobsuite": "mobsuite",
    "kleborate": "kleborate",
    "kaptive": "kaptive",
    "ectyper": "ectyper",
    "serotypefinder": "serotypefinder",
    "sccmecfinder": "sccmecfinder",
}

ABRICATE_NATIVE_DATABASES = {"amr", "vfdb", "plasmidfinder"}

OPTIONAL_TABLE_DATABASES = {
    "mobileelementfinder",
    "isfinder",
    "mobsuite",
    "defensefinder",
    "prophage",
    "kleborate",
    "kaptive",
    "ectyper",
    "serotypefinder",
    "sccmecfinder",
}


AMRFINDER_FEATURE_COLUMNS = [
    "Gene symbol",
    "Gene",
    "gene",
    "Element symbol",
    "Element name",
    "Sequence name",
    "Name of closest sequence",
    "HMM id",
]

AMRFINDER_CATEGORY_COLUMNS = [
    "Class",
    "Subclass",
    "Element subtype",
    "Element type",
    "Scope",
    "Method",
]


def detect_delimiter(path: Path) -> str:
    text = path.read_text(errors="ignore")
    first = next((line for line in text.splitlines() if line.strip() and not line.startswith("#")), "")
    return "\t" if first.count("\t") >= first.count(",") else ","


def read_table(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    delimiter = detect_delimiter(path)
    with path.open(newline="", errors="ignore") as handle:
        lines = [
            line
            for line in handle
            if line.strip() and (not line.lstrip().startswith("#") or line.lstrip().startswith("#FILE"))
        ]
    if not lines:
        return []
    reader = csv.DictReader(lines, delimiter=delimiter)
    if not reader.fieldnames:
        return []
    return [
        {str(key).strip(): str(value or "").strip() for key, value in row.items() if key is not None}
        for row in reader
    ]


def read_header(path: Path) -> list[str]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    delimiter = detect_delimiter(path)
    with path.open(newline="", errors="ignore") as handle:
        for line in handle:
            if line.strip() and (not line.lstrip().startswith("#") or line.lstrip().startswith("#FILE")):
                return [part.strip() for part in line.rstrip("\n").split(delimiter)]
    return []


def write_rows(path: Path, rows: list[dict[str, str]], fieldnames: list[str], delimiter: str = "\t") -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter=delimiter, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return str(path)


def first_value(row: dict[str, str], candidates: list[str], default: str = "") -> str:
    lowered = {key.lower(): value for key, value in row.items()}
    for name in candidates:
        if name in row and str(row[name]).strip():
            return str(row[name]).strip()
        value = lowered.get(name.lower())
        if value and str(value).strip():
            return str(value).strip()
    return default


def is_missing_value(value: str) -> bool:
    text = str(value or "").strip()
    return not text or text in {"-", ".", "?"} or text.lower() in {"na", "n/a", "nan", "none", "null", "unknown"}


def is_placeholder_mlst_feature(value: str) -> bool:
    text = str(value or "").strip()
    if is_missing_value(text):
        return True
    return re.fullmatch(r"[-_:\s]*ST[-_:\s]*", text, flags=re.IGNORECASE) is not None


def clean_sample_id(value: str) -> str:
    text = Path(str(value or "").strip()).name
    for suffix in FASTA_SUFFIXES:
        if text.lower().endswith(suffix):
            text = text[: -len(suffix)]
            break
    return text.strip()


def extract_accession(value: str) -> str:
    match = ACCESSION_RE.search(str(value or ""))
    return match.group(1) if match else ""


def load_sample_map(sample_dir: Path) -> dict[str, str]:
    candidates = [
        sample_dir / "metadata_output" / "sample_map.csv",
        sample_dir / "panr2_inputs" / "metadata" / "sample_map.csv",
    ]
    mapping: dict[str, str] = {}
    for path in candidates:
        for row in read_table(path):
            accession = first_value(row, ["Assembly Accession", "assembly_accession", "assembly"], "")
            if not accession:
                continue
            for key in ["sample_id", "sample", "sequence_file", "file", "#FILE"]:
                sample = first_value(row, [key], "")
                if sample:
                    mapping[clean_sample_id(sample)] = accession
                    mapping[sample] = accession
            mapping[accession] = accession
    return mapping


def load_metadata_accessions(sample_dir: Path) -> set[str]:
    accessions: set[str] = set()
    for path in [
        sample_dir / "metadata_output" / "ncbi_clean_qc_pass.csv",
        sample_dir / "metadata_output" / "ncbi_enriched.csv",
        sample_dir / "metadata_output" / "ncbi_clean.csv",
        sample_dir / "metadata_output" / "fetchm2_clean.csv",
        sample_dir / "panr2_inputs" / "metadata" / "ncbi_clean_qc_pass.csv",
        sample_dir / "panr2_inputs" / "metadata" / "ncbi_enriched.csv",
        sample_dir / "panr2_inputs" / "metadata" / "ncbi_clean.csv",
        sample_dir / "panr2_inputs" / "metadata" / "fetchm2_clean.csv",
    ]:
        for row in read_table(path):
            accession = first_value(row, ["Assembly Accession", "assembly_accession", "Assembly"], "")
            if accession:
                accessions.add(accession)
    return accessions


def resolve_assembly(sample: str, sample_map: dict[str, str]) -> str:
    sample_id = clean_sample_id(sample)
    return sample_map.get(sample) or sample_map.get(sample_id) or extract_accession(sample) or sample_id


def safe_number(value: str) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    try:
        return str(float(text))
    except ValueError:
        return text


def contract_row(
    sample: str,
    database: str,
    feature_id: str,
    feature_category: str = "",
    identity: str = "",
    coverage: str = "",
    contig: str = "",
    start: str = "",
    end: str = "",
    tool: str = "",
    tool_version: str = "",
    database_version: str = "",
    sample_map: dict[str, str] | None = None,
    **optional_values: str,
) -> dict[str, str]:
    sample_map = sample_map or {}
    sample_id = clean_sample_id(sample)
    row = {
        "sample_id": sample_id,
        "assembly_accession": resolve_assembly(sample, sample_map),
        "database": database,
        "feature_id": str(feature_id or "").strip(),
        "feature_category": str(feature_category or "").strip(),
        "presence": "1",
        "identity": safe_number(identity),
        "coverage": safe_number(coverage),
        "contig": str(contig or "").strip(),
        "start": str(start or "").strip(),
        "end": str(end or "").strip(),
        "tool": tool,
        "tool_version": tool_version,
        "database_version": database_version,
    }
    for column in OPTIONAL_CONTRACT_COLUMNS:
        row[column] = str(optional_values.get(column, "") or "").strip()
    return row


def database_from_abricate_path(path: Path, sample_dir: Path) -> str:
    parts = [part.lower() for part in path.relative_to(sample_dir).parts]
    for part in parts:
        if part in ABRICATE_DATABASE_MAP:
            return ABRICATE_DATABASE_MAP[part]
    stem = path.stem.lower()
    for key, value in ABRICATE_DATABASE_MAP.items():
        if key in stem:
            return value
    return "feature"


def tool_for_feature_table(database: str, row: dict[str, str]) -> str:
    explicit_tool = first_value(row, ["tool", "Tool", "TOOL"], "")
    if explicit_tool:
        return explicit_tool
    return "abricate" if database in ABRICATE_NATIVE_DATABASES else database


def parse_abricate_results(path: Path, sample_dir: Path, sample_map: dict[str, str]) -> list[dict[str, str]]:
    rows = []
    database = database_from_abricate_path(path, sample_dir)
    for row in read_table(path):
        sample = first_value(row, ["#FILE", "FILE", "file", "sample_id", "sample"], "")
        feature_id = first_value(row, ["GENE", "gene", "feature_id", "id"], "")
        if not sample or not feature_id:
            continue
        rows.append(
            contract_row(
                sample,
                database,
                feature_id,
                feature_category=first_value(row, ["RESISTANCE", "PRODUCT", "product", "category"], ""),
                identity=first_value(row, ["%IDENTITY", "IDENTITY", "identity"], ""),
                coverage=first_value(row, ["%COVERAGE", "COVERAGE", "coverage"], ""),
                contig=first_value(row, ["SEQUENCE", "contig", "sequence"], ""),
                start=first_value(row, ["START", "start"], ""),
                end=first_value(row, ["END", "end"], ""),
                tool=tool_for_feature_table(database, row),
                sample_map=sample_map,
                product=first_value(row, ["PRODUCT", "product"], ""),
                drug_class=first_value(row, ["RESISTANCE", "resistance"], ""),
                source_table=str(path),
                source_file=first_value(row, ["#FILE", "FILE", "file"], ""),
                source_database=first_value(row, ["DATABASE", "database"], database),
                raw_feature_id=feature_id,
                raw_category=first_value(row, ["RESISTANCE", "PRODUCT", "product", "category"], ""),
            )
        )
    return rows


def parse_amrfinder_tables(path: Path, sample_map: dict[str, str]) -> list[dict[str, str]]:
    rows = []
    for row in read_table(path):
        sample = first_value(
            row,
            ["sample_id", "sample", "file", "#FILE", "assembly", "assembly_accession", "Source file"],
            path.stem,
        )
        feature_id = first_value(row, AMRFINDER_FEATURE_COLUMNS, "")
        if not feature_id:
            continue
        rows.append(
            contract_row(
                sample,
                "amrfinderplus",
                feature_id,
                feature_category=first_value(row, AMRFINDER_CATEGORY_COLUMNS, ""),
                identity=first_value(row, ["% Identity to reference sequence", "% identity", "identity"], ""),
                coverage=first_value(row, ["% Coverage of reference sequence", "% coverage", "coverage"], ""),
                contig=first_value(row, ["Contig id", "contig", "Sequence id"], ""),
                start=first_value(row, ["Start", "start"], ""),
                end=first_value(row, ["Stop", "End", "end"], ""),
                tool="amrfinderplus",
                sample_map=sample_map,
                feature_name=first_value(row, ["Element name", "Sequence name", "Name of closest sequence"], ""),
                feature_description=first_value(row, ["HMM description", "Name of closest sequence"], ""),
                feature_subcategory=first_value(row, ["Subclass", "Element subtype"], ""),
                mechanism=first_value(row, ["Element type", "Scope"], ""),
                drug_class=first_value(row, ["Class"], ""),
                product=first_value(row, ["Sequence name", "Name of closest sequence"], ""),
                sequence_id=first_value(row, ["Protein identifier", "Contig id", "Sequence id"], ""),
                strand=first_value(row, ["Strand", "strand"], ""),
                source_table=str(path),
                source_file=first_value(row, ["sample_id", "sample", "file", "Source file"], path.stem),
                source_database="NCBI AMRFinderPlus",
                raw_feature_id=feature_id,
                raw_category=first_value(row, AMRFINDER_CATEGORY_COLUMNS, ""),
                raw_method=first_value(row, ["Method"], ""),
                evidence_type=first_value(row, ["Method", "Scope"], ""),
            )
        )
    return rows


def parse_mlst_tables(path: Path, sample_map: dict[str, str]) -> list[dict[str, str]]:
    header = set(read_header(path))
    if not header.intersection({"sequence_type", "Sequence Type", "st", "ST", "allele_profile", "feature_id"}):
        return parse_mlst_raw_table(path, sample_map)

    rows = []
    for row in read_table(path):
        sample = first_value(
            row,
            ["Assembly Accession", "assembly_accession", "sample_id", "sample", "file", "#FILE"],
            "",
        )
        if not sample:
            continue
        scheme = first_value(row, ["scheme", "Scheme"], "")
        sequence_type = first_value(row, ["sequence_type", "Sequence Type"], "")
        st = first_value(row, ["st", "ST"], "")
        if is_missing_value(scheme):
            scheme = ""
        if is_placeholder_mlst_feature(sequence_type):
            sequence_type = ""
        if is_missing_value(st):
            st = ""
        primary_feature = first_value(row, ["feature_id"], "") or sequence_type or (f"ST_{st}" if st else "")
        if primary_feature and not is_placeholder_mlst_feature(primary_feature):
            rows.append(
                contract_row(
                    sample,
                    "mlst",
                    primary_feature,
                    feature_category="sequence_type",
                    identity=first_value(row, ["identity"], "100"),
                    tool="mlst",
                    sample_map=sample_map,
                    feature_name=primary_feature,
                    feature_subcategory=scheme,
                    source_table=str(path),
                    source_file=sample,
                    raw_feature_id=primary_feature,
                    raw_category="sequence_type",
                    evidence_type="sequence_type_call",
                )
            )
        if st:
            st_feature = f"ST_{st}"
            if st_feature != primary_feature and not is_placeholder_mlst_feature(st_feature):
                rows.append(
                    contract_row(
                        sample,
                        "mlst",
                        st_feature,
                        feature_category="sequence_type",
                        identity=first_value(row, ["identity"], "100"),
                        tool="mlst",
                        sample_map=sample_map,
                        feature_name=st_feature,
                        feature_subcategory=scheme,
                        source_table=str(path),
                        source_file=sample,
                        raw_feature_id=st,
                        raw_category="sequence_type",
                        evidence_type="sequence_type_call",
                    )
                )
        allele_profile = first_value(row, ["allele_profile", "Allele profile"], "")
        for allele in [part.strip() for part in allele_profile.split(";") if part.strip()]:
            match = re.fullmatch(r"([^()]+)\(([^()]+)\)", allele)
            if not match:
                continue
            locus, allele_number = match.groups()
            if is_missing_value(locus) or is_missing_value(allele_number):
                continue
            allele_feature = f"{locus}_{allele_number}"
            rows.append(
                contract_row(
                    sample,
                    "mlst",
                    allele_feature,
                    feature_category="mlst_allele",
                    identity="100",
                    tool="mlst",
                    sample_map=sample_map,
                    feature_name=allele_feature,
                    feature_subcategory=scheme,
                    product=locus,
                    source_table=str(path),
                    source_file=sample,
                    raw_feature_id=allele,
                    raw_category="mlst_allele",
                    evidence_type="allele_call",
                )
            )
    return rows


def parse_mlst_raw_table(path: Path, sample_map: dict[str, str]) -> list[dict[str, str]]:
    """Parse Torsten Seemann mlst tabular output.

    The native `mlst` command writes rows without a header:
    file, scheme, ST, then optional allele calls.  Unsupported organisms often
    report `-` for scheme and ST; those are intentionally retained as raw run
    evidence but are not converted into biological feature calls.
    """
    rows: list[dict[str, str]] = []
    if not path.exists() or path.stat().st_size == 0:
        return rows
    for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = [part.strip() for part in line.split("\t")]
        if len(parts) < 3:
            continue
        sample, scheme, st = parts[:3]
        if not sample:
            continue
        if is_missing_value(scheme):
            scheme = ""
        if not is_missing_value(st):
            st_feature = f"ST_{st}"
            if not is_placeholder_mlst_feature(st_feature):
                rows.append(
                    contract_row(
                        sample,
                        "mlst",
                        st_feature,
                        feature_category="sequence_type",
                        identity="100",
                        tool="mlst",
                        sample_map=sample_map,
                        feature_name=st_feature,
                        feature_subcategory=scheme,
                        source_table=str(path),
                        source_file=sample,
                        raw_feature_id=st,
                        raw_category="sequence_type",
                        evidence_type="sequence_type_call",
                    )
                )
        for allele in parts[3:]:
            match = re.fullmatch(r"([^()]+)\(([^()]+)\)", allele)
            if not match:
                continue
            locus, allele_number = match.groups()
            if is_missing_value(locus) or is_missing_value(allele_number):
                continue
            allele_feature = f"{locus}_{allele_number}"
            rows.append(
                contract_row(
                    sample,
                    "mlst",
                    allele_feature,
                    feature_category="mlst_allele",
                    identity="100",
                    tool="mlst",
                    sample_map=sample_map,
                    feature_name=allele_feature,
                    feature_subcategory=scheme,
                    product=locus,
                    source_table=str(path),
                    source_file=sample,
                    raw_feature_id=allele,
                    raw_category="mlst_allele",
                    evidence_type="allele_call",
                )
            )
    return rows


def _find_suffix_value(row: dict[str, str], suffixes: list[str]) -> str:
    for key, value in row.items():
        lowered = key.lower()
        if any(lowered.endswith(suffix.lower()) for suffix in suffixes):
            text = str(value or "").strip()
            if not is_missing_value(text):
                return text
    return ""


def parse_kleborate_tables(path: Path, sample_map: dict[str, str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for row in read_table(path):
        sample = first_value(row, ["strain", "Input_file_name", "sample_id"], path.stem)
        if not sample:
            continue

        st_value = _find_suffix_value(row, ["__mlst__ST"]) or first_value(row, ["ST"], "")
        st_template = "{}" if str(st_value).upper().startswith("ST") else "ST_{}"
        feature_specs = [
            ("sequence_type", st_value, st_template),
            ("virulence_score", _find_suffix_value(row, ["__virulence_score__virulence_score"]) or first_value(row, ["virulence_score"], ""), "virulence_score_{}"),
            ("resistance_score", _find_suffix_value(row, ["__resistance_score__resistance_score"]) or first_value(row, ["resistance_score"], ""), "resistance_score_{}"),
            ("resistance_class_count", _find_suffix_value(row, ["__resistance_class_count__num_resistance_classes"]) or first_value(row, ["num_resistance_classes"], ""), "resistance_classes_{}"),
            ("k_locus", _find_suffix_value(row, ["__kaptive__K_locus"]), "{}"),
            ("k_type", _find_suffix_value(row, ["__kaptive__K_type"]), "{}"),
            ("o_locus", _find_suffix_value(row, ["__kaptive__O_locus"]), "{}"),
            ("o_type", _find_suffix_value(row, ["__kaptive__O_type"]), "{}"),
            ("yersiniabactin", _find_suffix_value(row, ["__ybst__Yersiniabactin"]), "yersiniabactin_{}"),
            ("colibactin", _find_suffix_value(row, ["__cbst__Colibactin"]), "colibactin_{}"),
            ("aerobactin", _find_suffix_value(row, ["__abst__Aerobactin"]), "aerobactin_{}"),
            ("salmochelin", _find_suffix_value(row, ["__smst__Salmochelin"]), "salmochelin_{}"),
            ("rmpadc", _find_suffix_value(row, ["__rmst__RmpADC"]), "rmpadc_{}"),
            ("wzi", _find_suffix_value(row, ["__wzi__wzi"]), "wzi_{}"),
        ]
        for category, value, template in feature_specs:
            if is_missing_value(value):
                continue
            clean_value = str(value).split(";")[0].strip()
            if is_missing_value(clean_value):
                continue
            feature_id = template.format(clean_value.replace(" ", "_"))
            if is_placeholder_mlst_feature(feature_id):
                continue
            rows.append(
                contract_row(
                    sample,
                    "kleborate",
                    feature_id,
                    feature_category=category,
                    tool="kleborate",
                    sample_map=sample_map,
                    feature_name=feature_id,
                    feature_subcategory=category,
                    source_table=str(path),
                    source_file=sample,
                    raw_feature_id=value,
                    raw_category=category,
                    evidence_type="kleborate_call",
                )
            )

        for key, value in row.items():
            if not key.endswith("__amr__Bla_chr") and not key.endswith("__amr__Bla_acquired"):
                continue
            for gene in re.split(r"[;,]", str(value or "")):
                gene = gene.strip().replace("^", "")
                if is_missing_value(gene):
                    continue
                rows.append(
                    contract_row(
                        sample,
                        "kleborate",
                        gene,
                        feature_category="amr_marker",
                        tool="kleborate",
                        sample_map=sample_map,
                        feature_name=gene,
                        source_table=str(path),
                        source_file=sample,
                        raw_feature_id=gene,
                        raw_category=key,
                        evidence_type="kleborate_amr_marker",
                    )
                )
    return rows


def _mobsuite_values(value: str) -> list[str]:
    values: list[str] = []
    for part in re.split(r"[;,]", str(value or "")):
        text = part.strip()
        if "|" in text:
            text = text.split("|")[-1].strip()
        text = re.sub(r"\s+", "_", text)
        if not is_missing_value(text) and text not in values:
            values.append(text)
    return values


def _mobsuite_sample(value: str, fallback: str) -> str:
    text = first_value({"sample": value}, ["sample"], fallback)
    if ":" in text:
        text = text.split(":", 1)[0]
    return clean_sample_id(text)


def _mobsuite_category(value: str, default: str = "mobsuite_feature") -> str:
    text = str(value or default).strip()
    if is_missing_value(text):
        text = default
    return re.sub(r"[^A-Za-z0-9_]+", "_", text).strip("_").lower()


def _append_mobsuite_feature(
    rows: list[dict[str, str]],
    row: dict[str, str],
    path: Path,
    sample: str,
    sample_map: dict[str, str],
    feature_id: str,
    category: str,
    evidence_type: str,
    raw_value: str = "",
) -> None:
    category = _mobsuite_category(category)
    if is_missing_value(feature_id) or is_missing_value(category):
        return
    rows.append(
        contract_row(
            sample,
            "mobsuite",
            feature_id,
            feature_category=category,
            identity=first_value(row, ["pident", "identity"], ""),
            coverage=first_value(row, ["qcovs", "qcovhsp", "coverage"], ""),
            contig=first_value(row, ["contig_id", "sseqid", "SEQUENCE", "sequence"], ""),
            start=first_value(row, ["contig_start", "sstart", "START", "start"], ""),
            end=first_value(row, ["contig_end", "send", "END", "end"], ""),
            tool="mobsuite",
            sample_map=sample_map,
            feature_name=feature_id,
            feature_subcategory=category,
            source_table=str(path),
            source_file=sample,
            source_database="MOB-suite",
            raw_feature_id=raw_value or feature_id,
            raw_category=category,
            evidence_type=evidence_type,
            confidence=first_value(row, ["mash_neighbor_distance", "evalue", "bitscore"], ""),
            notes=first_value(row, ["filtering_reason", "mash_nearest_neighbor"], ""),
        )
    )


def parse_mobsuite_tables(path: Path, sample_map: dict[str, str]) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for row in read_table(path):
        sample = _mobsuite_sample(first_value(row, ["sample_id", "sample", "file", "#FILE"], path.stem), path.stem)
        if not sample:
            continue

        generic_feature = first_value(row, ["GENE", "gene", "feature_id", "id"], "")
        if generic_feature:
            _append_mobsuite_feature(
                rows,
                row,
                path,
                sample,
                sample_map,
                generic_feature,
                first_value(row, ["PRODUCT", "RESISTANCE", "product", "category"], "mobsuite_feature"),
                "mobsuite_table_feature",
                raw_value=generic_feature,
            )
            continue

        biomarker = first_value(row, ["biomarker"], "")
        for feature in _mobsuite_values(first_value(row, ["qseqid"], "")):
            _append_mobsuite_feature(rows, row, path, sample, sample_map, feature, biomarker, "mobsuite_biomarker", raw_value=first_value(row, ["qseqid"], ""))

        for column, category, evidence_type in [
            ("rep_type(s)", "replicon", "mobsuite_contig_report"),
            ("relaxase_type(s)", "relaxase", "mobsuite_contig_report"),
            ("mpf_type", "mate_pair_formation", "mobsuite_contig_report"),
            ("orit_type(s)", "origin_of_transfer", "mobsuite_contig_report"),
            ("predicted_mobility", "mobility", "mobsuite_contig_report"),
            ("primary_cluster_id", "plasmid_cluster", "mobsuite_contig_report"),
            ("secondary_cluster_id", "plasmid_cluster_secondary", "mobsuite_contig_report"),
            ("predicted_host_range_overall_name", "predicted_host_range", "mobsuite_contig_report"),
            ("observed_host_range_ncbi_name", "observed_host_range", "mobsuite_contig_report"),
            ("mash_nearest_neighbor", "mash_neighbor", "mobsuite_contig_report"),
        ]:
            for feature in _mobsuite_values(first_value(row, [column], "")):
                _append_mobsuite_feature(rows, row, path, sample, sample_map, feature, category, evidence_type, raw_value=first_value(row, [column], ""))

        molecule_type = first_value(row, ["molecule_type"], "")
        if not is_missing_value(molecule_type):
            _append_mobsuite_feature(
                rows,
                row,
                path,
                sample,
                sample_map,
                f"molecule_type_{molecule_type}",
                "molecule_type",
                "mobsuite_contig_report",
                raw_value=molecule_type,
            )

        mge_feature = first_value(row, ["mge_type", "mge_subtype", "mge_id"], "")
        mge_category = first_value(row, ["mge_subtype"], "mobile_genetic_element")
        for feature in _mobsuite_values(mge_feature):
            _append_mobsuite_feature(rows, row, path, sample, sample_map, feature, mge_category, "mobsuite_mge_report", raw_value=mge_feature)
    return rows


def discover_raw_feature_databases(sample_dir: Path) -> set[str]:
    raw_databases: set[str] = set()
    abricate_dirs_by_database = {
        "amr": [
            sample_dir / "abricate",
            sample_dir / "amr",
            sample_dir / "tool_results" / "abricate" / "ncbi",
        ],
        "vfdb": [
            sample_dir / "vfdb",
            sample_dir / "tool_results" / "abricate" / "vfdb",
        ],
        "plasmidfinder": [
            sample_dir / "plasmidfinder",
            sample_dir / "tool_results" / "abricate" / "plasmidfinder",
        ],
    }
    for database, dirs in abricate_dirs_by_database.items():
        if any(path.exists() and any(child.is_file() for child in path.rglob("*")) for path in dirs):
            raw_databases.add(database)
    mlst_dirs = [
        sample_dir / "mlst",
        sample_dir / "tool_results" / "mlst",
    ]
    if any(path.exists() and any(child.is_file() for child in path.rglob("*")) for path in mlst_dirs):
        raw_databases.add("mlst")
    integronfinder_dirs = [
        sample_dir / "integronfinder",
        sample_dir / "tool_results" / "integronfinder",
    ]
    if any(path.exists() and any(child.is_file() for child in path.rglob("*")) for path in integronfinder_dirs):
        raw_databases.add("integronfinder")
    for database in OPTIONAL_TABLE_DATABASES:
        dirs = [
            sample_dir / database,
            sample_dir / database / "tables",
            sample_dir / "tool_results" / database,
        ]
        for directory in dirs:
            if directory.exists() and any(path.is_file() for path in directory.rglob("*")):
                raw_databases.add(database)
                break
    return raw_databases


def optional_table_paths(sample_dir: Path) -> list[Path]:
    paths: list[Path] = []
    skip_suffixes = {
        "_collection_status.tsv",
        "_database_setup_status.tsv",
        "_warning.txt",
        "module_status.tsv",
        "all_features.tsv",
    }
    for database in OPTIONAL_TABLE_DATABASES:
        for directory in [
            sample_dir / database,
            sample_dir / database / "tables",
            sample_dir / "tool_results" / database,
        ]:
            if not directory.exists():
                continue
            for path in sorted(directory.rglob("*")):
                if not path.is_file() or path.suffix.lower() not in {".tsv", ".tab", ".csv"}:
                    continue
                if any(part in {"analysis", "figures", "merged_output", "panr2_inputs", "raw"} for part in path.parts):
                    continue
                if any(path.name.endswith(suffix) for suffix in skip_suffixes):
                    continue
                if path.name.endswith("_tidy.csv"):
                    continue
                if path.name.endswith(".features.tsv"):
                    continue
                paths.append(path)
    return paths


def discover_feature_rows(sample_dir: Path) -> list[dict[str, str]]:
    sample_map = load_sample_map(sample_dir)
    rows: list[dict[str, str]] = []
    seen_paths: set[Path] = set()

    native_handoff_dirs = [
        sample_dir / "tool_results" / "integronfinder" / "panr2_inputs",
        sample_dir / "tool_results" / "mobileelementfinder" / "panr2_inputs",
    ]
    for handoff_dir in native_handoff_dirs:
        for path in sorted(handoff_dir.glob("*_results.*")):
            if path.suffix.lower() not in {".tab", ".tsv", ".csv"}:
                continue
            resolved = path.resolve()
            if resolved in seen_paths:
                continue
            seen_paths.add(resolved)
            rows.extend(parse_abricate_results(path, sample_dir, sample_map))

    for path in sorted(sample_dir.rglob("*_results.*")):
        if path.suffix.lower() not in {".tab", ".tsv", ".csv"}:
            continue
        if any(part in {"analysis", "figures", "merged_output", "panr2_inputs"} for part in path.parts):
            continue
        resolved = path.resolve()
        if resolved in seen_paths:
            continue
        seen_paths.add(resolved)
        rows.extend(parse_abricate_results(path, sample_dir, sample_map))

    for path in optional_table_paths(sample_dir):
        resolved = path.resolve()
        if resolved in seen_paths:
            continue
        seen_paths.add(resolved)
        if "kleborate" in path.parts:
            rows.extend(parse_kleborate_tables(path, sample_map))
        elif "mobsuite" in path.parts:
            rows.extend(parse_mobsuite_tables(path, sample_map))
        else:
            rows.extend(parse_abricate_results(path, sample_dir, sample_map))

    for path in sorted((sample_dir / "amrfinderplus").rglob("*")):
        if path.is_file() and path.suffix.lower() in {".tsv", ".tab", ".csv"}:
            rows.extend(parse_amrfinder_tables(path, sample_map))

    for mlst_dir in [sample_dir / "mlst", sample_dir / "tool_results" / "mlst"]:
        for path in sorted(mlst_dir.rglob("*")):
            if path.is_file() and path.suffix.lower() in {".tsv", ".tab", ".csv"}:
                rows.extend(parse_mlst_tables(path, sample_map))

    unique: dict[tuple[str, str, str, str, str, str], dict[str, str]] = {}
    for row in rows:
        if not row["feature_id"]:
            continue
        key = (
            row["sample_id"],
            row["assembly_accession"],
            row["database"],
            row["feature_id"],
            row["contig"],
            row["start"],
        )
        unique[key] = row
    return list(unique.values())


def write_feature_tables(sample_dir: Path, out_dir: Path) -> dict[str, str]:
    rows = discover_feature_rows(sample_dir)
    raw_databases = discover_raw_feature_databases(sample_dir)
    features_dir = out_dir / "features"
    features_dir.mkdir(parents=True, exist_ok=True)

    by_database: dict[str, list[dict[str, str]]] = {}
    for row in rows:
        by_database.setdefault(row["database"], []).append(row)

    written: dict[str, str] = {}
    for database, database_rows in sorted(by_database.items()):
        path = features_dir / f"{database}.features.tsv"
        write_rows(path, database_rows, FEATURE_COLUMNS)
        written[database] = str(path)

    for database in sorted(raw_databases - set(by_database)):
        path = features_dir / f"{database}.features.tsv"
        write_rows(path, [], FEATURE_COLUMNS)
        written[database] = str(path)

    all_path = features_dir / "all_features.tsv"
    write_rows(all_path, rows, FEATURE_COLUMNS)
    written["all_features"] = str(all_path)
    return written


def _is_numeric_or_blank(value: str) -> bool:
    text = str(value or "").strip()
    if not text:
        return True
    try:
        float(text)
        return True
    except ValueError:
        return False


def _is_int_or_blank(value: str) -> bool:
    text = str(value or "").strip()
    if not text:
        return True
    try:
        int(float(text))
        return True
    except ValueError:
        return False


def _as_float(value: str) -> float | None:
    try:
        text = str(value or "").strip()
        return float(text) if text else None
    except ValueError:
        return None


def _as_int(value: str) -> int | None:
    try:
        text = str(value or "").strip()
        return int(float(text)) if text else None
    except ValueError:
        return None


def _feature_key(row: dict[str, str]) -> tuple[str, str, str, str, str, str, str]:
    return (
        row.get("sample_id", ""),
        row.get("assembly_accession", ""),
        row.get("database", ""),
        row.get("feature_id", ""),
        row.get("contig", ""),
        row.get("start", ""),
        row.get("end", ""),
    )


def validate_contract_row(row: dict[str, str], metadata_accessions: set[str]) -> list[str]:
    errors = []
    for column in ["sample_id", "assembly_accession", "database", "feature_id", "tool"]:
        if not str(row.get(column, "")).strip():
            errors.append(f"{column}_empty")
    if row.get("presence", "") not in {"0", "1", 0, 1}:
        errors.append("presence_not_0_or_1")
    for column in ["identity", "coverage"]:
        if not _is_numeric_or_blank(row.get(column, "")):
            errors.append(f"{column}_not_numeric")
    for column in ["start", "end"]:
        if not _is_int_or_blank(row.get(column, "")):
            errors.append(f"{column}_not_integer")
    start = _as_int(row.get("start", ""))
    end = _as_int(row.get("end", ""))
    if start is not None and end is not None and start > end:
        errors.append("start_greater_than_end")
    database = str(row.get("database", "")).strip()
    if database and database not in KNOWN_DATABASES:
        errors.append("database_not_in_vocabulary")
    accession = str(row.get("assembly_accession", "")).strip()
    if metadata_accessions and accession and accession not in metadata_accessions:
        errors.append("assembly_accession_unmatched")
    return errors


def validate_feature_tables(sample_dir: Path, out_dir: Path) -> dict[str, str]:
    metadata_accessions = load_metadata_accessions(sample_dir)
    features_dir = out_dir / "features"
    manifest_dir = out_dir / "manifest"
    manifest_dir.mkdir(parents=True, exist_ok=True)

    validation_rows = []
    unmatched_rows = []
    invalid_rows = []
    duplicate_rows = []
    total_rows = 0
    files_checked = 0
    databases_seen: set[str] = set()
    samples_seen: set[str] = set()
    duplicate_counter: Counter[tuple[str, str, str, str, str, str, str]] = Counter()

    feature_paths = sorted(features_dir.glob("*.features.tsv"))
    for path in feature_paths:
        header = read_header(path)
        rows = read_table(path)
        files_checked += 1
        missing_columns = [col for col in CONTRACT_COLUMNS if col not in header]
        status = "fail" if missing_columns else "ok"
        validation_rows.append(
            {
                "file": str(path),
                "status": status,
                "rows": str(len(rows)),
                "missing_columns": ",".join(missing_columns),
            }
        )

    all_features_path = features_dir / "all_features.tsv"
    row_level_paths = [all_features_path] if all_features_path.exists() else feature_paths
    for path in row_level_paths:
        rows = read_table(path)
        for row in rows:
            total_rows += 1
            duplicate_counter[_feature_key(row)] += 1
            database = row.get("database", "")
            sample = row.get("sample_id", "")
            accession = row.get("assembly_accession", "")
            if database:
                databases_seen.add(database)
            if sample:
                samples_seen.add(sample)
            errors = validate_contract_row(row, metadata_accessions)
            if errors:
                invalid = dict(row)
                invalid["validation_errors"] = ";".join(errors)
                invalid_rows.append(invalid)
            if metadata_accessions and accession and accession not in metadata_accessions:
                unmatched_rows.append(row)

    for key, count in duplicate_counter.items():
        if count > 1:
            duplicate_rows.append({
                "sample_id": key[0],
                "assembly_accession": key[1],
                "database": key[2],
                "feature_id": key[3],
                "contig": key[4],
                "start": key[5],
                "end": key[6],
                "duplicate_count": str(count),
            })

    report_path = manifest_dir / "schema_validation_report.csv"
    if not validation_rows:
        validation_rows.append({"file": str(features_dir), "status": "warning", "rows": "0", "missing_columns": ""})
    write_rows(report_path, validation_rows, ["file", "status", "rows", "missing_columns"], delimiter=",")

    unmatched_path = manifest_dir / "unmatched_features.csv"
    write_rows(unmatched_path, unmatched_rows, FEATURE_COLUMNS, delimiter=",")

    duplicate_path = manifest_dir / "duplicate_features.csv"
    write_rows(
        duplicate_path,
        duplicate_rows,
        ["sample_id", "assembly_accession", "database", "feature_id", "contig", "start", "end", "duplicate_count"],
        delimiter=",",
    )

    invalid_path = manifest_dir / "invalid_feature_rows.csv"
    write_rows(invalid_path, invalid_rows, FEATURE_COLUMNS + ["validation_errors"], delimiter=",")

    summary_path = manifest_dir / "schema_validation_summary.txt"
    summary_path.write_text(
        "\n".join(
            [
                f"feature_files_checked={files_checked}",
                f"feature_rows={total_rows}",
                f"databases_seen={','.join(sorted(databases_seen))}",
                f"samples_seen={len(samples_seen)}",
                f"metadata_accessions={len(metadata_accessions)}",
                f"unmatched_feature_rows={len(unmatched_rows)}",
                f"invalid_feature_rows={len(invalid_rows)}",
                f"duplicate_feature_rows={len(duplicate_rows)}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    return {
        "schema_validation_report": str(report_path),
        "unmatched_features": str(unmatched_path),
        "duplicate_features": str(duplicate_path),
        "invalid_feature_rows": str(invalid_path),
        "schema_validation_summary": str(summary_path),
    }


def load_metadata_rows(sample_dir: Path) -> list[dict[str, str]]:
    for path in [
        sample_dir / "metadata_output" / "ncbi_clean_qc_pass.csv",
        sample_dir / "metadata_output" / "ncbi_clean.csv",
        sample_dir / "metadata_output" / "fetchm2_clean.csv",
        sample_dir / "metadata_output" / "ncbi_enriched.csv",
        sample_dir / "panr2_inputs" / "metadata" / "ncbi_clean_qc_pass.csv",
        sample_dir / "panr2_inputs" / "metadata" / "ncbi_clean.csv",
        sample_dir / "panr2_inputs" / "metadata" / "fetchm2_clean.csv",
        sample_dir / "panr2_inputs" / "metadata" / "ncbi_enriched.csv",
    ]:
        rows = read_table(path)
        if rows:
            return rows
    return []


def metadata_accession(row: dict[str, str]) -> str:
    return first_value(row, ["Assembly Accession", "assembly_accession", "Assembly"], "")


def normalize_metadata_rows(metadata_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    normalized = []
    for row in metadata_rows:
        out = {
            "assembly_accession": metadata_accession(row),
            "sample_id": first_value(row, ["sample_id", "Assembly Accession", "assembly_accession"], metadata_accession(row)),
        }
        for alias, candidates in METADATA_ALIASES.items():
            out[alias] = first_value(row, candidates, "")
        normalized.append(out)
    return normalized


def _classify_metadata_column(values: list[str]) -> str:
    non_missing = [value for value in values if value]
    if not non_missing:
        return "mostly_missing"
    numeric = sum(1 for value in non_missing if _is_numeric_or_blank(value))
    years = sum(1 for value in non_missing if re.fullmatch(r"\d{4}", value))
    unique = len(set(non_missing))
    if years / len(non_missing) >= 0.8:
        return "year"
    if numeric / len(non_missing) >= 0.8:
        return "numeric"
    if unique == len(non_missing) and unique > 10:
        return "identifier"
    return "categorical"


def _float_or_none(value: str | None) -> float | None:
    text = str(value or "").strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _format_fraction(value: float | None) -> str:
    return "" if value is None else f"{value:.4f}"


def _flag_string(flags: list[str]) -> str:
    return ";".join(sorted({flag for flag in flags if flag}))


def _interpretation_label(effect_size: str, q_value: str = "", warning_flags: str = "") -> str:
    flags = {flag for flag in str(warning_flags or "").split(";") if flag}
    effect = _float_or_none(effect_size) or 0.0
    q = _float_or_none(q_value)
    if "single_bioproject_dominance" in flags or "bioproject_bias_warning" in flags:
        return "bioproject_bias_warning"
    if "small_group" in flags or "low_sample_warning" in flags:
        return "low_sample_warning"
    if "sparse_metadata_warning" in flags:
        return "sparse_metadata_warning"
    if q is not None and q <= 0.05 and effect >= 0.5:
        return "strong_supported"
    if effect >= 0.5 and "screening_no_p_value" not in flags:
        return "moderate_supported"
    if effect >= 0.5:
        return "exploratory"
    return "exploratory"


def _bioproject_by_sample(metadata_rows: list[dict[str, str]]) -> dict[str, str]:
    normalized = normalize_metadata_rows(metadata_rows)
    return {
        row["assembly_accession"]: row.get("bioproject", "")
        for row in normalized
        if row.get("assembly_accession") and row.get("bioproject")
    }


def _top_project_summary(samples: set[str], project_by_sample: dict[str, str]) -> dict[str, str]:
    projects = [project_by_sample.get(sample, "") for sample in samples if project_by_sample.get(sample, "")]
    counts = Counter(projects)
    top_project, top_count = counts.most_common(1)[0] if counts else ("", 0)
    fraction = top_count / len(projects) if projects else None
    status = "WARNING_DOMINATED" if fraction is not None and len(projects) >= 3 and fraction >= 0.8 else "PASS"
    warning = "single_bioproject_dominance" if status == "WARNING_DOMINATED" else ""
    return {
        "samples_evaluated": str(len(projects)),
        "n_bioprojects": str(len(counts)),
        "largest_bioproject": top_project,
        "largest_bioproject_count": str(top_count),
        "largest_bioproject_fraction": _format_fraction(fraction),
        "status": status,
        "warning": warning,
    }


def _dominance_summary(samples: set[str], values_by_sample: dict[str, str]) -> dict[str, str]:
    values = [values_by_sample.get(sample, "") for sample in samples if values_by_sample.get(sample, "")]
    counts = Counter(values)
    top_value, top_count = counts.most_common(1)[0] if counts else ("", 0)
    fraction = top_count / len(values) if values else None
    return {
        "samples_evaluated": str(len(values)),
        "unique_values": str(len(counts)),
        "dominant_value": top_value,
        "dominant_count": str(top_count),
        "dominant_fraction": _format_fraction(fraction),
    }


def _dominance_warning(summary: dict[str, str], flag: str, min_samples: int = 3, threshold: float = 0.8) -> str:
    samples = int(summary.get("samples_evaluated", "0") or 0)
    fraction = _float_or_none(summary.get("dominant_fraction", ""))
    if samples >= min_samples and fraction is not None and fraction >= threshold:
        return flag
    return ""


def _st_by_sample(rows: list[dict[str, str]]) -> dict[str, str]:
    st = {}
    for row in rows:
        if row.get("database") != "mlst" or row.get("presence", "1") != "1":
            continue
        feature_id = row.get("feature_id", "")
        category = row.get("feature_category", "")
        sample = row.get("assembly_accession", "") or row.get("sample_id", "")
        if sample and (category == "sequence_type" or feature_id.startswith("ST_")) and feature_id.startswith("ST_"):
            st.setdefault(sample, feature_id)
    return st


def _ani_cluster_by_sample(out_dir: Path) -> dict[str, str]:
    cluster_by_sample = {}
    for path in [
        out_dir / "ani" / "analysis" / "duplicate_clusters.csv",
        out_dir / "ani" / "duplicate_clusters.csv",
        out_dir.parent / "ani" / "analysis" / "duplicate_clusters.csv",
        out_dir.parent / "ani" / "duplicate_clusters.csv",
    ]:
        for row in read_table(path):
            raw_sample = first_value(row, ["genome", "assembly_accession", "sample_id"], "")
            cluster = first_value(row, ["ani_cluster", "cluster", "feature_id"], "")
            if not raw_sample or not cluster:
                continue
            sample = extract_accession(raw_sample) or clean_sample_id(raw_sample)
            cluster_by_sample[sample] = cluster
    return cluster_by_sample


def _metadata_group_by_sample(metadata_rows: list[dict[str, str]], column: str) -> dict[str, str]:
    normalized = normalize_metadata_rows(metadata_rows)
    return {
        row["assembly_accession"]: row.get(column, "")
        for row in normalized
        if row.get("assembly_accession") and row.get(column)
    }


def write_metadata_analysis(sample_dir: Path, out_dir: Path) -> dict[str, str]:
    metadata_rows = load_metadata_rows(sample_dir)
    analysis_dir = out_dir / "metadata_feature_analysis"
    prevalence_dir = analysis_dir / "prevalence_tables"
    prevalence_dir.mkdir(parents=True, exist_ok=True)

    normalized = normalize_metadata_rows(metadata_rows)
    normalized_fields = ["sample_id", "assembly_accession", *METADATA_ALIASES.keys()]
    normalized_path = write_rows(analysis_dir / "metadata_normalized_for_analysis.tsv", normalized, normalized_fields)

    audit_rows = []
    eligibility_rows = []
    recommended: dict[str, list[str]] = defaultdict(list)
    if metadata_rows:
        metadata_columns = sorted({column for row in metadata_rows for column in row})
        total = len(metadata_rows)
        for column in metadata_columns:
            values = [str(row.get(column, "") or "").strip() for row in metadata_rows]
            non_missing_values = [value for value in values if value and value.lower() not in {"nan", "none", "unknown"}]
            counts = Counter(non_missing_values)
            top_value, top_count = counts.most_common(1)[0] if counts else ("", 0)
            data_type = _classify_metadata_column(non_missing_values)
            missing_fraction = 1 - (len(non_missing_values) / total if total else 0)
            largest_group_fraction = top_count / len(non_missing_values) if non_missing_values else 0
            unique_values = len(counts)
            eligible = (
                data_type in {"categorical", "year", "numeric"}
                and len(non_missing_values) >= 5
                and unique_values >= 2
                and largest_group_fraction < 0.95
            )
            reason = "eligible" if eligible else "sparse_constant_identifier_or_dominated"
            audit_rows.append({
                "column": column,
                "standardized_name": next((alias for alias, candidates in METADATA_ALIASES.items() if column in candidates), ""),
                "data_type": data_type,
                "non_missing_count": str(len(non_missing_values)),
                "missing_fraction": f"{missing_fraction:.4f}",
                "unique_values": str(unique_values),
                "top_value": top_value,
                "top_value_fraction": f"{largest_group_fraction:.4f}",
                "recommended_for_analysis": str(eligible).lower(),
                "reason": reason,
            })
            eligibility_rows.append({
                "metadata_column": column,
                "data_type": data_type,
                "non_missing_count": str(len(non_missing_values)),
                "missing_fraction": f"{missing_fraction:.4f}",
                "unique_values": str(unique_values),
                "largest_group_fraction": f"{largest_group_fraction:.4f}",
                "eligible": str(eligible).lower(),
                "reason": reason,
            })
            if eligible:
                group = next((alias for alias, candidates in METADATA_ALIASES.items() if column in candidates), "other")
                recommended[group].append(column)

    audit_path = write_rows(
        analysis_dir / "fetchm2_metadata_audit.tsv",
        audit_rows,
        [
            "column",
            "standardized_name",
            "data_type",
            "non_missing_count",
            "missing_fraction",
            "unique_values",
            "top_value",
            "top_value_fraction",
            "recommended_for_analysis",
            "reason",
        ],
    )
    eligibility_path = write_rows(
        analysis_dir / "metadata_column_eligibility.tsv",
        eligibility_rows,
        [
            "metadata_column",
            "data_type",
            "non_missing_count",
            "missing_fraction",
            "unique_values",
            "largest_group_fraction",
            "eligible",
            "reason",
        ],
    )
    usability_path = write_rows(
        analysis_dir / "metadata_usability_summary.tsv",
        [
            {
                "metadata_column": row["metadata_column"],
                "non_missing_count": row["non_missing_count"],
                "missing_fraction": row["missing_fraction"],
                "unique_values": row["unique_values"],
                "largest_group_fraction": row["largest_group_fraction"],
                "recommended_for_analysis": row["eligible"],
                "reason": row["reason"],
            }
            for row in eligibility_rows
        ],
        [
            "metadata_column",
            "non_missing_count",
            "missing_fraction",
            "unique_values",
            "largest_group_fraction",
            "recommended_for_analysis",
            "reason",
        ],
    )
    recommended_path = analysis_dir / "recommended_metadata_columns.txt"
    recommended_path.parent.mkdir(parents=True, exist_ok=True)
    with recommended_path.open("w", encoding="utf-8") as handle:
        for group in sorted(set(METADATA_ALIASES) | set(recommended)):
            handle.write(f"[{group}]\n")
            for column in sorted(recommended.get(group, [])):
                handle.write(f"{column}\n")
            handle.write("\n")
    return {
        "metadata_normalized": normalized_path,
        "metadata_audit": audit_path,
        "metadata_column_eligibility": eligibility_path,
        "metadata_usability_summary": usability_path,
        "recommended_metadata_columns": str(recommended_path),
    }


def _feature_label(row: dict[str, str]) -> str:
    return f"{row.get('database', '')}|{row.get('feature_id', '')}"


def feature_presence(rows: list[dict[str, str]]) -> dict[tuple[str, str], set[str]]:
    presence: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in rows:
        if row.get("presence", "1") != "1":
            continue
        database = row.get("database", "")
        feature_id = row.get("feature_id", "")
        sample = row.get("assembly_accession", "") or row.get("sample_id", "")
        if database and feature_id and sample:
            presence[(database, feature_id)].add(sample)
    return presence


def _supporting_samples_for_finding(
    finding: dict[str, str],
    normalized_metadata_rows: list[dict[str, str]],
    presence: dict[tuple[str, str], set[str]],
) -> set[str]:
    metadata_column = finding.get("metadata_column", "")
    metadata_value = finding.get("metadata_value", "")
    group_samples = {
        row["assembly_accession"]
        for row in normalized_metadata_rows
        if row.get("assembly_accession") and (not metadata_column or row.get(metadata_column, "") == metadata_value)
    }
    if finding.get("feature_id") == "database_burden":
        return group_samples
    feature_samples = presence.get((finding.get("database", ""), finding.get("feature_id", "")), set())
    return feature_samples & group_samples if group_samples else set(feature_samples)


def _annotate_top_findings(
    top_findings: list[dict[str, str]],
    normalized_metadata_rows: list[dict[str, str]],
    presence: dict[tuple[str, str], set[str]],
    project_by_sample: dict[str, str],
    st_by_sample: dict[str, str] | None = None,
    ani_cluster_by_sample: dict[str, str] | None = None,
) -> list[dict[str, str]]:
    st_by_sample = st_by_sample or {}
    ani_cluster_by_sample = ani_cluster_by_sample or {}
    annotated = []
    for finding in top_findings:
        support_samples = _supporting_samples_for_finding(finding, normalized_metadata_rows, presence)
        project_summary = _top_project_summary(support_samples, project_by_sample)
        st_summary = _dominance_summary(support_samples, st_by_sample)
        ani_summary = _dominance_summary(support_samples, ani_cluster_by_sample)
        flags = []
        if finding.get("warning"):
            flags.extend(flag.strip() for flag in finding["warning"].split(";") if flag.strip())
        if not finding.get("p_value") and not finding.get("q_value"):
            flags.append("screening_no_p_value")
        if 0 < len(support_samples) < 5:
            flags.append("low_sample_warning")
        if project_summary["warning"]:
            flags.append(project_summary["warning"])
        st_warning = _dominance_warning(st_summary, "single_ST_dominance")
        ani_warning = _dominance_warning(ani_summary, "single_ani_cluster_dominance")
        if st_warning:
            flags.append(st_warning)
        if ani_warning:
            flags.append(ani_warning)
        lineage_flags = [flag for flag in [st_warning, ani_warning] if flag]
        if not st_by_sample and not ani_cluster_by_sample:
            lineage_flags.append("insufficient_lineage_data")
        elif st_by_sample and len({value for value in st_by_sample.values() if value}) <= 1 and len(st_by_sample) >= 3:
            lineage_flags.append("low_lineage_diversity")
        flags.extend(lineage_flags)
        warning_flags = _flag_string(flags)
        out = dict(finding)
        out.update({
            "supporting_samples": str(len(support_samples)),
            "largest_bioproject": project_summary["largest_bioproject"],
            "largest_bioproject_fraction": project_summary["largest_bioproject_fraction"],
            "dominant_ST": st_summary["dominant_value"],
            "dominant_ST_fraction": st_summary["dominant_fraction"],
            "dominant_ani_cluster": ani_summary["dominant_value"],
            "dominant_ani_cluster_fraction": ani_summary["dominant_fraction"],
            "lineage_warning_flags": _flag_string(lineage_flags),
            "warning_flags": warning_flags,
            "interpretation_label": _interpretation_label(finding.get("effect_size", ""), finding.get("q_value", ""), warning_flags),
        })
        annotated.append(out)
    return annotated


def write_bioproject_bias_report(
    rows: list[dict[str, str]],
    metadata_rows: list[dict[str, str]],
    top_findings: list[dict[str, str]],
    out_dir: Path,
) -> str:
    analysis_dir = out_dir / "metadata_feature_analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    normalized = normalize_metadata_rows(metadata_rows)
    project_by_sample = _bioproject_by_sample(metadata_rows)
    presence = feature_presence(rows)
    output_rows = []

    all_samples = {row.get("assembly_accession", "") for row in normalized if row.get("assembly_accession")}
    summary = _top_project_summary(all_samples, project_by_sample)
    output_rows.append({
        "row_type": "dataset_summary",
        "database": "",
        "feature_id": "",
        "metadata_column": "bioproject",
        "metadata_value": "",
        **summary,
    })

    for (database, feature_id), present_samples in sorted(presence.items()):
        if len(present_samples) < 3:
            continue
        summary = _top_project_summary(present_samples, project_by_sample)
        output_rows.append({
            "row_type": "feature_project_dominance",
            "database": database,
            "feature_id": feature_id,
            "metadata_column": "bioproject",
            "metadata_value": "",
            **summary,
        })

    for finding in top_findings:
        support_samples = _supporting_samples_for_finding(finding, normalized, presence)
        summary = _top_project_summary(support_samples, project_by_sample)
        output_rows.append({
            "row_type": "top_finding_project_support",
            "database": finding.get("database", ""),
            "feature_id": finding.get("feature_id", ""),
            "metadata_column": finding.get("metadata_column", ""),
            "metadata_value": finding.get("metadata_value", ""),
            **summary,
        })

    return write_rows(
        analysis_dir / "bioproject_bias_report.tsv",
        output_rows,
        [
            "row_type",
            "database",
            "feature_id",
            "metadata_column",
            "metadata_value",
            "samples_evaluated",
            "n_bioprojects",
            "largest_bioproject",
            "largest_bioproject_count",
            "largest_bioproject_fraction",
            "status",
            "warning",
        ],
    )


def write_lineage_analysis(
    rows: list[dict[str, str]],
    metadata_rows: list[dict[str, str]],
    top_findings: list[dict[str, str]],
    out_dir: Path,
) -> dict[str, str]:
    analysis_dir = out_dir / "metadata_feature_analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    samples = sorted({metadata_accession(row) for row in metadata_rows if metadata_accession(row)} | {row.get("assembly_accession", "") for row in rows if row.get("assembly_accession", "")})
    st_by_sample = _st_by_sample(rows)
    ani_by_sample = _ani_cluster_by_sample(out_dir)
    project_by_sample = _bioproject_by_sample(metadata_rows)
    presence = feature_presence(rows)
    normalized = normalize_metadata_rows(metadata_rows)

    lineage_rows = []
    for sample in samples:
        lineage_rows.append({
            "assembly_accession": sample,
            "mlst_ST": st_by_sample.get(sample, ""),
            "ani_cluster": ani_by_sample.get(sample, ""),
            "bioproject": project_by_sample.get(sample, ""),
            "lineage_data_status": "available" if st_by_sample.get(sample) or ani_by_sample.get(sample) else "missing",
        })
    lineage_summary_path = write_rows(
        analysis_dir / "lineage_summary.tsv",
        lineage_rows,
        ["assembly_accession", "mlst_ST", "ani_cluster", "bioproject", "lineage_data_status"],
    )

    burden_rows = []
    for lineage_type, lineage_map in [("mlst_ST", st_by_sample), ("ani_cluster", ani_by_sample)]:
        grouped: dict[str, set[str]] = defaultdict(set)
        for sample, value in lineage_map.items():
            if value:
                grouped[value].add(sample)
        for lineage_value, members in sorted(grouped.items()):
            for (database, feature_id), present_samples in sorted(presence.items()):
                present = len(present_samples & members)
                if present == 0:
                    continue
                burden_rows.append({
                    "lineage_type": lineage_type,
                    "lineage_value": lineage_value,
                    "database": database,
                    "feature_id": feature_id,
                    "lineage_sample_count": str(len(members)),
                    "present_count": str(present),
                    "prevalence": f"{present / len(members):.4f}" if members else "0.0000",
                })
    lineage_burden_path = write_rows(
        analysis_dir / "lineage_feature_burden.tsv",
        burden_rows,
        ["lineage_type", "lineage_value", "database", "feature_id", "lineage_sample_count", "present_count", "prevalence"],
    )

    confounding_rows = []
    for metadata_column in ["country", "host", "sample_type", "isolation_source", "environment_medium", "collection_year", "bioproject"]:
        metadata_map = _metadata_group_by_sample(metadata_rows, metadata_column)
        for lineage_type, lineage_map in [("mlst_ST", st_by_sample), ("ani_cluster", ani_by_sample)]:
            grouped: dict[str, set[str]] = defaultdict(set)
            for sample, lineage_value in lineage_map.items():
                if lineage_value and metadata_map.get(sample):
                    grouped[lineage_value].add(sample)
            for lineage_value, members in sorted(grouped.items()):
                summary = _dominance_summary(members, metadata_map)
                warning = _dominance_warning(summary, "metadata_lineage_confounding")
                confounding_rows.append({
                    "metadata_column": metadata_column,
                    "lineage_type": lineage_type,
                    "lineage_value": lineage_value,
                    "samples_evaluated": summary["samples_evaluated"],
                    "dominant_metadata_value": summary["dominant_value"],
                    "dominant_metadata_fraction": summary["dominant_fraction"],
                    "status": "WARNING_CONFOUNDED" if warning else "PASS",
                    "warning": warning,
                })
    lineage_confounding_path = write_rows(
        analysis_dir / "lineage_metadata_confounding.tsv",
        confounding_rows,
        ["metadata_column", "lineage_type", "lineage_value", "samples_evaluated", "dominant_metadata_value", "dominant_metadata_fraction", "status", "warning"],
    )

    warning_rows = []
    for finding in top_findings:
        warning_rows.append({
            "finding_type": finding.get("finding_type", ""),
            "database": finding.get("database", ""),
            "feature_id": finding.get("feature_id", ""),
            "metadata_column": finding.get("metadata_column", ""),
            "metadata_value": finding.get("metadata_value", ""),
            "supporting_samples": finding.get("supporting_samples", ""),
            "dominant_ST": finding.get("dominant_ST", ""),
            "dominant_ST_fraction": finding.get("dominant_ST_fraction", ""),
            "dominant_ani_cluster": finding.get("dominant_ani_cluster", ""),
            "dominant_ani_cluster_fraction": finding.get("dominant_ani_cluster_fraction", ""),
            "lineage_warning_flags": finding.get("lineage_warning_flags", ""),
            "interpretation": "lineage_context_required" if finding.get("lineage_warning_flags") else "no_lineage_warning",
        })
    lineage_warnings_path = write_rows(
        analysis_dir / "lineage_adjusted_warnings.tsv",
        warning_rows,
        [
            "finding_type",
            "database",
            "feature_id",
            "metadata_column",
            "metadata_value",
            "supporting_samples",
            "dominant_ST",
            "dominant_ST_fraction",
            "dominant_ani_cluster",
            "dominant_ani_cluster_fraction",
            "lineage_warning_flags",
            "interpretation",
        ],
    )
    return {
        "lineage_summary": lineage_summary_path,
        "lineage_feature_burden": lineage_burden_path,
        "lineage_metadata_confounding": lineage_confounding_path,
        "lineage_adjusted_warnings": lineage_warnings_path,
    }


def write_feature_matrices(rows: list[dict[str, str]], metadata_rows: list[dict[str, str]], out_dir: Path, max_features: int = 300) -> dict[str, str]:
    matrix_dir = out_dir / "feature_matrices"
    matrix_dir.mkdir(parents=True, exist_ok=True)
    samples = sorted({metadata_accession(row) for row in metadata_rows if metadata_accession(row)} | {row.get("assembly_accession", "") for row in rows if row.get("assembly_accession", "")})
    presence = feature_presence(rows)
    ordered_features = sorted(presence, key=lambda item: (-len(presence[item]), item[0], item[1]))[:max_features]

    def write_matrix(path: Path, selected_features: list[tuple[str, str]]) -> str:
        fields = ["assembly_accession"] + [f"{database}|{feature_id}" for database, feature_id in selected_features]
        output_rows = []
        for sample in samples:
            out = {"assembly_accession": sample}
            for database, feature_id in selected_features:
                out[f"{database}|{feature_id}"] = "1" if sample in presence[(database, feature_id)] else "0"
            output_rows.append(out)
        return write_rows(path, output_rows, fields)

    outputs = {"all_feature_matrix": write_matrix(matrix_dir / "all_features_presence_absence.tsv", ordered_features)}
    by_database: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for feature in ordered_features:
        by_database[feature[0]].append(feature)
    for database, features in sorted(by_database.items()):
        outputs[f"{database}_matrix"] = write_matrix(matrix_dir / f"{database}_presence_absence.tsv", features)
    return outputs


def _feature_sets_by_sample(rows: list[dict[str, str]]) -> dict[str, set[tuple[str, str]]]:
    sample_features: dict[str, set[tuple[str, str]]] = defaultdict(set)
    for row in rows:
        if row.get("presence", "1") != "1":
            continue
        sample = row.get("assembly_accession", "") or row.get("sample_id", "")
        database = row.get("database", "")
        feature_id = row.get("feature_id", "")
        if sample and database and feature_id:
            sample_features[sample].add((database, feature_id))
    return sample_features


def _shannon(values: list[int]) -> float:
    total = sum(values)
    if total == 0:
        return 0.0
    score = 0.0
    for value in values:
        if value <= 0:
            continue
        p = value / total
        score -= p * math.log(p)
    return score


def _simpson(values: list[int]) -> float:
    total = sum(values)
    if total == 0:
        return 0.0
    return 1.0 - sum((value / total) ** 2 for value in values if value > 0)


def write_diversity_analysis(
    rows: list[dict[str, str]],
    metadata_rows: list[dict[str, str]],
    out_dir: Path,
    core_feature_threshold: float = 0.95,
    rare_feature_threshold: float = 0.05,
) -> dict[str, str]:
    diversity_dir = out_dir / "diversity"
    diversity_dir.mkdir(parents=True, exist_ok=True)
    samples = sorted({metadata_accession(row) for row in metadata_rows if metadata_accession(row)} | {row.get("assembly_accession", "") for row in rows if row.get("assembly_accession", "")})
    sample_features = _feature_sets_by_sample(rows)
    databases = sorted({row.get("database", "") for row in rows if row.get("database")})
    normalized = {row["assembly_accession"]: row for row in normalize_metadata_rows(metadata_rows) if row.get("assembly_accession")}

    richness_rows = []
    for sample in samples:
        features = sample_features.get(sample, set())
        counts = Counter(database for database, _feature in features)
        out = {
            "assembly_accession": sample,
            "total_feature_richness": str(len(features)),
            "database_count": str(len(counts)),
            "shannon_database_diversity": f"{_shannon(list(counts.values())):.6f}",
            "simpson_database_diversity": f"{_simpson(list(counts.values())):.6f}",
        }
        for metadata_column in METADATA_ALIASES:
            out[metadata_column] = normalized.get(sample, {}).get(metadata_column, "")
        richness_rows.append(out)
    richness_path = write_rows(
        diversity_dir / "feature_richness_by_sample.tsv",
        richness_rows,
        ["assembly_accession", "total_feature_richness", "database_count", "shannon_database_diversity", "simpson_database_diversity", *METADATA_ALIASES.keys()],
    )

    database_rows = []
    for sample in samples:
        features = sample_features.get(sample, set())
        for database in databases:
            count = len({feature for feature_db, feature in features if feature_db == database})
            database_rows.append({"assembly_accession": sample, "database": database, "feature_richness": str(count)})
    database_diversity_path = write_rows(
        diversity_dir / "database_diversity_by_sample.tsv",
        database_rows,
        ["assembly_accession", "database", "feature_richness"],
    )

    all_features = sorted({feature for features in sample_features.values() for feature in features})
    sample_count = len(samples)
    class_rows = []
    for database, feature_id in all_features:
        present = sum(1 for sample in samples if (database, feature_id) in sample_features.get(sample, set()))
        prevalence = present / sample_count if sample_count else 0.0
        if prevalence >= core_feature_threshold:
            feature_class = "core"
        elif prevalence < rare_feature_threshold:
            feature_class = "rare"
        else:
            feature_class = "accessory"
        class_rows.append({
            "database": database,
            "feature_id": feature_id,
            "present_count": str(present),
            "sample_count": str(sample_count),
            "prevalence": f"{prevalence:.6f}",
            "feature_class": feature_class,
            "core_threshold": f"{core_feature_threshold:.4f}",
            "rare_threshold": f"{rare_feature_threshold:.4f}",
        })
    core_accessory_path = write_rows(
        diversity_dir / "core_accessory_rare_features.tsv",
        class_rows,
        ["database", "feature_id", "present_count", "sample_count", "prevalence", "feature_class", "core_threshold", "rare_threshold"],
    )

    jaccard_rows = []
    for sample_a in samples:
        features_a = sample_features.get(sample_a, set())
        for sample_b in samples:
            features_b = sample_features.get(sample_b, set())
            union = features_a | features_b
            intersection = features_a & features_b
            jaccard = len(intersection) / len(union) if union else 1.0
            jaccard_rows.append({
                "sample_a": sample_a,
                "sample_b": sample_b,
                "shared_features": str(len(intersection)),
                "union_features": str(len(union)),
                "jaccard_similarity": f"{jaccard:.6f}",
                "jaccard_distance": f"{1 - jaccard:.6f}",
            })
    jaccard_path = write_rows(
        diversity_dir / "jaccard_distance_matrix.tsv",
        jaccard_rows,
        ["sample_a", "sample_b", "shared_features", "union_features", "jaccard_similarity", "jaccard_distance"],
    )

    accumulation_rows = []
    seen: set[tuple[str, str]] = set()
    for index, sample in enumerate(samples, start=1):
        seen.update(sample_features.get(sample, set()))
        accumulation_rows.append({
            "sample_order": str(index),
            "assembly_accession": sample,
            "cumulative_unique_features": str(len(seen)),
            "new_features_added": str(len(seen) - int(accumulation_rows[-1]["cumulative_unique_features"]) if accumulation_rows else len(seen)),
        })
    accumulation_path = write_rows(
        diversity_dir / "pan_feature_accumulation.tsv",
        accumulation_rows,
        ["sample_order", "assembly_accession", "cumulative_unique_features", "new_features_added"],
    )
    return {
        "feature_richness_by_sample": richness_path,
        "database_diversity_by_sample": database_diversity_path,
        "jaccard_distance_matrix": jaccard_path,
        "core_accessory_rare_features": core_accessory_path,
        "pan_feature_accumulation": accumulation_path,
    }


def write_statistical_summary(out_dir: Path) -> dict[str, str]:
    analysis_dir = out_dir / "metadata_feature_analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    feature_associations = read_table(analysis_dir / "feature_metadata_associations.tsv")
    burden_associations = read_table(analysis_dir / "database_burden_metadata_associations.tsv")
    category_associations = read_table(analysis_dir / "category_metadata_associations.tsv")
    top_findings = read_table(analysis_dir / "top_findings.tsv")
    metadata_columns = {
        row.get("metadata_column", "")
        for table in [feature_associations, burden_associations, category_associations, top_findings]
        for row in table
        if row.get("metadata_column")
    }
    warning_text = ";".join(
        ";".join([row.get("warning", ""), row.get("warning_flags", ""), row.get("lineage_warning_flags", "")])
        for table in [feature_associations, burden_associations, category_associations, top_findings]
        for row in table
    )
    q_values = [
        _float_or_none(row.get("q_value", ""))
        for row in feature_associations + top_findings
        if _float_or_none(row.get("q_value", "")) is not None
    ]
    rows = [
        {"metric": "features_tested", "value": str(len({(row.get("database", ""), row.get("feature_id", "")) for row in feature_associations if row.get("feature_id")})), "message": "Unique database-feature pairs in metadata-feature screening."},
        {"metric": "metadata_columns_tested", "value": str(len(metadata_columns)), "message": "Metadata columns represented in association summaries."},
        {"metric": "tests_performed", "value": str(len(feature_associations) + len(burden_associations) + len(category_associations)), "message": "Screening tests/summaries written across feature, burden, and category association tables."},
        {"metric": "top_findings_generated", "value": str(len(top_findings)), "message": "Top exploratory findings emitted for report navigation."},
        {"metric": "q_le_0_05_findings", "value": str(sum(1 for value in q_values if value is not None and value <= 0.05)), "message": "Findings with q<=0.05 when q-values are available."},
        {"metric": "small_group_warnings", "value": str(warning_text.count("small_group") + warning_text.count("low_sample_warning")), "message": "Small-group or low-sample warning flags."},
        {"metric": "bioproject_warnings", "value": str(warning_text.count("single_bioproject_dominance") + warning_text.count("bioproject_bias_warning")), "message": "BioProject/study-dominance warning flags."},
        {"metric": "lineage_warnings", "value": str(warning_text.count("single_ST_dominance") + warning_text.count("single_ani_cluster_dominance") + warning_text.count("metadata_lineage_confounding")), "message": "Lineage/ST/ANI warning flags."},
        {"metric": "insufficient_lineage_warnings", "value": str(warning_text.count("insufficient_lineage_data")), "message": "Findings or summaries where lineage context was unavailable."},
    ]
    path = write_rows(analysis_dir / "statistical_summary.tsv", rows, ["metric", "value", "message"])
    return {"statistical_summary": path}


def write_feature_eligibility_and_prevalence(
    rows: list[dict[str, str]],
    metadata_rows: list[dict[str, str]],
    out_dir: Path,
    top_n_features_per_database: int = 25,
) -> dict[str, str]:
    analysis_dir = out_dir / "metadata_feature_analysis"
    prevalence_dir = analysis_dir / "prevalence_tables"
    prevalence_dir.mkdir(parents=True, exist_ok=True)
    samples = sorted({metadata_accession(row) for row in metadata_rows if metadata_accession(row)} | {row.get("assembly_accession", "") for row in rows if row.get("assembly_accession", "")})
    sample_count = len(samples)
    presence = feature_presence(rows)
    category_by_feature = {}
    for row in rows:
        category_by_feature[(row.get("database", ""), row.get("feature_id", ""))] = row.get("feature_category", "")

    eligibility_rows = []
    for (database, feature_id), present_samples in sorted(presence.items()):
        present = len(present_samples)
        absent = max(sample_count - present, 0)
        eligible = present >= 2 and absent >= 2
        eligibility_rows.append({
            "database": database,
            "feature_id": feature_id,
            "feature_category": category_by_feature.get((database, feature_id), ""),
            "present_count": str(present),
            "absent_count": str(absent),
            "prevalence": f"{(present / sample_count if sample_count else 0):.4f}",
            "eligible": str(eligible).lower(),
            "reason": "eligible" if eligible else "too_rare_or_too_common",
        })
    feature_eligibility_path = write_rows(
        analysis_dir / "feature_eligibility.tsv",
        eligibility_rows,
        ["database", "feature_id", "feature_category", "present_count", "absent_count", "prevalence", "eligible", "reason"],
    )
    top_feature_rows = []
    by_database: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in eligibility_rows:
        by_database[row["database"]].append(row)
    for database, database_rows in sorted(by_database.items()):
        ranked = sorted(
            database_rows,
            key=lambda row: (-_float_or_none(row.get("prevalence", "")) if _float_or_none(row.get("prevalence", "")) is not None else 0, row["feature_id"]),
        )[:max(top_n_features_per_database, 0)]
        for rank, row in enumerate(ranked, start=1):
            top_feature_rows.append({
                "database": database,
                "rank": str(rank),
                "feature_id": row["feature_id"],
                "feature_category": row["feature_category"],
                "present_count": row["present_count"],
                "absent_count": row["absent_count"],
                "prevalence": row["prevalence"],
            })
    top_features_path = write_rows(
        analysis_dir / "top_features_by_database.tsv",
        top_feature_rows,
        ["database", "rank", "feature_id", "feature_category", "present_count", "absent_count", "prevalence"],
    )

    normalized = normalize_metadata_rows(metadata_rows)
    metadata_by_accession = {row["assembly_accession"]: row for row in normalized if row.get("assembly_accession")}
    prevalence_outputs = {}
    feature_metadata_association_rows = []
    for metadata_column in METADATA_ALIASES:
        group_values = [row.get(metadata_column, "") for row in normalized if row.get(metadata_column, "")]
        if len(set(group_values)) < 2:
            continue
        rows_out = []
        groups = sorted(set(group_values))
        group_samples = {
            group: {row["assembly_accession"] for row in normalized if row.get(metadata_column) == group and row.get("assembly_accession")}
            for group in groups
        }
        for (database, feature_id), present_samples in sorted(presence.items()):
            group_prevalence = []
            for group, members in group_samples.items():
                if not members:
                    continue
                present = len(present_samples & members)
                prevalence = present / len(members)
                group_prevalence.append((group, len(members), present, prevalence))
                rows_out.append({
                    "metadata_column": metadata_column,
                    "metadata_value": group,
                    "database": database,
                    "feature_id": feature_id,
                    "feature_category": category_by_feature.get((database, feature_id), ""),
                    "n_group": str(len(members)),
                    "present_count": str(present),
                    "prevalence": f"{prevalence:.4f}",
                })
            eligible_group_prevalence = [item for item in group_prevalence if item[1] >= 3]
            if len(eligible_group_prevalence) >= 2:
                top_group = max(eligible_group_prevalence, key=lambda item: item[3])
                low_group = min(eligible_group_prevalence, key=lambda item: item[3])
                feature_metadata_association_rows.append({
                    "database": database,
                    "feature_id": feature_id,
                    "feature_category": category_by_feature.get((database, feature_id), ""),
                    "metadata_column": metadata_column,
                    "metadata_type": "categorical_or_alias",
                    "test_used": "prevalence_range_screen",
                    "n_total": str(sum(item[1] for item in eligible_group_prevalence)),
                    "n_present": str(sum(item[2] for item in eligible_group_prevalence)),
                    "n_absent": str(sum(item[1] - item[2] for item in eligible_group_prevalence)),
                    "groups_tested": str(len(eligible_group_prevalence)),
                    "effect_size": f"{top_group[3] - low_group[3]:.4f}",
                    "odds_ratio": "",
                    "p_value": "",
                    "q_value": "",
                    "top_enriched_group": top_group[0],
                    "top_enriched_group_prevalence": f"{top_group[3]:.4f}",
                    "lowest_group": low_group[0],
                    "lowest_group_prevalence": f"{low_group[3]:.4f}",
                    "status": "screening_summary",
                    "warning": "small_group" if min(item[1] for item in eligible_group_prevalence) < 5 else "",
                })
        if rows_out:
            path = prevalence_dir / f"all_databases__by__{metadata_column}.tsv"
            prevalence_outputs[f"prevalence_by_{metadata_column}"] = write_rows(
                path,
                rows_out,
                ["metadata_column", "metadata_value", "database", "feature_id", "feature_category", "n_group", "present_count", "prevalence"],
            )
    feature_metadata_association_rows = sorted(
        feature_metadata_association_rows,
        key=lambda row: (-float(row["effect_size"]), row["metadata_column"], row["database"], row["feature_id"]),
    )
    feature_metadata_association_path = write_rows(
        analysis_dir / "feature_metadata_associations.tsv",
        feature_metadata_association_rows,
        [
            "database",
            "feature_id",
            "feature_category",
            "metadata_column",
            "metadata_type",
            "test_used",
            "n_total",
            "n_present",
            "n_absent",
            "groups_tested",
            "effect_size",
            "odds_ratio",
            "p_value",
            "q_value",
            "top_enriched_group",
            "top_enriched_group_prevalence",
            "lowest_group",
            "lowest_group_prevalence",
            "status",
            "warning",
        ],
    )

    burden_rows = []
    features_by_sample_database: dict[tuple[str, str], set[str]] = defaultdict(set)
    categories_by_sample_database: dict[tuple[str, str], set[str]] = defaultdict(set)
    features_by_sample_database_category: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    for row in rows:
        sample = row.get("assembly_accession", "") or row.get("sample_id", "")
        database = row.get("database", "")
        feature_id = row.get("feature_id", "")
        category = row.get("feature_category", "")
        if sample and database and feature_id:
            features_by_sample_database[(sample, database)].add(feature_id)
            if category:
                categories_by_sample_database[(sample, database)].add(category)
                features_by_sample_database_category[(sample, database, category)].add(feature_id)
    databases = sorted({row.get("database", "") for row in rows if row.get("database")})
    for sample in samples:
        metadata = metadata_by_accession.get(sample, {})
        for database in databases:
            key = (sample, database)
            out = {
                "assembly_accession": sample,
                "database": database,
                "unique_feature_count": str(len(features_by_sample_database.get(key, set()))),
                "category_count": str(len(categories_by_sample_database.get(key, set()))),
            }
            for metadata_column in METADATA_ALIASES:
                out[metadata_column] = metadata.get(metadata_column, "")
            burden_rows.append(out)
    burden_fields = ["assembly_accession", "database", "unique_feature_count", "category_count", *METADATA_ALIASES.keys()]
    burden_path = write_rows(analysis_dir / "database_burden_by_sample.tsv", burden_rows, burden_fields)

    burden_association_rows = []
    for metadata_column in METADATA_ALIASES:
        for database in databases:
            grouped: dict[str, list[int]] = defaultdict(list)
            for row in burden_rows:
                if row["database"] != database or not row.get(metadata_column):
                    continue
                grouped[row[metadata_column]].append(int(row["unique_feature_count"]))
            eligible = [(group, values) for group, values in grouped.items() if len(values) >= 3]
            if len(eligible) < 2:
                continue
            means = [(group, len(values), sum(values) / len(values)) for group, values in eligible]
            top = max(means, key=lambda item: item[2])
            low = min(means, key=lambda item: item[2])
            burden_association_rows.append({
                "database": database,
                "metadata_column": metadata_column,
                "test_used": "mean_burden_range_screen",
                "groups_tested": str(len(means)),
                "top_group": top[0],
                "top_group_n": str(top[1]),
                "top_group_mean_unique_features": f"{top[2]:.4f}",
                "lowest_group": low[0],
                "lowest_group_n": str(low[1]),
                "lowest_group_mean_unique_features": f"{low[2]:.4f}",
                "effect_size": f"{top[2] - low[2]:.4f}",
                "warning": "small_group" if min(item[1] for item in means) < 5 else "",
            })
    burden_association_rows = sorted(
        burden_association_rows,
        key=lambda row: (-float(row["effect_size"]), row["metadata_column"], row["database"]),
    )
    burden_association_path = write_rows(
        analysis_dir / "database_burden_metadata_associations.tsv",
        burden_association_rows,
        [
            "database",
            "metadata_column",
            "test_used",
            "groups_tested",
            "top_group",
            "top_group_n",
            "top_group_mean_unique_features",
            "lowest_group",
            "lowest_group_n",
            "lowest_group_mean_unique_features",
            "effect_size",
            "warning",
        ],
    )

    category_burden_rows = []
    for sample in samples:
        metadata = metadata_by_accession.get(sample, {})
        for sample_key, database, category in sorted(features_by_sample_database_category):
            if sample_key != sample:
                continue
            out = {
                "assembly_accession": sample,
                "database": database,
                "feature_category": category,
                "unique_feature_count": str(len(features_by_sample_database_category[(sample, database, category)])),
            }
            for metadata_column in METADATA_ALIASES:
                out[metadata_column] = metadata.get(metadata_column, "")
            category_burden_rows.append(out)
    category_burden_fields = ["assembly_accession", "database", "feature_category", "unique_feature_count", *METADATA_ALIASES.keys()]
    category_burden_path = write_rows(analysis_dir / "category_burden_by_sample.tsv", category_burden_rows, category_burden_fields)

    category_association_rows = []
    for metadata_column in METADATA_ALIASES:
        grouped_values: dict[tuple[str, str, str], list[int]] = defaultdict(list)
        for row in category_burden_rows:
            group = row.get(metadata_column, "")
            if not group:
                continue
            grouped_values[(row["database"], row["feature_category"], group)].append(int(row["unique_feature_count"]))
        by_category: dict[tuple[str, str], dict[str, list[int]]] = defaultdict(dict)
        for (database, category, group), values in grouped_values.items():
            by_category[(database, category)][group] = values
        for (database, category), group_map in by_category.items():
            eligible = [(group, values) for group, values in group_map.items() if len(values) >= 3]
            if len(eligible) < 2:
                continue
            means = [(group, len(values), sum(values) / len(values)) for group, values in eligible]
            top = max(means, key=lambda item: item[2])
            low = min(means, key=lambda item: item[2])
            category_association_rows.append({
                "database": database,
                "feature_category": category,
                "metadata_column": metadata_column,
                "test_used": "mean_category_burden_range_screen",
                "groups_tested": str(len(means)),
                "top_group": top[0],
                "top_group_mean_unique_features": f"{top[2]:.4f}",
                "lowest_group": low[0],
                "lowest_group_mean_unique_features": f"{low[2]:.4f}",
                "effect_size": f"{top[2] - low[2]:.4f}",
                "warning": "small_group" if min(item[1] for item in means) < 5 else "",
            })
    category_association_rows = sorted(
        category_association_rows,
        key=lambda row: (-float(row["effect_size"]), row["metadata_column"], row["database"], row["feature_category"]),
    )
    category_association_path = write_rows(
        analysis_dir / "category_metadata_associations.tsv",
        category_association_rows,
        [
            "database",
            "feature_category",
            "metadata_column",
            "test_used",
            "groups_tested",
            "top_group",
            "top_group_mean_unique_features",
            "lowest_group",
            "lowest_group_mean_unique_features",
            "effect_size",
            "warning",
        ],
    )

    top_findings = []
    for row in feature_metadata_association_rows:
        if float(row["effect_size"]) >= 0.5:
            top_findings.append({
                "finding_type": "feature_metadata_association",
                "database": row["database"],
                "feature_id": row["feature_id"],
                "metadata_column": row["metadata_column"],
                "metadata_value": row["top_enriched_group"],
                "effect_size": row["effect_size"],
                "p_value": row.get("p_value", ""),
                "q_value": row.get("q_value", ""),
                "warning": row.get("warning", ""),
                "message": f"{row['database']} feature {row['feature_id']} was enriched in {row['metadata_column']}={row['top_enriched_group']} compared with {row['lowest_group']} (prevalence difference {row['effect_size']}).",
            })
    for row in burden_association_rows:
        if float(row["effect_size"]) > 0:
            top_findings.append({
                "finding_type": "database_burden_metadata_association",
                "database": row["database"],
                "feature_id": "database_burden",
                "metadata_column": row["metadata_column"],
                "metadata_value": row["top_group"],
                "effect_size": row["effect_size"],
                "p_value": "",
                "q_value": "",
                "warning": row.get("warning", ""),
                "message": f"{row['database']} feature burden was highest in {row['metadata_column']}={row['top_group']} compared with {row['lowest_group']} (mean difference {row['effect_size']}).",
            })
    for key, path in prevalence_outputs.items():
        for row in read_table(Path(path)):
            n_group = int(row.get("n_group", "0") or 0)
            if n_group >= 3 and float(row.get("prevalence", "0") or 0) >= 0.8:
                top_findings.append({
                    "finding_type": "metadata_feature_prevalence",
                    "database": row["database"],
                    "feature_id": row["feature_id"],
                    "metadata_column": row["metadata_column"],
                    "metadata_value": row["metadata_value"],
                    "effect_size": row["prevalence"],
                    "p_value": "",
                    "q_value": "",
                    "warning": "small_group" if n_group < 5 else "",
                    "message": f"{row['database']} feature {row['feature_id']} was common in {row['metadata_column']}={row['metadata_value']} (prevalence {row['prevalence']}).",
                })
    top_findings = sorted(top_findings, key=lambda row: float(row["effect_size"]), reverse=True)[:50]
    st_by_sample = _st_by_sample(rows)
    ani_cluster_by_sample = _ani_cluster_by_sample(out_dir)
    top_findings = _annotate_top_findings(
        top_findings,
        normalized,
        presence,
        _bioproject_by_sample(metadata_rows),
        st_by_sample=st_by_sample,
        ani_cluster_by_sample=ani_cluster_by_sample,
    )
    top_findings_path = write_rows(
        analysis_dir / "top_findings.tsv",
        top_findings,
        [
            "finding_type",
            "database",
            "feature_id",
            "metadata_column",
            "metadata_value",
            "effect_size",
            "p_value",
            "q_value",
            "supporting_samples",
            "largest_bioproject",
            "largest_bioproject_fraction",
            "dominant_ST",
            "dominant_ST_fraction",
            "dominant_ani_cluster",
            "dominant_ani_cluster_fraction",
            "lineage_warning_flags",
            "warning_flags",
            "interpretation_label",
            "message",
        ],
    )
    top_findings_md = analysis_dir / "top_findings.md"
    with top_findings_md.open("w", encoding="utf-8") as handle:
        handle.write("# Top Metadata-Feature Findings\n\n")
        handle.write("These are screening summaries only. They show sample-level enrichment and do not prove causality, physical linkage, or transfer.\n\n")
        for row in top_findings[:20]:
            warning_note = f" [{row['interpretation_label']}]" if row.get("interpretation_label") else ""
            handle.write(f"- {row['message']}{warning_note}\n")
    bioproject_bias_path = write_bioproject_bias_report(rows, metadata_rows, top_findings, out_dir)
    lineage_outputs = write_lineage_analysis(rows, metadata_rows, top_findings, out_dir)
    statistical_outputs = write_statistical_summary(out_dir)

    return {
        "feature_eligibility": feature_eligibility_path,
        "top_features_by_database": top_features_path,
        "feature_metadata_associations": feature_metadata_association_path,
        "database_burden_by_sample": burden_path,
        "database_burden_metadata_associations": burden_association_path,
        "category_burden_by_sample": category_burden_path,
        "category_metadata_associations": category_association_path,
        "top_findings": top_findings_path,
        "top_findings_md": str(top_findings_md),
        "bioproject_bias_report": bioproject_bias_path,
        **lineage_outputs,
        **statistical_outputs,
        **prevalence_outputs,
    }


def _fisher_like_phi(a: int, b: int, c: int, d: int) -> float:
    denominator = math.sqrt((a + b) * (c + d) * (a + c) * (b + d))
    return ((a * d) - (b * c)) / denominator if denominator else 0.0


def _log_comb(n: int, k: int) -> float:
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _hypergeom_probability(x: int, row_total: int, col_total: int, n_total: int) -> float:
    log_p = (
        _log_comb(col_total, x)
        + _log_comb(n_total - col_total, row_total - x)
        - _log_comb(n_total, row_total)
    )
    return math.exp(log_p) if math.isfinite(log_p) else 0.0


def fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> float:
    n_total = a + b + c + d
    if n_total == 0:
        return 1.0
    row_total = a + b
    col_total = a + c
    min_x = max(0, row_total - (n_total - col_total))
    max_x = min(row_total, col_total)
    observed = _hypergeom_probability(a, row_total, col_total, n_total)
    p_value = 0.0
    for x in range(min_x, max_x + 1):
        probability = _hypergeom_probability(x, row_total, col_total, n_total)
        if probability <= observed + 1e-12:
            p_value += probability
    return min(max(p_value, 0.0), 1.0)


def add_bh_qvalues(rows: list[dict[str, str]], p_column: str = "p_value", q_column: str = "q_value") -> None:
    indexed = []
    for index, row in enumerate(rows):
        p_value = _as_float(row.get(p_column, ""))
        if p_value is not None:
            indexed.append((index, min(max(p_value, 0.0), 1.0)))
    if not indexed:
        return
    indexed.sort(key=lambda item: item[1], reverse=True)
    total = len(indexed)
    running_min = 1.0
    qvalues: dict[int, float] = {}
    for rank_from_largest, (index, p_value) in enumerate(indexed, start=1):
        rank = total - rank_from_largest + 1
        running_min = min(running_min, p_value * total / rank)
        qvalues[index] = min(running_min, 1.0)
    for index, q_value in qvalues.items():
        rows[index][q_column] = f"{q_value:.6g}"


AMR_CONTEXT_DATABASES = {"amr", "amrfinderplus"}
MGE_CONTEXT_DATABASES = {"isfinder", "mobileelementfinder", "integronfinder"}
PLASMID_CONTEXT_DATABASES = {"plasmidfinder", "mobsuite"}
INTEGRON_CONTEXT_DATABASES = {"integronfinder"}


def _feature_interval(row: dict[str, str]) -> tuple[int, int] | None:
    start = _as_int(row.get("start", ""))
    end = _as_int(row.get("end", ""))
    if start is None or end is None:
        return None
    if start > end:
        start, end = end, start
    return start, end


def _feature_distance(row_a: dict[str, str], row_b: dict[str, str]) -> tuple[str, str]:
    interval_a = _feature_interval(row_a)
    interval_b = _feature_interval(row_b)
    if interval_a is None or interval_b is None:
        return "", "level_2_same_contig_coordinates_missing"
    a_start, a_end = interval_a
    b_start, b_end = interval_b
    if max(a_start, b_start) <= min(a_end, b_end):
        return "0", "level_4_same_contig_overlapping"
    distance = min(abs(b_start - a_end), abs(a_start - b_end))
    if distance <= 10000:
        return str(distance), "level_3_same_contig_within_10kb"
    return str(distance), "level_2_same_contig"


def _context_feature_rows(rows: list[dict[str, str]], databases: set[str]) -> list[dict[str, str]]:
    return [
        row
        for row in rows
        if row.get("presence", "1") == "1"
        and row.get("database") in databases
        and row.get("assembly_accession")
        and row.get("contig")
        and row.get("feature_id")
    ]


def _same_contig_pairs(
    rows: list[dict[str, str]],
    left_databases: set[str],
    right_databases: set[str],
    context: str,
) -> list[dict[str, str]]:
    left_rows = _context_feature_rows(rows, left_databases)
    right_rows = _context_feature_rows(rows, right_databases)
    by_key: dict[tuple[str, str], tuple[list[dict[str, str]], list[dict[str, str]]]] = defaultdict(lambda: ([], []))
    for row in left_rows:
        by_key[(row["assembly_accession"], row["contig"])][0].append(row)
    for row in right_rows:
        by_key[(row["assembly_accession"], row["contig"])][1].append(row)

    pair_rows = []
    seen: set[tuple[str, str, str, str, str, str, str]] = set()
    for (assembly, contig), (left, right) in sorted(by_key.items()):
        for row_a in left:
            for row_b in right:
                if row_a is row_b:
                    continue
                if (
                    row_a.get("database") == row_b.get("database")
                    and row_a.get("feature_id") == row_b.get("feature_id")
                    and row_a.get("start") == row_b.get("start")
                    and row_a.get("end") == row_b.get("end")
                ):
                    continue
                distance, level = _feature_distance(row_a, row_b)
                key = (
                    assembly,
                    contig,
                    row_a.get("database", ""),
                    row_a.get("feature_id", ""),
                    row_b.get("database", ""),
                    row_b.get("feature_id", ""),
                    distance,
                )
                if key in seen:
                    continue
                seen.add(key)
                pair_rows.append({
                    "assembly_accession": assembly,
                    "contig": contig,
                    "context": context,
                    "feature_a_database": row_a.get("database", ""),
                    "feature_a_id": row_a.get("feature_id", ""),
                    "feature_a_category": row_a.get("feature_category", ""),
                    "feature_a_start": row_a.get("start", ""),
                    "feature_a_end": row_a.get("end", ""),
                    "feature_b_database": row_b.get("database", ""),
                    "feature_b_id": row_b.get("feature_id", ""),
                    "feature_b_category": row_b.get("feature_category", ""),
                    "feature_b_start": row_b.get("start", ""),
                    "feature_b_end": row_b.get("end", ""),
                    "distance_bp": distance,
                    "interpretation_level": level,
                    "evidence_level": level,
                    "interpretation_warning": "Same-contig/proximity evidence is stronger than same-genome co-occurrence but does not prove transfer, expression, phenotype, or plasmid localization.",
                })
    return pair_rows


def _normalize_amr_symbol(value: str) -> str:
    text = str(value or "").strip().lower()
    text = re.sub(r"[^a-z0-9]+", "", text)
    return text


def _category_terms(row: dict[str, str]) -> set[str]:
    raw_values = [
        row.get("feature_category", ""),
        row.get("drug_class", ""),
        row.get("feature_subcategory", ""),
        row.get("raw_category", ""),
    ]
    terms = set()
    for value in raw_values:
        for part in re.split(r"[;,/|]", str(value or "").lower()):
            clean = re.sub(r"\s+", "_", part.strip())
            if clean and clean not in {"na", "n/a", "none", "unknown", "-"}:
                terms.add(clean)
    return terms


def _amr_calls_by_normalized_symbol(rows: list[dict[str, str]], database: str) -> dict[str, list[dict[str, str]]]:
    calls: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        if row.get("database") != database or row.get("presence", "1") != "1":
            continue
        feature_id = row.get("feature_id", "")
        normalized = _normalize_amr_symbol(feature_id)
        sample = row.get("assembly_accession", "") or row.get("sample_id", "")
        if feature_id and normalized and sample:
            calls[normalized].append(row)
    return calls


def write_cross_database_outputs(rows: list[dict[str, str]], metadata_rows: list[dict[str, str]], out_dir: Path, max_features: int = 300) -> dict[str, str]:
    cross_dir = out_dir / "cross_database"
    cross_dir.mkdir(parents=True, exist_ok=True)
    samples = sorted({metadata_accession(row) for row in metadata_rows if metadata_accession(row)} | {row.get("assembly_accession", "") for row in rows if row.get("assembly_accession", "")})
    sample_count = len(samples)
    presence = feature_presence(rows)
    ordered_features = sorted(presence, key=lambda item: (-len(presence[item]), item[0], item[1]))[:max_features]
    category_by_feature = {(row.get("database", ""), row.get("feature_id", "")): row.get("feature_category", "") for row in rows}

    cooccurrence_rows = []
    for feature_a, feature_b in itertools.combinations(ordered_features, 2):
        a_samples = presence[feature_a]
        b_samples = presence[feature_b]
        n_both = len(a_samples & b_samples)
        n_a_only = len(a_samples - b_samples)
        n_b_only = len(b_samples - a_samples)
        n_neither = max(sample_count - n_both - n_a_only - n_b_only, 0)
        union = len(a_samples | b_samples)
        odds_denominator = n_a_only * n_b_only
        odds_ratio = ((n_both * n_neither) / odds_denominator) if odds_denominator else ""
        p_value = fisher_exact_two_sided(n_both, n_a_only, n_b_only, n_neither)
        cooccurrence_rows.append({
            "feature_a_database": feature_a[0],
            "feature_a_id": feature_a[1],
            "feature_a_category": category_by_feature.get(feature_a, ""),
            "feature_b_database": feature_b[0],
            "feature_b_id": feature_b[1],
            "feature_b_category": category_by_feature.get(feature_b, ""),
            "n_total": str(sample_count),
            "n_both_present": str(n_both),
            "n_a_only": str(n_a_only),
            "n_b_only": str(n_b_only),
            "n_neither": str(n_neither),
            "jaccard": f"{(n_both / union if union else 0):.4f}",
            "phi": f"{_fisher_like_phi(n_both, n_a_only, n_b_only, n_neither):.4f}",
            "odds_ratio": f"{odds_ratio:.4f}" if isinstance(odds_ratio, float) else "",
            "p_value": f"{p_value:.6g}",
            "q_value": "",
            "status": "screening_exact_fisher_bh_fdr",
        })
    add_bh_qvalues(cooccurrence_rows)
    cooccurrence_path = write_rows(
        cross_dir / "feature_cooccurrence.tsv",
        cooccurrence_rows,
        [
            "feature_a_database",
            "feature_a_id",
            "feature_a_category",
            "feature_b_database",
            "feature_b_id",
            "feature_b_category",
            "n_total",
            "n_both_present",
            "n_a_only",
            "n_b_only",
            "n_neither",
            "jaccard",
            "phi",
            "odds_ratio",
            "p_value",
            "q_value",
            "status",
        ],
    )

    samples_by_database: dict[str, set[str]] = defaultdict(set)
    for (database, _feature), members in presence.items():
        samples_by_database[database].update(members)
    database_rows = []
    for database_a, database_b in itertools.combinations(sorted(samples_by_database), 2):
        a = samples_by_database[database_a]
        b = samples_by_database[database_b]
        database_rows.append({
            "database_a": database_a,
            "database_b": database_b,
            "n_total": str(sample_count),
            "n_both_present": str(len(a & b)),
            "n_a_only": str(len(a - b)),
            "n_b_only": str(len(b - a)),
            "n_neither": str(max(sample_count - len(a | b), 0)),
            "jaccard": f"{(len(a & b) / len(a | b) if a or b else 0):.4f}",
        })
    database_summary_path = write_rows(
        cross_dir / "database_cooccurrence_summary.tsv",
        database_rows,
        ["database_a", "database_b", "n_total", "n_both_present", "n_a_only", "n_b_only", "n_neither", "jaccard"],
    )

    proximity_fields = [
        "assembly_accession",
        "contig",
        "context",
        "feature_a_database",
        "feature_a_id",
        "feature_a_category",
        "feature_a_start",
        "feature_a_end",
        "feature_b_database",
        "feature_b_id",
        "feature_b_category",
        "feature_b_start",
        "feature_b_end",
        "distance_bp",
        "interpretation_level",
        "evidence_level",
        "interpretation_warning",
    ]
    amr_mge_same_contig_rows = _same_contig_pairs(rows, AMR_CONTEXT_DATABASES, MGE_CONTEXT_DATABASES, "amr_mge")
    amr_plasmid_same_contig_rows = _same_contig_pairs(rows, AMR_CONTEXT_DATABASES, PLASMID_CONTEXT_DATABASES, "amr_plasmid")
    amr_integron_same_contig_rows = _same_contig_pairs(rows, AMR_CONTEXT_DATABASES, INTEGRON_CONTEXT_DATABASES, "amr_integron")
    feature_proximity_all_rows = sorted(
        amr_mge_same_contig_rows + amr_plasmid_same_contig_rows + amr_integron_same_contig_rows,
        key=lambda row: (
            row["assembly_accession"],
            row["contig"],
            row["context"],
            row["feature_a_database"],
            row["feature_a_id"],
            row["feature_b_database"],
            row["feature_b_id"],
            _as_int(row["distance_bp"]) if row["distance_bp"] else 10**18,
        ),
    )
    proximity_feature_counts: Counter[tuple[str, str]] = Counter()
    for row in feature_proximity_all_rows:
        proximity_feature_counts[(row["feature_a_database"], row["feature_a_id"])] += 1
        proximity_feature_counts[(row["feature_b_database"], row["feature_b_id"])] += 1
    proximity_allowed_features = {
        feature
        for feature, _count in sorted(
            proximity_feature_counts.items(),
            key=lambda item: (-item[1], item[0][0], item[0][1]),
        )[:max_features]
    }
    if proximity_allowed_features and len(proximity_feature_counts) > max_features:
        feature_proximity_rows = [
            row for row in feature_proximity_all_rows
            if (row["feature_a_database"], row["feature_a_id"]) in proximity_allowed_features
            and (row["feature_b_database"], row["feature_b_id"]) in proximity_allowed_features
        ]
    else:
        feature_proximity_rows = feature_proximity_all_rows
    amr_mge_same_contig_path = write_rows(cross_dir / "amr_mge_same_contig.tsv", amr_mge_same_contig_rows, proximity_fields)
    amr_plasmid_same_contig_path = write_rows(cross_dir / "amr_plasmid_same_contig.tsv", amr_plasmid_same_contig_rows, proximity_fields)
    amr_integron_same_contig_path = write_rows(cross_dir / "amr_integron_same_contig.tsv", amr_integron_same_contig_rows, proximity_fields)
    feature_proximity_all_path = write_rows(cross_dir / "feature_proximity_all.tsv", feature_proximity_all_rows, proximity_fields)
    feature_proximity_path = write_rows(cross_dir / "feature_proximity.tsv", feature_proximity_rows, proximity_fields)
    samples_with_same_contig = {
        row["assembly_accession"]
        for row in feature_proximity_all_rows
        if row.get("assembly_accession")
    }

    by_sample_database: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in rows:
        sample = row.get("assembly_accession", "") or row.get("sample_id", "")
        if sample and row.get("database") and row.get("feature_id"):
            by_sample_database[(sample, row["database"])].add(row["feature_id"])
    context_rows = []
    for sample in samples:
        amr = by_sample_database.get((sample, "amr"), set()) | by_sample_database.get((sample, "amrfinderplus"), set())
        plasmid = by_sample_database.get((sample, "plasmidfinder"), set()) | by_sample_database.get((sample, "mobsuite"), set())
        mge = (
            by_sample_database.get((sample, "isfinder"), set())
            | by_sample_database.get((sample, "mobileelementfinder"), set())
            | by_sample_database.get((sample, "integronfinder"), set())
        )
        context_rows.append({
            "assembly_accession": sample,
            "amr_features": ";".join(sorted(amr)),
            "plasmid_features": ";".join(sorted(plasmid)),
            "mge_features": ";".join(sorted(mge)),
            "amr_count": str(len(amr)),
            "plasmid_count": str(len(plasmid)),
            "mge_count": str(len(mge)),
            "same_contig_evidence": "yes" if sample in samples_with_same_contig else "no",
        })
    amr_mge_path = write_rows(
        cross_dir / "amr_mge_context.tsv",
        context_rows,
        ["assembly_accession", "amr_features", "plasmid_features", "mge_features", "amr_count", "plasmid_count", "mge_count", "same_contig_evidence"],
    )
    amr_plasmid_path = write_rows(
        cross_dir / "amr_plasmid_context.tsv",
        context_rows,
        ["assembly_accession", "amr_features", "plasmid_features", "mge_features", "amr_count", "plasmid_count", "mge_count", "same_contig_evidence"],
    )

    abricate_calls = _amr_calls_by_normalized_symbol(rows, "amr")
    amrfinder_calls = _amr_calls_by_normalized_symbol(rows, "amrfinderplus")
    concordance_rows = []
    for normalized_feature in sorted(set(abricate_calls) | set(amrfinder_calls)):
        abricate_rows = abricate_calls.get(normalized_feature, [])
        amrfinder_rows = amrfinder_calls.get(normalized_feature, [])
        abricate_features = sorted({row.get("feature_id", "") for row in abricate_rows if row.get("feature_id")})
        amrfinder_features = sorted({row.get("feature_id", "") for row in amrfinder_rows if row.get("feature_id")})
        abricate_samples = {row.get("assembly_accession", "") or row.get("sample_id", "") for row in abricate_rows}
        amrfinder_samples = {row.get("assembly_accession", "") or row.get("sample_id", "") for row in amrfinder_rows}
        shared_samples = abricate_samples & amrfinder_samples
        if shared_samples:
            status = "called_by_both"
        elif abricate_rows and amrfinder_rows:
            status = "same_symbol_different_samples"
        elif abricate_rows:
            status = "abricate_only"
        else:
            status = "amrfinderplus_only"
        concordance_rows.append({
            "feature_id": ";".join(abricate_features or amrfinder_features),
            "sample_id": "",
            "normalized_feature_id": normalized_feature,
            "abricate_feature_ids": ";".join(abricate_features),
            "amrfinderplus_feature_ids": ";".join(amrfinder_features),
            "abricate_present": str(bool(abricate_rows)).lower(),
            "amrfinderplus_present": str(bool(amrfinder_rows)).lower(),
            "abricate_sample_count": str(len(abricate_samples)),
            "amrfinderplus_sample_count": str(len(amrfinder_samples)),
            "shared_sample_count": str(len(shared_samples)),
            "status": status,
            "possible_match_basis": "normalized_gene_symbol",
            "interpretation_note": "Gene-symbol matching is normalized for punctuation/case; inspect raw tool outputs for naming differences.",
        })

    abricate_by_sample_class: dict[tuple[str, str], set[str]] = defaultdict(set)
    amrfinder_by_sample_class: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in rows:
        database = row.get("database", "")
        if database not in {"amr", "amrfinderplus"} or row.get("presence", "1") != "1":
            continue
        sample = row.get("assembly_accession", "") or row.get("sample_id", "")
        if not sample:
            continue
        for term in _category_terms(row):
            if database == "amr":
                abricate_by_sample_class[(sample, term)].add(row.get("feature_id", ""))
            else:
                amrfinder_by_sample_class[(sample, term)].add(row.get("feature_id", ""))
    for sample_class in sorted(set(abricate_by_sample_class) & set(amrfinder_by_sample_class)):
        sample, category = sample_class
        abricate_features = sorted(abricate_by_sample_class[sample_class])
        amrfinder_features = sorted(amrfinder_by_sample_class[sample_class])
        shared_normalized = {
            _normalize_amr_symbol(feature)
            for feature in abricate_features
        } & {
            _normalize_amr_symbol(feature)
            for feature in amrfinder_features
        }
        if shared_normalized:
            continue
        concordance_rows.append({
            "feature_id": f"class:{category}",
            "sample_id": sample,
            "normalized_feature_id": f"class:{category}",
            "abricate_feature_ids": ";".join(abricate_features),
            "amrfinderplus_feature_ids": ";".join(amrfinder_features),
            "abricate_present": "true",
            "amrfinderplus_present": "true",
            "abricate_sample_count": "1",
            "amrfinderplus_sample_count": "1",
            "shared_sample_count": "1",
            "status": "possible_class_match",
            "possible_match_basis": "same_sample_same_drug_or_feature_class",
            "interpretation_note": "Tools reported different feature names in the same sample but shared a drug/class label; treat as possible concordance, not exact agreement.",
        })
    concordance_path = write_rows(
        cross_dir / "amrfinder_abricate_concordance.tsv",
        concordance_rows,
        [
            "feature_id",
            "sample_id",
            "normalized_feature_id",
            "abricate_feature_ids",
            "amrfinderplus_feature_ids",
            "abricate_present",
            "amrfinderplus_present",
            "abricate_sample_count",
            "amrfinderplus_sample_count",
            "shared_sample_count",
            "status",
            "possible_match_basis",
            "interpretation_note",
        ],
    )
    return {
        "feature_cooccurrence": cooccurrence_path,
        "database_cooccurrence_summary": database_summary_path,
        "amr_mge_context": amr_mge_path,
        "amr_plasmid_context": amr_plasmid_path,
        "amr_mge_same_contig": amr_mge_same_contig_path,
        "amr_plasmid_same_contig": amr_plasmid_same_contig_path,
        "amr_integron_same_contig": amr_integron_same_contig_path,
        "feature_proximity": feature_proximity_path,
        "feature_proximity_all": feature_proximity_all_path,
        "amrfinder_abricate_concordance": concordance_path,
    }


def write_feature_completeness_audit(sample_dir: Path, out_dir: Path) -> dict[str, str]:
    manifest_dir = out_dir / "manifest"
    features_dir = out_dir / "features"
    rows = []
    databases = sorted(KNOWN_DATABASES - {"custom", "ani", "assembly_qc", "quast", "mash"})
    feature_counts = {}
    sample_counts = {}
    unique_counts = {}
    for path in features_dir.glob("*.features.tsv"):
        if path.name == "all_features.tsv":
            continue
        database = path.name.replace(".features.tsv", "")
        table_rows = read_table(path)
        feature_counts[database] = len(table_rows)
        sample_counts[database] = len({row.get("assembly_accession", "") or row.get("sample_id", "") for row in table_rows})
        unique_counts[database] = len({row.get("feature_id", "") for row in table_rows if row.get("feature_id")})
    for database in databases:
        module_dirs = [sample_dir / database]
        if database == "amr":
            module_dirs.extend([sample_dir / "abricate", sample_dir / "tool_results" / "abricate" / "ncbi"])
        elif database in {"vfdb", "plasmidfinder", "isfinder"}:
            module_dirs.append(sample_dir / "tool_results" / "abricate" / database)
        elif database in {"integronfinder", "mobileelementfinder", "mlst"}:
            module_dirs.append(sample_dir / "tool_results" / database)
        raw_dirs = [
            module_dir
            for module_dir in module_dirs
            if module_dir.exists() and any(path.is_file() for path in module_dir.rglob("*"))
        ]
        raw_output_found = bool(raw_dirs)
        module_dir = raw_dirs[0] if raw_dirs else module_dirs[0]
        feature_table = features_dir / f"{database}.features.tsv"
        feature_table_found = feature_table.exists()
        feature_rows = feature_counts.get(database, 0)
        if feature_table_found and feature_rows > 0:
            status = "PASS"
            message = "feature table contains rows"
        elif raw_output_found and not feature_table_found:
            status = "FAIL_MISSING_FEATURE_TABLE"
            message = "raw output found but standardized feature table is missing"
        elif feature_table_found and feature_rows == 0:
            status = "WARNING_EMPTY"
            message = "feature table exists but has zero rows"
        else:
            status = "SKIPPED_NOT_ENABLED"
            message = "module output not found"
        rows.append({
            "database": database,
            "expected_from_profile": "unknown_from_export_context",
            "module_enabled": str(raw_output_found or feature_table_found).lower(),
            "raw_output_found": str(raw_output_found).lower(),
            "feature_table_found": str(feature_table_found).lower(),
            "feature_rows": str(feature_rows),
            "unique_features": str(unique_counts.get(database, 0)),
            "samples_with_features": str(sample_counts.get(database, 0)),
            "samples_processed": str(sample_counts.get(database, 0)),
            "status": status,
            "message": message,
        })
    audit_path = write_rows(
        manifest_dir / "feature_completeness_audit.tsv",
        rows,
        [
            "database",
            "expected_from_profile",
            "module_enabled",
            "raw_output_found",
            "feature_table_found",
            "feature_rows",
            "unique_features",
            "samples_with_features",
            "samples_processed",
            "status",
            "message",
        ],
    )
    native_status_rows = []
    for status_file in sorted(sample_dir.rglob("module_status.tsv")):
        if out_dir in status_file.parents:
            continue
        for row in read_table(status_file):
            native_status_rows.append({
                "module": row.get("module", status_file.parent.name),
                "enabled": row.get("enabled", ""),
                "started": row.get("started", ""),
                "completed": row.get("completed", ""),
                "status": row.get("status", ""),
                "samples_input": row.get("samples_input", ""),
                "samples_processed": row.get("samples_processed", ""),
                "samples_failed": row.get("samples_failed", ""),
                "raw_tables_created": row.get("raw_tables_created", ""),
                "feature_rows_created": row.get("feature_rows_created", ""),
                "unique_features_created": row.get("unique_features_created", ""),
                "output_dir": row.get("output_dir", str(status_file.parent)),
                "message": row.get("message", ""),
            })
    audit_status_rows = [
        {
            "module": row["database"],
            "enabled": row["module_enabled"],
            "started": row["raw_output_found"],
            "completed": row["feature_table_found"],
            "status": row["status"],
            "samples_input": row["samples_processed"],
            "samples_processed": row["samples_processed"],
            "samples_failed": "",
            "raw_tables_created": row["raw_output_found"],
            "feature_rows_created": row["feature_rows"],
            "unique_features_created": row["unique_features"],
            "output_dir": str(sample_dir / row["database"]),
            "message": row["message"],
        }
        for row in rows
    ]
    module_status_path = write_rows(
        manifest_dir / "module_status_summary.tsv",
        native_status_rows + audit_status_rows,
        [
            "module",
            "enabled",
            "started",
            "completed",
            "status",
            "samples_input",
            "samples_processed",
            "samples_failed",
            "raw_tables_created",
            "feature_rows_created",
            "unique_features_created",
            "output_dir",
            "message",
        ],
    )
    return {"feature_completeness_audit": audit_path, "module_status_summary": module_status_path}


def write_report_controls(
    out_dir: Path,
    rows: list[dict[str, str]],
    metadata_rows: list[dict[str, str]],
    large_dataset: bool = False,
    report_mode: str = "publication",
    max_features_heatmap: int = 300,
    max_features_network: int = 300,
    max_metadata_columns: int = 80,
    top_n_features_per_database: int = 25,
    skip_heavy_interactive_plots: bool = False,
) -> dict[str, str]:
    manifest_dir = out_dir / "manifest"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    databases = sorted({row.get("database", "") for row in rows if row.get("database")})
    feature_count = len(rows)
    unique_feature_count = len({(row.get("database", ""), row.get("feature_id", "")) for row in rows if row.get("database") and row.get("feature_id")})
    sample_count = len({metadata_accession(row) for row in metadata_rows if metadata_accession(row)} | {row.get("assembly_accession", "") for row in rows if row.get("assembly_accession")})
    warning = ""
    if large_dataset or feature_count > 10000 or unique_feature_count > max_features_network:
        warning = "large_dataset_summary_mode"
    controls_rows = [
        {
            "setting": "large_dataset",
            "value": str(large_dataset).lower(),
            "message": "Large-dataset safeguards enabled." if large_dataset else "Standard report limits.",
        },
        {"setting": "report_mode", "value": report_mode, "message": "Handoff HTML density preset."},
        {"setting": "max_features_heatmap", "value": str(max_features_heatmap), "message": "Feature cap for handoff presence/absence matrices."},
        {"setting": "max_features_network", "value": str(max_features_network), "message": "Feature cap for report-facing co-occurrence/proximity summaries; complete proximity evidence is preserved as feature_proximity_all.tsv."},
        {"setting": "max_metadata_columns", "value": str(max_metadata_columns), "message": "Metadata rows shown in compact HTML report tables."},
        {"setting": "top_n_features_per_database", "value": str(top_n_features_per_database), "message": "Top prevalent features summarized per database."},
        {"setting": "skip_heavy_interactive_plots", "value": str(skip_heavy_interactive_plots).lower(), "message": "Heavy interactive plots are deprioritized/skipped when supported."},
        {"setting": "samples", "value": str(sample_count), "message": "Samples represented in metadata or feature tables."},
        {"setting": "feature_rows", "value": str(feature_count), "message": "Total standardized feature rows."},
        {"setting": "unique_features", "value": str(unique_feature_count), "message": "Unique database-feature pairs."},
        {"setting": "databases", "value": ",".join(databases), "message": "Databases with standardized feature rows."},
        {"setting": "report_warning", "value": warning, "message": "Summary/report-level warning generated from dataset size and configured limits."},
    ]
    controls_path = write_rows(
        manifest_dir / "report_controls.tsv",
        controls_rows,
        ["setting", "value", "message"],
    )
    return {"report_controls": controls_path}


def write_feature_contract_manifest(out_dir: Path) -> dict[str, str]:
    manifest_dir = out_dir / "manifest"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    contract = {
        "contract_name": "PanR2 standardized feature table",
        "contract_version": CONTRACT_VERSION,
        "schema_version": SCHEMA_VERSION,
        "required_columns": CONTRACT_COLUMNS,
        "optional_columns": OPTIONAL_CONTRACT_COLUMNS,
        "all_columns": FEATURE_COLUMNS,
        "known_databases": sorted(KNOWN_DATABASES),
        "allowed_values": FEATURE_CONTRACT_ALLOWED_VALUES,
        "strict_downstream_layer": "panr2_inputs/features/*.features.tsv",
        "complete_merged_table": "panr2_inputs/features/all_features.tsv",
        "backward_compatibility": (
            "v0.3.x and v0.4.0 feature tables remain valid when all required "
            "columns are present; optional columns should be preserved when supplied."
        ),
    }
    path = manifest_dir / "feature_contract.json"
    path.write_text(json.dumps(contract, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"feature_contract_manifest": str(path)}


def _relative_link(target: Path, base: Path) -> str:
    try:
        return target.relative_to(base).as_posix()
    except ValueError:
        return target.as_posix()


def _html_table(rows: list[dict[str, str]], fields: list[str], max_rows: int = 25) -> str:
    if not rows:
        return "<p>No rows available.</p>"
    header = "".join(f"<th>{html.escape(field)}</th>" for field in fields)
    body_rows = []
    for row in rows[:max_rows]:
        body_rows.append(
            "<tr>"
            + "".join(f"<td>{html.escape(str(row.get(field, '') or ''))}</td>" for field in fields)
            + "</tr>"
        )
    return f"<table><thead><tr>{header}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"


COUNTRY_COORDS = {
    "argentina": (-38.42, -63.62),
    "australia": (-25.27, 133.78),
    "austria": (47.52, 14.55),
    "bangladesh": (23.68, 90.36),
    "belgium": (50.50, 4.47),
    "brazil": (-14.24, -51.93),
    "canada": (56.13, -106.35),
    "chile": (-35.68, -71.54),
    "china": (35.86, 104.20),
    "colombia": (4.57, -74.30),
    "denmark": (56.26, 9.50),
    "egypt": (26.82, 30.80),
    "finland": (61.92, 25.75),
    "france": (46.23, 2.21),
    "germany": (51.17, 10.45),
    "ghana": (7.95, -1.02),
    "greece": (39.07, 21.82),
    "india": (20.59, 78.96),
    "indonesia": (-0.79, 113.92),
    "iran": (32.43, 53.69),
    "iraq": (33.22, 43.68),
    "ireland": (53.41, -8.24),
    "israel": (31.05, 34.85),
    "italy": (41.87, 12.57),
    "japan": (36.20, 138.25),
    "kenya": (-0.02, 37.91),
    "malaysia": (4.21, 101.98),
    "mexico": (23.63, -102.55),
    "nepal": (28.39, 84.12),
    "netherlands": (52.13, 5.29),
    "new zealand": (-40.90, 174.89),
    "nigeria": (9.08, 8.68),
    "norway": (60.47, 8.47),
    "pakistan": (30.38, 69.35),
    "peru": (-9.19, -75.02),
    "philippines": (12.88, 121.77),
    "poland": (51.92, 19.15),
    "portugal": (39.40, -8.22),
    "russia": (61.52, 105.32),
    "saudi arabia": (23.89, 45.08),
    "singapore": (1.35, 103.82),
    "south africa": (-30.56, 22.94),
    "south korea": (35.91, 127.77),
    "spain": (40.46, -3.75),
    "sri lanka": (7.87, 80.77),
    "sweden": (60.13, 18.64),
    "switzerland": (46.82, 8.23),
    "taiwan": (23.70, 120.96),
    "thailand": (15.87, 100.99),
    "turkey": (38.96, 35.24),
    "united kingdom": (55.38, -3.44),
    "uk": (55.38, -3.44),
    "united states": (37.09, -95.71),
    "usa": (37.09, -95.71),
    "vietnam": (14.06, 108.28),
}


BASIC_DATASET_FIELDS = [
    "assembly_accession",
    "sample_id",
    "organism_name",
    "taxid",
    "assembly_level",
    "refseq_category",
    "genome_representation",
    "assembly_release_date",
    "bioproject",
    "biosample",
    "country",
    "continent",
    "subcontinent",
    "collection_year",
    "host",
    "isolation_source",
    "sample_type",
    "environment_medium",
    "submitter",
    "ftp_path_refseq",
    "ftp_path_genbank",
    "qc_status",
    "qc_pass",
    "qc_fail_reasons",
    "genome_size",
    "contig_count",
    "n50",
    "gc_percent",
    "checkm2_completeness",
    "checkm2_contamination",
    "quast_status",
    "ani_status",
    "ani_species_match",
    "mash_status",
    "amr_gene_count",
    "amrfinderplus_gene_count",
    "vfdb_gene_count",
    "plasmidfinder_replicon_count",
    "integronfinder_feature_count",
    "mlst_feature_count",
    "mobsuite_feature_count",
    "genomad_region_count",
    "defensefinder_feature_count",
    "mobileelementfinder_feature_count",
    "amr_genes",
    "amrfinderplus_genes",
    "drug_classes",
    "resistance_mechanisms",
    "vfdb_genes",
    "vfdb_categories",
    "plasmid_replicons",
    "integron_features",
    "mlst_ST",
    "mobsuite_plasmid_types",
    "genomad_regions",
    "defense_systems",
    "mobile_elements",
    "ani_cluster",
    "mash_cluster",
    "dominant_lineage_label",
    "features_detected_databases",
    "modules_run",
    "modules_failed",
    "modules_warning",
    "panresistome_version",
    "feature_contract_version",
]


def _clean_country(value: str) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    return text.split(":", 1)[0].strip()


def _accession_keys(value: str) -> set[str]:
    text = str(value or "").strip()
    if not text:
        return set()
    stem = Path(text).stem
    values = {text, stem}
    if stem.endswith("_genomic"):
        values.add(stem[:-8])
    parts = stem.split("_")
    if len(parts) >= 2 and parts[0] in {"GCF", "GCA"}:
        values.add("_".join(parts[:2]))
        values.add("_".join(parts[:2]).split(".", 1)[0])
    return {value.lower() for value in values if value}


def _index_by_accession(rows: list[dict[str, str]], candidates: list[str]) -> dict[str, dict[str, str]]:
    index = {}
    for row in rows:
        value = first_value(row, candidates, "")
        for key in _accession_keys(value):
            index[key] = row
    return index


def _lookup_by_accession(index: dict[str, dict[str, str]], value: str) -> dict[str, str]:
    for key in _accession_keys(value):
        if key in index:
            return index[key]
    return {}


def _join_values(values: set[str], limit: int = 80) -> str:
    cleaned = sorted({value for value in values if value and not is_missing_value(value)})
    if len(cleaned) > limit:
        return ";".join(cleaned[:limit] + [f"...{len(cleaned) - limit}_more"])
    return ";".join(cleaned)


def _float_text(value: str) -> str:
    numeric = _float_or_none(value)
    if numeric is None:
        return ""
    return f"{numeric:g}"


def _module_summary(out_dir: Path) -> tuple[str, str, str]:
    run = set()
    failed = set()
    warning = set()
    paths = [
        out_dir / "manifest" / "database_setup_status.tsv",
        out_dir / "manifest" / "module_status_summary.tsv",
        out_dir / "manifest" / "native_runner_merge_audit.tsv",
        out_dir / "manifest" / "native_runner_module_status.tsv",
    ]
    for path in [
        path for path in paths if path.exists()
    ]:
        rows = read_table(path)
        run_source = path.name != "database_setup_status.tsv"
        for row in rows:
            name = first_value(row, ["database_or_tool", "module", "database"], "")
            status = first_value(row, ["status"], "").upper()
            if not name:
                continue
            skipped_status = any(token in status for token in ["SKIPPED", "NOT_REQUESTED", "NOT_RUN", "NOT_FOUND", "UNAVAILABLE"])
            if run_source and status and not skipped_status:
                run.add(name)
            if "FAIL" in status or "ERROR" in status:
                failed.add(name)
            if "WARN" in status:
                warning.add(name)
    return ";".join(sorted(run)), ";".join(sorted(failed)), ";".join(sorted(warning))


def write_basic_enriched_dataset(
    sample_dir: Path,
    out_dir: Path,
    output_dir: Path,
    pipeline_version: str = "",
) -> dict[str, str]:
    output_dir.mkdir(parents=True, exist_ok=True)
    metadata_rows = load_metadata_rows(sample_dir)
    normalized_metadata = normalize_metadata_rows(metadata_rows)
    metadata_index = _index_by_accession(metadata_rows, ["Assembly Accession", "assembly_accession", "sequence_accession", "sample_id"])
    normalized_index = _index_by_accession(normalized_metadata, ["assembly_accession", "sample_id"])
    all_features = read_table(out_dir / "features" / "all_features.tsv")
    ani_rows = read_table(out_dir / "ani" / "analysis" / "panr2_ani_summary.csv") or read_table(sample_dir / "ani" / "analysis" / "panr2_ani_summary.csv")
    mash_rows = read_table(sample_dir / "mash" / "analysis" / "closest_mash_neighbor.csv")
    ani_index = _index_by_accession(ani_rows, ["assembly_accession", "sample_id"])
    mash_index = _index_by_accession(mash_rows, ["query"])
    modules_run, modules_failed, modules_warning = _module_summary(out_dir)

    samples = sorted({
        metadata_accession(row) for row in metadata_rows if metadata_accession(row)
    } | {
        row.get("assembly_accession", "") for row in all_features if row.get("assembly_accession")
    })

    by_sample_database: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    detected_databases_by_sample: dict[str, set[str]] = defaultdict(set)
    for feature in all_features:
        sample = feature.get("assembly_accession", "") or feature.get("sample_id", "")
        database = feature.get("database", "")
        if not sample or not database or feature.get("presence", "1") == "0":
            continue
        by_sample_database[(sample, database)].append(feature)
        detected_databases_by_sample[sample].add(database)

    def features(sample: str, database: str) -> list[dict[str, str]]:
        return by_sample_database.get((sample, database), [])

    def feature_ids(sample: str, database: str, categories: set[str] | None = None) -> set[str]:
        values = set()
        for row in features(sample, database):
            if categories and row.get("feature_category", "") not in categories:
                continue
            values.add(row.get("feature_id", ""))
        return values

    rows = []
    for sample in samples:
        raw_metadata = _lookup_by_accession(metadata_index, sample)
        norm_metadata = _lookup_by_accession(normalized_index, sample)
        ani_row = _lookup_by_accession(ani_index, sample)
        mash_row = _lookup_by_accession(mash_index, sample)
        mlst_st = feature_ids(sample, "mlst", {"sequence_type"})
        mlst_st = {value for value in mlst_st if not value.endswith("ST-") and value != "ST_-"}
        ani_cluster = first_value(ani_row, ["feature_id"], "")
        mash_cluster = ""
        if mash_row:
            neighbor = first_value(mash_row, ["reference"], "")
            distance = first_value(mash_row, ["mash_distance"], "")
            mash_cluster = f"nearest:{neighbor};distance:{distance}" if neighbor else ""
        qc_status = first_value(raw_metadata, ["combined_qc_status", "sequence_qc_status"], "")
        row = {
            "assembly_accession": sample,
            "sample_id": first_value(raw_metadata, ["sample_id", "sequence_accession", "Assembly Name"], sample),
            "organism_name": first_value(raw_metadata, ["Organism Name", "organism_name", "Species"], ""),
            "taxid": first_value(raw_metadata, ["Organism Taxonomic ID", "TaxID", "taxid"], ""),
            "assembly_level": first_value(raw_metadata, ["Assembly Level", "assembly_level"], ""),
            "refseq_category": first_value(raw_metadata, ["RefSeq Category", "refseq_category"], ""),
            "genome_representation": first_value(raw_metadata, ["Genome Representation", "genome_representation"], ""),
            "assembly_release_date": first_value(raw_metadata, ["Assembly Release Date", "assembly_release_date"], ""),
            "bioproject": first_value(raw_metadata, ["Assembly BioProject Accession", "bioproject"], ""),
            "biosample": first_value(raw_metadata, ["Assembly BioSample Accession", "biosample"], ""),
            "country": _clean_country(first_value(norm_metadata, ["country"], "") or first_value(raw_metadata, ["Country", "Geographic Location"], "")),
            "continent": first_value(norm_metadata, ["continent"], "") or first_value(raw_metadata, ["Continent"], ""),
            "subcontinent": first_value(norm_metadata, ["subcontinent"], "") or first_value(raw_metadata, ["Subcontinent"], ""),
            "collection_year": first_value(norm_metadata, ["collection_year"], "") or first_value(raw_metadata, ["Collection_Year", "Collection Date"], ""),
            "host": first_value(norm_metadata, ["host"], "") or first_value(raw_metadata, ["Host_SD", "Host"], ""),
            "isolation_source": first_value(norm_metadata, ["isolation_source"], "") or first_value(raw_metadata, ["Isolation_Source_SD", "Isolation Source"], ""),
            "sample_type": first_value(norm_metadata, ["sample_type"], "") or first_value(raw_metadata, ["Sample_Type_SD", "Sample Type"], ""),
            "environment_medium": first_value(norm_metadata, ["environment_medium"], "") or first_value(raw_metadata, ["Environment_Medium_SD", "Environment Medium"], ""),
            "submitter": first_value(raw_metadata, ["Assembly Submitter", "submitter"], ""),
            "ftp_path_refseq": first_value(raw_metadata, ["FTP Path RefSeq", "ftp_path_refseq"], ""),
            "ftp_path_genbank": first_value(raw_metadata, ["FTP Path GenBank", "ftp_path_genbank"], ""),
            "qc_status": qc_status,
            "qc_pass": "true" if qc_status.upper() == "PASS" else "false",
            "qc_fail_reasons": first_value(raw_metadata, ["combined_qc_fail_reasons", "sequence_qc_fail_reasons"], ""),
            "genome_size": first_value(raw_metadata, ["sequence_total_length", "Assembly Stats Total Sequence Length"], ""),
            "contig_count": first_value(raw_metadata, ["sequence_num_contigs", "Assembly Stats Number of Scaffolds"], ""),
            "n50": first_value(raw_metadata, ["sequence_n50", "Assembly Stats Contig N50", "Assembly Stats Scaffold N50"], ""),
            "gc_percent": first_value(raw_metadata, ["sequence_gc_percent", "checkm2_gc_percent"], ""),
            "checkm2_completeness": first_value(raw_metadata, ["checkm2_completeness", "CheckM completeness"], ""),
            "checkm2_contamination": first_value(raw_metadata, ["checkm2_contamination", "CheckM contamination"], ""),
            "quast_status": "PASS" if (out_dir / "assembly_qc" / "analysis" / "panr2_quast_summary.csv").exists() else "",
            "ani_status": "PASS" if ani_row else "",
            "ani_species_match": _float_text(first_value(ani_row, ["identity"], "")),
            "mash_status": "PASS" if mash_row else "",
            "amr_gene_count": str(len(feature_ids(sample, "amr"))),
            "amrfinderplus_gene_count": str(len(feature_ids(sample, "amrfinderplus"))),
            "vfdb_gene_count": str(len(feature_ids(sample, "vfdb"))),
            "plasmidfinder_replicon_count": str(len(feature_ids(sample, "plasmidfinder"))),
            "integronfinder_feature_count": str(len(feature_ids(sample, "integronfinder"))),
            "mlst_feature_count": str(len(feature_ids(sample, "mlst"))),
            "mobsuite_feature_count": str(len(feature_ids(sample, "mobsuite"))),
            "genomad_region_count": str(len(feature_ids(sample, "prophage")) + len(feature_ids(sample, "genomad"))),
            "defensefinder_feature_count": str(len(feature_ids(sample, "defensefinder"))),
            "mobileelementfinder_feature_count": str(len(feature_ids(sample, "mobileelementfinder"))),
            "amr_genes": _join_values(feature_ids(sample, "amr")),
            "amrfinderplus_genes": _join_values(feature_ids(sample, "amrfinderplus")),
            "drug_classes": _join_values({row.get("drug_class") or row.get("feature_category", "") for db in ["amr", "amrfinderplus"] for row in features(sample, db)}),
            "resistance_mechanisms": _join_values({row.get("mechanism", "") for db in ["amr", "amrfinderplus"] for row in features(sample, db)}),
            "vfdb_genes": _join_values(feature_ids(sample, "vfdb")),
            "vfdb_categories": _join_values({row.get("feature_category", "") for row in features(sample, "vfdb")}),
            "plasmid_replicons": _join_values(feature_ids(sample, "plasmidfinder")),
            "integron_features": _join_values(feature_ids(sample, "integronfinder")),
            "mlst_ST": _join_values(mlst_st),
            "mobsuite_plasmid_types": _join_values(feature_ids(sample, "mobsuite")),
            "genomad_regions": _join_values(feature_ids(sample, "prophage") | feature_ids(sample, "genomad")),
            "defense_systems": _join_values(feature_ids(sample, "defensefinder")),
            "mobile_elements": _join_values(feature_ids(sample, "mobileelementfinder")),
            "ani_cluster": ani_cluster,
            "mash_cluster": mash_cluster,
            "dominant_lineage_label": _join_values(mlst_st, limit=3) or ani_cluster or mash_cluster,
            "features_detected_databases": ";".join(sorted(detected_databases_by_sample.get(sample, set()))),
            "modules_run": modules_run,
            "modules_failed": modules_failed,
            "modules_warning": modules_warning,
            "panresistome_version": pipeline_version,
            "feature_contract_version": CONTRACT_VERSION,
        }
        rows.append(row)

    csv_path = output_dir / "enriched_genome_dataset.csv"
    tsv_path = output_dir / "enriched_genome_dataset.tsv"
    write_rows(csv_path, rows, BASIC_DATASET_FIELDS, delimiter=",")
    write_rows(tsv_path, rows, BASIC_DATASET_FIELDS, delimiter="\t")
    return {"basic_enriched_csv": str(csv_path), "basic_enriched_tsv": str(tsv_path)}


def _country_xy(country: str, width: int, height: int) -> tuple[float, float] | None:
    coords = COUNTRY_COORDS.get(_clean_country(country).lower())
    if not coords:
        return None
    lat, lon = coords
    x = (lon + 180.0) / 360.0 * width
    y = (90.0 - lat) / 180.0 * height
    return x, y


def _svg_geographic_map(rows: list[dict[str, str]], title: str) -> str:
    width, height = 960, 480
    points = []
    for row in rows:
        xy = _country_xy(row.get("country", ""), width, height)
        if not xy:
            continue
        prevalence = _float_or_none(row.get("prevalence", "")) or 0.0
        total = int(_float_or_none(row.get("total_genomes", "")) or 0)
        radius = max(5, min(28, 4 + math.sqrt(max(total, 1)) * 4))
        red = int(230 * prevalence)
        blue = int(200 * (1 - prevalence))
        fill = f"rgb({red},80,{blue})"
        label = f"{row.get('country', '')}: {row.get('positive_genomes', '0')}/{row.get('total_genomes', '0')} ({prevalence * 100:.1f}%)"
        points.append(
            f"<circle cx='{xy[0]:.1f}' cy='{xy[1]:.1f}' r='{radius:.1f}' fill='{fill}' "
            "fill-opacity='0.75' stroke='#1f2933' stroke-width='1'>"
            f"<title>{html.escape(label)}</title></circle>"
            f"<text x='{xy[0] + radius + 3:.1f}' y='{xy[1] + 4:.1f}' font-size='11' fill='#1f2933'>{html.escape(row.get('country', ''))}</text>"
        )
    grid = []
    for lon in range(-120, 181, 60):
        x = (lon + 180) / 360 * width
        grid.append(f"<line x1='{x:.1f}' y1='0' x2='{x:.1f}' y2='{height}' stroke='#d9e2ec' stroke-width='1'/>")
    for lat in range(-60, 91, 30):
        y = (90 - lat) / 180 * height
        grid.append(f"<line x1='0' y1='{y:.1f}' x2='{width}' y2='{y:.1f}' stroke='#d9e2ec' stroke-width='1'/>")
    return (
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height + 58}' viewBox='0 0 {width} {height + 58}'>"
        "<rect width='100%' height='100%' fill='#f8fafc'/>"
        f"<text x='20' y='28' font-size='20' font-family='Arial' font-weight='700' fill='#102a43'>{html.escape(title)}</text>"
        f"<g transform='translate(0,48)'><rect x='0' y='0' width='{width}' height='{height}' fill='#eff6ff' stroke='#bcccdc'/>"
        + "".join(grid)
        + "".join(points)
        + "</g></svg>\n"
    )


def _write_png(path: Path, width: int, height: int, pixels: list[bytearray]) -> None:
    def chunk(kind: bytes, data: bytes) -> bytes:
        return (
            struct.pack(">I", len(data))
            + kind
            + data
            + struct.pack(">I", zlib.crc32(kind + data) & 0xFFFFFFFF)
        )

    raw = b"".join(b"\x00" + bytes(row) for row in pixels)
    payload = (
        b"\x89PNG\r\n\x1a\n"
        + chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
        + chunk(b"IDAT", zlib.compress(raw, 9))
        + chunk(b"IEND", b"")
    )
    path.write_bytes(payload)


def _pdf_escape(value: str) -> str:
    return str(value or "").replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def _write_simple_pdf(path: Path, title: str, lines: list[str]) -> None:
    """Write a dependency-free PDF companion for report figures.

    The publication-quality vector artifact remains the SVG; this lightweight PDF
    gives users a portable manuscript/supplement placeholder without requiring
    cairo, matplotlib, or browser-based rendering on remote machines.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    content_lines = ["BT", "/F1 18 Tf", "72 760 Td", f"({_pdf_escape(title)}) Tj", "/F1 10 Tf"]
    y_step = 18
    for line in lines[:34]:
        content_lines.append(f"0 -{y_step} Td")
        content_lines.append(f"({_pdf_escape(line)[:150]}) Tj")
    content_lines.append("ET")
    stream = "\n".join(content_lines).encode("latin-1", errors="replace")
    objects = [
        b"1 0 obj\n<< /Type /Catalog /Pages 2 0 R >>\nendobj\n",
        b"2 0 obj\n<< /Type /Pages /Kids [3 0 R] /Count 1 >>\nendobj\n",
        b"3 0 obj\n<< /Type /Page /Parent 2 0 R /MediaBox [0 0 612 792] /Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>\nendobj\n",
        b"4 0 obj\n<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>\nendobj\n",
        b"5 0 obj\n<< /Length " + str(len(stream)).encode("ascii") + b" >>\nstream\n" + stream + b"\nendstream\nendobj\n",
    ]
    payload = bytearray(b"%PDF-1.4\n")
    offsets = [0]
    for obj in objects:
        offsets.append(len(payload))
        payload.extend(obj)
    xref_offset = len(payload)
    payload.extend(f"xref\n0 {len(objects) + 1}\n".encode("ascii"))
    payload.extend(b"0000000000 65535 f \n")
    for offset in offsets[1:]:
        payload.extend(f"{offset:010d} 00000 n \n".encode("ascii"))
    payload.extend(
        f"trailer\n<< /Size {len(objects) + 1} /Root 1 0 R >>\nstartxref\n{xref_offset}\n%%EOF\n".encode("ascii")
    )
    path.write_bytes(bytes(payload))


def _write_zip_bundle(path: Path, files: list[Path], base_dir: Path) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for file_path in files:
            if file_path.exists() and file_path.is_file():
                archive.write(file_path, arcname=str(file_path.relative_to(base_dir)))
    return str(path)


def _geographic_map_png(rows: list[dict[str, str]], path: Path) -> None:
    width, map_height, header = 960, 480, 48
    height = map_height + header + 10
    pixels = [bytearray([248, 250, 252] * width) for _ in range(height)]

    def set_pixel(x: int, y: int, color: tuple[int, int, int]) -> None:
        if 0 <= x < width and 0 <= y < height:
            idx = x * 3
            pixels[y][idx:idx + 3] = bytes(color)

    def draw_rect(x0: int, y0: int, x1: int, y1: int, color: tuple[int, int, int]) -> None:
        for y in range(max(0, y0), min(height, y1)):
            row = pixels[y]
            for x in range(max(0, x0), min(width, x1)):
                idx = x * 3
                row[idx:idx + 3] = bytes(color)

    def draw_circle(cx: float, cy: float, radius: float, color: tuple[int, int, int]) -> None:
        r2 = radius * radius
        for y in range(int(cy - radius) - 1, int(cy + radius) + 2):
            for x in range(int(cx - radius) - 1, int(cx + radius) + 2):
                if (x - cx) * (x - cx) + (y - cy) * (y - cy) <= r2:
                    set_pixel(x, y, color)

    draw_rect(0, header, width, header + map_height, (239, 246, 255))
    for lon in range(-120, 181, 60):
        x = int((lon + 180) / 360 * width)
        for y in range(header, header + map_height):
            set_pixel(x, y, (217, 226, 236))
    for lat in range(-60, 91, 30):
        y = header + int((90 - lat) / 180 * map_height)
        for x in range(width):
            set_pixel(x, y, (217, 226, 236))

    for row in rows:
        xy = _country_xy(row.get("country", ""), width, map_height)
        if not xy:
            continue
        prevalence = _float_or_none(row.get("prevalence", "")) or 0.0
        total = int(_float_or_none(row.get("total_genomes", "")) or 0)
        radius = max(5, min(28, 4 + math.sqrt(max(total, 1)) * 4))
        color = (int(230 * prevalence), 80, int(200 * (1 - prevalence)))
        draw_circle(xy[0], xy[1] + header, radius, color)
    _write_png(path, width, height, pixels)


def write_important_geographic_outputs(sample_dir: Path, out_dir: Path, important_dir: Path, top_n: int = 20) -> dict[str, str]:
    key_tables = important_dir / "key_tables"
    figures = important_dir / "figures"
    key_tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    metadata_rows = normalize_metadata_rows(load_metadata_rows(sample_dir))
    features = read_table(out_dir / "features" / "all_features.tsv")
    samples = sorted({row.get("assembly_accession", "") for row in metadata_rows if row.get("assembly_accession")})
    metadata_by_sample = {row["assembly_accession"]: row for row in metadata_rows if row.get("assembly_accession")}
    presence = feature_presence(features)
    database_presence: dict[str, set[str]] = defaultdict(set)
    for feature in features:
        sample = feature.get("assembly_accession", "")
        database = feature.get("database", "")
        if sample and database and feature.get("presence", "1") != "0":
            database_presence[database].add(sample)

    feature_rank = sorted(presence.items(), key=lambda item: (-len(item[1]), item[0][0], item[0][1]))
    selected_features = feature_rank[:max(top_n, 0)]
    rows = []

    def add_rows(mode: str, database: str, feature_id: str, present_samples: set[str]):
        for year_mode in ["all_years", "by_year"]:
            groups: dict[tuple[str, str], set[str]] = defaultdict(set)
            positives: dict[tuple[str, str], set[str]] = defaultdict(set)
            for sample in samples:
                meta = metadata_by_sample.get(sample, {})
                country = _clean_country(meta.get("country", ""))
                if not country:
                    continue
                year = meta.get("collection_year", "") if year_mode == "by_year" else "all"
                key = (country, year or "unknown")
                groups[key].add(sample)
                if sample in present_samples:
                    positives[key].add(sample)
            for (country, year), members in sorted(groups.items()):
                positive = len(positives.get((country, year), set()))
                total = len(members)
                meta_example = next((metadata_by_sample.get(sample, {}) for sample in members), {})
                prevalence = positive / total if total else 0.0
                warnings = []
                if total < 5:
                    warnings.append("small_sample_size")
                if not _country_xy(country, 960, 480):
                    warnings.append("missing_map_coordinate")
                rows.append({
                    "mode": mode,
                    "database": database,
                    "feature_id": feature_id,
                    "country": country,
                    "continent": meta_example.get("continent", ""),
                    "subcontinent": meta_example.get("subcontinent", ""),
                    "collection_year": year,
                    "total_genomes": str(total),
                    "positive_genomes": str(positive),
                    "prevalence": f"{prevalence:.4f}",
                    "prevalence_percent": f"{prevalence * 100:.1f}",
                    "warning_flags": ";".join(warnings),
                })

    for database, present_samples in sorted(database_presence.items()):
        add_rows("database_burden", database, "__any_feature__", present_samples)
    for (database, feature_id), present_samples in selected_features:
        add_rows("feature", database, feature_id, present_samples)

    fields = [
        "mode",
        "database",
        "feature_id",
        "country",
        "continent",
        "subcontinent",
        "collection_year",
        "total_genomes",
        "positive_genomes",
        "prevalence",
        "prevalence_percent",
        "warning_flags",
    ]
    data_path = key_tables / "geographic_distribution.tsv"
    write_rows(data_path, rows, fields)
    (figures / "geographic_distribution.data.tsv").write_text(data_path.read_text(encoding="utf-8"), encoding="utf-8")

    first_key = next((row for row in rows if row["collection_year"] == "all"), None)
    first_rows = [
        row for row in rows
        if first_key and row["mode"] == first_key["mode"] and row["database"] == first_key["database"] and row["feature_id"] == first_key["feature_id"] and row["collection_year"] == "all"
    ]
    svg_path = figures / "geographic_distribution_map.svg"
    svg_path.write_text(_svg_geographic_map(first_rows, "Geographic Distribution"), encoding="utf-8")
    png_path = figures / "geographic_distribution_map.png"
    _geographic_map_png(first_rows, png_path)

    datasets = json.dumps(rows)
    html_path = figures / "geographic_distribution_map.html"
    html_path.write_text(
        f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Geographic Distribution</title>
<style>
body {{ font-family: Arial, sans-serif; color: #1f2933; margin: 1.5rem; }}
label {{ font-weight: 700; margin-right: 0.4rem; }}
select {{ margin: 0 1rem 1rem 0; padding: 0.35rem; }}
#map svg {{ max-width: 100%; height: auto; border: 1px solid #d9e2ec; }}
.warning {{ background: #fff7ed; border-left: 4px solid #c2410c; padding: 0.75rem; margin: 1rem 0; }}
</style></head><body>
<h1>Geographic Distribution</h1>
<p>Select a database/feature and year scope. Prevalence reflects this dataset only, not global prevalence.</p>
<div class="warning">Small country/year groups are flagged in the data table. Interpret geographic patterns with BioProject, lineage, and collection-year bias in mind.</div>
<label for="feature">Feature</label><select id="feature"></select>
<label for="year">Year</label><select id="year"></select>
<div id="summary"></div>
<div id="map"></div>
<p><a href="geographic_distribution_map.png">Download initial PNG</a> | <a href="geographic_distribution_map.svg">Download initial SVG</a> | <a href="geographic_distribution.data.tsv">Download plotted data TSV</a></p>
<script>
const rows = {datasets};
const coords = {json.dumps(COUNTRY_COORDS)};
const width = 960, height = 480;
function cleanCountry(value) {{ return (value || '').split(':')[0].trim(); }}
function key(row) {{ return row.mode + ' | ' + row.database + ' | ' + row.feature_id; }}
function xy(country) {{
  const item = coords[cleanCountry(country).toLowerCase()];
  if (!item) return null;
  const lat = item[0], lon = item[1];
  return [((lon + 180) / 360) * width, ((90 - lat) / 180) * height];
}}
const featureSelect = document.getElementById('feature');
const yearSelect = document.getElementById('year');
for (const value of [...new Set(rows.map(key))].sort()) {{
  const opt = document.createElement('option'); opt.value = value; opt.textContent = value; featureSelect.appendChild(opt);
}}
function updateYears() {{
  const selected = featureSelect.value;
  yearSelect.innerHTML = '';
  const years = [...new Set(rows.filter(r => key(r) === selected).map(r => r.collection_year))].sort();
  for (const year of years) {{ const opt = document.createElement('option'); opt.value = year; opt.textContent = year; yearSelect.appendChild(opt); }}
  if (years.includes('all')) yearSelect.value = 'all';
}}
function render() {{
  const selected = featureSelect.value, year = yearSelect.value;
  const active = rows.filter(r => key(r) === selected && r.collection_year === year);
  const total = active.reduce((a, r) => a + Number(r.total_genomes || 0), 0);
  const positive = active.reduce((a, r) => a + Number(r.positive_genomes || 0), 0);
  document.getElementById('summary').innerHTML = `<p><strong>${{selected}}</strong>: ${{positive}} positive observations across ${{total}} country-level genome observations.</p>`;
  let svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${{width}}" height="${{height + 45}}" viewBox="0 0 ${{width}} ${{height + 45}}">`;
  svg += `<rect width="100%" height="100%" fill="#f8fafc"/><text x="20" y="28" font-size="20" font-family="Arial" font-weight="700" fill="#102a43">Geographic Distribution</text>`;
  svg += `<g transform="translate(0,45)"><rect x="0" y="0" width="${{width}}" height="${{height}}" fill="#eff6ff" stroke="#bcccdc"/>`;
  for (let lon = -120; lon <= 180; lon += 60) {{ const x = ((lon + 180) / 360) * width; svg += `<line x1="${{x}}" y1="0" x2="${{x}}" y2="${{height}}" stroke="#d9e2ec"/>`; }}
  for (let lat = -60; lat <= 90; lat += 30) {{ const y = ((90 - lat) / 180) * height; svg += `<line x1="0" y1="${{y}}" x2="${{width}}" y2="${{y}}" stroke="#d9e2ec"/>`; }}
  for (const row of active) {{
    const point = xy(row.country); if (!point) continue;
    const prevalence = Number(row.prevalence || 0), total = Number(row.total_genomes || 0);
    const radius = Math.max(5, Math.min(28, 4 + Math.sqrt(Math.max(total, 1)) * 4));
    const fill = `rgb(${{Math.round(230 * prevalence)}},80,${{Math.round(200 * (1 - prevalence))}})`;
    const label = `${{row.country}}: ${{row.positive_genomes}}/${{row.total_genomes}} (${{row.prevalence_percent}}%)`;
    svg += `<circle cx="${{point[0].toFixed(1)}}" cy="${{point[1].toFixed(1)}}" r="${{radius.toFixed(1)}}" fill="${{fill}}" fill-opacity="0.75" stroke="#1f2933"><title>${{label}}</title></circle>`;
    svg += `<text x="${{(point[0] + radius + 3).toFixed(1)}}" y="${{(point[1] + 4).toFixed(1)}}" font-size="11" fill="#1f2933">${{row.country}}</text>`;
  }}
  svg += '</g></svg>';
  document.getElementById('map').innerHTML = svg;
}}
featureSelect.addEventListener('change', () => {{ updateYears(); render(); }});
yearSelect.addEventListener('change', render);
updateYears(); render();
</script></body></html>
""",
        encoding="utf-8",
    )
    return {
        "important_geographic_distribution": str(data_path),
        "important_geographic_map_html": str(html_path),
        "important_geographic_map_svg": str(svg_path),
        "important_geographic_map_png": str(png_path),
    }


def _write_bar_svg(path: Path, rows: list[dict[str, str]], title: str, label_field: str, value_field: str, x_label: str = "") -> None:
    width = 960
    row_height = 28
    top = 58
    left = 250
    plot_width = 620
    height = max(180, top + row_height * max(len(rows), 1) + 30)
    values = [_float_or_none(row.get(value_field, "")) or 0.0 for row in rows]
    max_value = max(values) if values else 1.0
    if max_value <= 0:
        max_value = 1.0
    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='#f8fafc'/>",
        f"<text x='20' y='30' font-family='Arial' font-size='20' font-weight='700' fill='#102a43'>{html.escape(title)}</text>",
    ]
    for idx, row in enumerate(rows):
        y = top + idx * row_height
        value = _float_or_none(row.get(value_field, "")) or 0.0
        bar_width = value / max_value * plot_width
        label = row.get(label_field, "")
        if len(label) > 38:
            label = label[:35] + "..."
        parts.append(f"<text x='20' y='{y + 17}' font-family='Arial' font-size='12' fill='#1f2933'>{html.escape(label)}</text>")
        parts.append(f"<rect x='{left}' y='{y}' width='{bar_width:.1f}' height='18' fill='#0f766e'/>")
        parts.append(f"<text x='{left + bar_width + 8:.1f}' y='{y + 14}' font-family='Arial' font-size='12' fill='#1f2933'>{html.escape(row.get(value_field, ''))}</text>")
    if x_label:
        parts.append(f"<text x='{left}' y='{height - 10}' font-family='Arial' font-size='12' fill='#52606d'>{html.escape(x_label)}</text>")
    parts.append("</svg>\n")
    path.write_text("".join(parts), encoding="utf-8")


def _write_bar_png(path: Path, rows: list[dict[str, str]], value_field: str) -> None:
    width = 960
    row_height = 28
    top = 58
    left = 250
    plot_width = 620
    height = max(180, top + row_height * max(len(rows), 1) + 30)
    pixels = [bytearray([248, 250, 252] * width) for _ in range(height)]

    def rect(x0: int, y0: int, x1: int, y1: int, color: tuple[int, int, int]) -> None:
        for y in range(max(0, y0), min(height, y1)):
            row_pixels = pixels[y]
            for x in range(max(0, x0), min(width, x1)):
                idx = x * 3
                row_pixels[idx:idx + 3] = bytes(color)

    values = [_float_or_none(row.get(value_field, "")) or 0.0 for row in rows]
    max_value = max(values) if values else 1.0
    if max_value <= 0:
        max_value = 1.0
    for idx, value in enumerate(values):
        y = top + idx * row_height
        bar_width = int(value / max_value * plot_width)
        rect(left, y, left + bar_width, y + 18, (15, 118, 110))
    _write_png(path, width, height, pixels)


def _write_line_svg(path: Path, rows: list[dict[str, str]], title: str, label_field: str, value_field: str) -> None:
    width, height = 960, 420
    left, top, plot_width, plot_height = 70, 60, 820, 280
    values = [_float_or_none(row.get(value_field, "")) or 0.0 for row in rows]
    max_value = max(values) if values else 1.0
    if max_value <= 0:
        max_value = 1.0
    points = []
    denom = max(len(values) - 1, 1)
    for idx, value in enumerate(values):
        x = left + idx / denom * plot_width
        y = top + plot_height - value / max_value * plot_height
        points.append((x, y, value))
    point_text = " ".join(f"{x:.1f},{y:.1f}" for x, y, _ in points)
    circles = "".join(f"<circle cx='{x:.1f}' cy='{y:.1f}' r='4' fill='#0f766e'><title>{value:g}</title></circle>" for x, y, value in points)
    path.write_text(
        f"""<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>
<rect width='100%' height='100%' fill='#f8fafc'/>
<text x='20' y='30' font-family='Arial' font-size='20' font-weight='700' fill='#102a43'>{html.escape(title)}</text>
<line x1='{left}' y1='{top + plot_height}' x2='{left + plot_width}' y2='{top + plot_height}' stroke='#9fb3c8'/>
<line x1='{left}' y1='{top}' x2='{left}' y2='{top + plot_height}' stroke='#9fb3c8'/>
<polyline points='{point_text}' fill='none' stroke='#0f766e' stroke-width='3'/>
{circles}
<text x='{left}' y='{height - 30}' font-family='Arial' font-size='12' fill='#52606d'>{html.escape(label_field)}</text>
</svg>
""",
        encoding="utf-8",
    )


def _write_line_png(path: Path, rows: list[dict[str, str]], value_field: str) -> None:
    width, height = 960, 420
    left, top, plot_width, plot_height = 70, 60, 820, 280
    pixels = [bytearray([248, 250, 252] * width) for _ in range(height)]

    def set_pixel(x: int, y: int, color: tuple[int, int, int]) -> None:
        if 0 <= x < width and 0 <= y < height:
            idx = x * 3
            pixels[y][idx:idx + 3] = bytes(color)

    def line(x0: float, y0: float, x1: float, y1: float, color: tuple[int, int, int]) -> None:
        steps = int(max(abs(x1 - x0), abs(y1 - y0), 1))
        for i in range(steps + 1):
            t = i / steps
            set_pixel(int(x0 + (x1 - x0) * t), int(y0 + (y1 - y0) * t), color)

    values = [_float_or_none(row.get(value_field, "")) or 0.0 for row in rows]
    max_value = max(values) if values else 1.0
    if max_value <= 0:
        max_value = 1.0
    points = []
    denom = max(len(values) - 1, 1)
    for idx, value in enumerate(values):
        x = left + idx / denom * plot_width
        y = top + plot_height - value / max_value * plot_height
        points.append((x, y))
    for a, b in zip(points, points[1:]):
        line(a[0], a[1], b[0], b[1], (15, 118, 110))
    _write_png(path, width, height, pixels)


def _write_heatmap_svg(path: Path, rows: list[dict[str, str]], title: str, row_field: str, column_field: str, value_field: str) -> None:
    row_labels = list(dict.fromkeys(row.get(row_field, "") for row in rows if row.get(row_field, "")))
    column_labels = sorted({row.get(column_field, "") for row in rows if row.get(column_field, "")})
    values = {(row.get(row_field, ""), row.get(column_field, "")): _float_or_none(row.get(value_field, "")) or 0.0 for row in rows}
    cell_w = 68
    cell_h = 22
    left = 260
    top = 72
    width = max(720, left + cell_w * max(len(column_labels), 1) + 40)
    height = max(180, top + cell_h * max(len(row_labels), 1) + 40)
    max_value = max(values.values()) if values else 1.0
    if max_value <= 0:
        max_value = 1.0
    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='#f8fafc'/>",
        f"<text x='20' y='30' font-family='Arial' font-size='20' font-weight='700' fill='#102a43'>{html.escape(title)}</text>",
    ]
    for cidx, column in enumerate(column_labels):
        x = left + cidx * cell_w
        parts.append(f"<text x='{x + 6}' y='{top - 12}' font-family='Arial' font-size='11' fill='#52606d'>{html.escape(column)}</text>")
    for ridx, label in enumerate(row_labels):
        y = top + ridx * cell_h
        short_label = label if len(label) <= 42 else label[:39] + "..."
        parts.append(f"<text x='20' y='{y + 15}' font-family='Arial' font-size='11' fill='#1f2933'>{html.escape(short_label)}</text>")
        for cidx, column in enumerate(column_labels):
            x = left + cidx * cell_w
            value = values.get((label, column), 0.0)
            intensity = value / max_value
            red = int(232 - 160 * intensity)
            green = int(246 - 82 * intensity)
            blue = int(255 - 105 * intensity)
            parts.append(f"<rect x='{x}' y='{y}' width='{cell_w - 2}' height='{cell_h - 2}' fill='rgb({red},{green},{blue})' stroke='#d9e2ec'><title>{html.escape(label)} {html.escape(column)}: {value:g}</title></rect>")
    parts.append("</svg>\n")
    path.write_text("".join(parts), encoding="utf-8")


def _write_heatmap_png(path: Path, rows: list[dict[str, str]], row_field: str, column_field: str, value_field: str) -> None:
    row_labels = list(dict.fromkeys(row.get(row_field, "") for row in rows if row.get(row_field, "")))
    column_labels = sorted({row.get(column_field, "") for row in rows if row.get(column_field, "")})
    values = {(row.get(row_field, ""), row.get(column_field, "")): _float_or_none(row.get(value_field, "")) or 0.0 for row in rows}
    cell_w = 68
    cell_h = 22
    left = 260
    top = 72
    width = max(720, left + cell_w * max(len(column_labels), 1) + 40)
    height = max(180, top + cell_h * max(len(row_labels), 1) + 40)
    max_value = max(values.values()) if values else 1.0
    if max_value <= 0:
        max_value = 1.0
    pixels = [bytearray([248, 250, 252] * width) for _ in range(height)]

    def rect(x0: int, y0: int, x1: int, y1: int, color: tuple[int, int, int]) -> None:
        for y in range(max(0, y0), min(height, y1)):
            row_pixels = pixels[y]
            for x in range(max(0, x0), min(width, x1)):
                idx = x * 3
                row_pixels[idx:idx + 3] = bytes(color)

    for ridx, label in enumerate(row_labels):
        y = top + ridx * cell_h
        for cidx, column in enumerate(column_labels):
            x = left + cidx * cell_w
            intensity = values.get((label, column), 0.0) / max_value
            red = int(232 - 160 * intensity)
            green = int(246 - 82 * intensity)
            blue = int(255 - 105 * intensity)
            rect(x, y, x + cell_w - 2, y + cell_h - 2, (red, green, blue))
    _write_png(path, width, height, pixels)


def _write_temporal_series_svg(path: Path, rows: list[dict[str, str]], title: str) -> None:
    width, height = 960, 420
    left, top, plot_width, plot_height = 78, 58, 812, 280
    years = [_float_or_none(row.get("collection_year", "")) for row in rows]
    values = [_float_or_none(row.get("prevalence_percent", "")) for row in rows]
    pairs = [(year, value, row) for year, value, row in zip(years, values, rows) if year is not None and value is not None]
    if not pairs:
        pairs = [(0.0, 0.0, {})]
    min_year = min(year for year, _, _ in pairs)
    max_year = max(year for year, _, _ in pairs)
    if min_year == max_year:
        max_year = min_year + 1
    max_value = max(max(value for _, value, _ in pairs), 100.0)
    points = []
    circles = []
    labels = []
    for year, value, row in pairs:
        x = left + (year - min_year) / (max_year - min_year) * plot_width
        y = top + plot_height - value / max_value * plot_height
        points.append(f"{x:.1f},{y:.1f}")
        label = f"{int(year)}: {row.get('prevalence_percent', '0')}% ({row.get('positive_genomes', '0')}/{row.get('total_genomes', '0')})"
        circles.append(f"<circle cx='{x:.1f}' cy='{y:.1f}' r='5' fill='#0f766e'><title>{html.escape(label)}</title></circle>")
        labels.append(f"<text x='{x - 16:.1f}' y='{top + plot_height + 22}' font-family='Arial' font-size='11' fill='#52606d'>{int(year)}</text>")
    path.write_text(
        f"""<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>
<rect width='100%' height='100%' fill='#f8fafc'/>
<text x='20' y='30' font-family='Arial' font-size='20' font-weight='700' fill='#102a43'>{html.escape(title)}</text>
<line x1='{left}' y1='{top + plot_height}' x2='{left + plot_width}' y2='{top + plot_height}' stroke='#9fb3c8'/>
<line x1='{left}' y1='{top}' x2='{left}' y2='{top + plot_height}' stroke='#9fb3c8'/>
<text x='20' y='{top + 18}' font-family='Arial' font-size='12' fill='#52606d'>Prevalence %</text>
<polyline points='{" ".join(points)}' fill='none' stroke='#0f766e' stroke-width='3'/>
{"".join(circles)}
{"".join(labels)}
</svg>
""",
        encoding="utf-8",
    )


def _write_temporal_series_png(path: Path, rows: list[dict[str, str]]) -> None:
    width, height = 960, 420
    left, top, plot_width, plot_height = 78, 58, 812, 280
    pixels = [bytearray([248, 250, 252] * width) for _ in range(height)]

    def set_pixel(x: int, y: int, color: tuple[int, int, int]) -> None:
        if 0 <= x < width and 0 <= y < height:
            idx = x * 3
            pixels[y][idx:idx + 3] = bytes(color)

    def line(x0: float, y0: float, x1: float, y1: float, color: tuple[int, int, int]) -> None:
        steps = int(max(abs(x1 - x0), abs(y1 - y0), 1))
        for i in range(steps + 1):
            t = i / steps
            set_pixel(int(x0 + (x1 - x0) * t), int(y0 + (y1 - y0) * t), color)

    def circle(cx: float, cy: float, radius: float, color: tuple[int, int, int]) -> None:
        r2 = radius * radius
        for y in range(int(cy - radius) - 1, int(cy + radius) + 2):
            for x in range(int(cx - radius) - 1, int(cx + radius) + 2):
                if (x - cx) ** 2 + (y - cy) ** 2 <= r2:
                    set_pixel(x, y, color)

    pairs = []
    for row in rows:
        year = _float_or_none(row.get("collection_year", ""))
        value = _float_or_none(row.get("prevalence_percent", ""))
        if year is not None and value is not None:
            pairs.append((year, value))
    if not pairs:
        pairs = [(0.0, 0.0)]
    min_year = min(year for year, _ in pairs)
    max_year = max(year for year, _ in pairs)
    if min_year == max_year:
        max_year = min_year + 1
    max_value = max(max(value for _, value in pairs), 100.0)
    line(left, top + plot_height, left + plot_width, top + plot_height, (159, 179, 200))
    line(left, top, left, top + plot_height, (159, 179, 200))
    points = []
    for year, value in pairs:
        x = left + (year - min_year) / (max_year - min_year) * plot_width
        y = top + plot_height - value / max_value * plot_height
        points.append((x, y))
        circle(x, y, 5, (15, 118, 110))
    for a, b in zip(points, points[1:]):
        line(a[0], a[1], b[0], b[1], (15, 118, 110))
    _write_png(path, width, height, pixels)


def _write_temporal_slope_svg(path: Path, rows: list[dict[str, str]], title: str) -> None:
    width = 960
    top = 62
    left_x = 250
    right_x = 710
    row_h = 28
    height = max(220, top + row_h * max(len(rows), 1) + 50)
    values = []
    for row in rows:
        first = _float_or_none(row.get("first_year_prevalence_percent", "")) or 0.0
        last = _float_or_none(row.get("last_year_prevalence_percent", "")) or 0.0
        values.extend([first, last])
    max_value = max(max(values), 100.0) if values else 100.0

    def y_for(idx: int, value: float) -> float:
        base = top + idx * row_h + 14
        return base - (value / max_value) * 8

    def color(label: str) -> str:
        if "increasing" in label:
            return "#0f766e"
        if "decreasing" in label:
            return "#b91c1c"
        return "#64748b"

    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='#f8fafc'/>",
        f"<text x='20' y='30' font-family='Arial' font-size='20' font-weight='700' fill='#102a43'>{html.escape(title)}</text>",
        f"<text x='{left_x - 24}' y='52' font-family='Arial' font-size='12' fill='#52606d'>First year</text>",
        f"<text x='{right_x - 20}' y='52' font-family='Arial' font-size='12' fill='#52606d'>Last year</text>",
    ]
    for idx, row in enumerate(rows):
        first = _float_or_none(row.get("first_year_prevalence_percent", "")) or 0.0
        last = _float_or_none(row.get("last_year_prevalence_percent", "")) or 0.0
        y1 = y_for(idx, first)
        y2 = y_for(idx, last)
        label = f"{row.get('database', '')}:{row.get('feature_id', '')}"
        short_label = label if len(label) <= 34 else label[:31] + "..."
        stroke = color(row.get("trend_label", ""))
        parts.append(f"<text x='20' y='{top + idx * row_h + 18}' font-family='Arial' font-size='11' fill='#1f2933'>{html.escape(short_label)}</text>")
        parts.append(f"<line x1='{left_x}' y1='{y1:.1f}' x2='{right_x}' y2='{y2:.1f}' stroke='{stroke}' stroke-width='2.4'><title>{html.escape(label)}: {first:.1f}% to {last:.1f}%</title></line>")
        parts.append(f"<circle cx='{left_x}' cy='{y1:.1f}' r='4' fill='{stroke}'/>")
        parts.append(f"<circle cx='{right_x}' cy='{y2:.1f}' r='4' fill='{stroke}'/>")
        parts.append(f"<text x='{left_x + 8}' y='{y1 + 4:.1f}' font-family='Arial' font-size='10' fill='#52606d'>{first:.1f}%</text>")
        parts.append(f"<text x='{right_x + 8}' y='{y2 + 4:.1f}' font-family='Arial' font-size='10' fill='#52606d'>{last:.1f}%</text>")
    parts.append("</svg>\n")
    path.write_text("".join(parts), encoding="utf-8")


def _write_temporal_slope_png(path: Path, rows: list[dict[str, str]]) -> None:
    width = 960
    top = 62
    left_x = 250
    right_x = 710
    row_h = 28
    height = max(220, top + row_h * max(len(rows), 1) + 50)
    pixels = [bytearray([248, 250, 252] * width) for _ in range(height)]

    def set_pixel(x: int, y: int, color: tuple[int, int, int]) -> None:
        if 0 <= x < width and 0 <= y < height:
            idx = x * 3
            pixels[y][idx:idx + 3] = bytes(color)

    def line(x0: float, y0: float, x1: float, y1: float, color: tuple[int, int, int]) -> None:
        steps = int(max(abs(x1 - x0), abs(y1 - y0), 1))
        for i in range(steps + 1):
            t = i / steps
            set_pixel(int(x0 + (x1 - x0) * t), int(y0 + (y1 - y0) * t), color)

    def color(label: str) -> tuple[int, int, int]:
        if "increasing" in label:
            return (15, 118, 110)
        if "decreasing" in label:
            return (185, 28, 28)
        return (100, 116, 139)

    for idx, row in enumerate(rows):
        first = _float_or_none(row.get("first_year_prevalence_percent", "")) or 0.0
        last = _float_or_none(row.get("last_year_prevalence_percent", "")) or 0.0
        y1 = top + idx * row_h + 14 - first / 100.0 * 8
        y2 = top + idx * row_h + 14 - last / 100.0 * 8
        line(left_x, y1, right_x, y2, color(row.get("trend_label", "")))
    _write_png(path, width, height, pixels)


def _safe_filename(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value or "").strip())
    return text.strip("_") or "unknown"


def _cooccurrence_direction(phi: float, q_value: float | None, min_abs_phi: float = 0.2, q_threshold: float = 0.05) -> tuple[str, str]:
    if q_value is None or q_value > q_threshold or abs(phi) < min_abs_phi:
        return "not_significant", "not_significant"
    if phi > 0:
        return "positive", "significant_positive"
    if phi < 0:
        return "negative", "significant_negative"
    return "not_significant", "not_significant"


def _context_evidence_label(value: str) -> str:
    text = str(value or "").lower()
    if "overlap" in text or "adjacent" in text or "level_4" in text:
        return "overlap_or_adjacent"
    if "within_10kb" in text or "level_3" in text:
        return "within_10kb"
    if "within_50kb" in text:
        return "within_50kb"
    if "same_contig" in text or "level_2" in text:
        return "same_contig"
    if "same_genome" in text:
        return "same_genome"
    return "unknown"


def _database_color(database: str) -> str:
    palette = {
        "amr": "#dc2626",
        "amrfinderplus": "#ef4444",
        "vfdb": "#16a34a",
        "plasmidfinder": "#7c3aed",
        "integronfinder": "#ea580c",
        "mlst": "#2563eb",
        "prophage": "#0891b2",
        "genomad": "#0891b2",
        "mobsuite": "#9333ea",
        "isfinder": "#ca8a04",
        "mobileelementfinder": "#ca8a04",
        "defensefinder": "#475569",
    }
    return palette.get(str(database or "").lower(), "#64748b")


def _cooccurrence_cell_color(row: dict[str, str]) -> str:
    label = row.get("significance_label", "")
    phi = _float_or_none(row.get("phi_correlation", "")) or 0.0
    intensity = min(abs(phi), 1.0)
    if label != "significant_positive" and label != "significant_negative":
        return "#f1f5f9"
    if label == "significant_positive":
        red = 254
        green = int(226 - 140 * intensity)
        blue = int(226 - 140 * intensity)
        return f"rgb({red},{green},{blue})"
    red = int(219 - 150 * intensity)
    green = int(234 - 120 * intensity)
    blue = 254
    return f"rgb({red},{green},{blue})"


def _write_cooccurrence_heatmap_svg(path: Path, rows: list[dict[str, str]], title: str) -> None:
    x_features = list(dict.fromkeys(row.get("feature_a_id", "") for row in rows if row.get("feature_a_id")))
    y_features = list(dict.fromkeys(row.get("feature_b_id", "") for row in rows if row.get("feature_b_id")))
    cell = 34
    left = 230
    top = 170
    width = max(820, left + cell * max(len(x_features), 1) + 70)
    height = max(320, top + cell * max(len(y_features), 1) + 50)
    by_pair = {(row.get("feature_a_id", ""), row.get("feature_b_id", "")): row for row in rows}
    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='#f8fafc'/>",
        f"<text x='20' y='32' font-family='Arial' font-size='20' font-weight='700' fill='#102a43'>{html.escape(title)}</text>",
        "<text x='20' y='60' font-family='Arial' font-size='12' fill='#52606d'>Red = significant positive association; blue = significant negative association; gray = not significant or below support threshold.</text>",
    ]
    for x_idx, feature in enumerate(x_features):
        x = left + x_idx * cell + 18
        label = feature if len(feature) <= 22 else feature[:19] + "..."
        parts.append(f"<text x='{x}' y='{top - 8}' font-family='Arial' font-size='10' fill='#1f2933' transform='rotate(-60 {x} {top - 8})'>{html.escape(label)}</text>")
    for y_idx, feature in enumerate(y_features):
        y = top + y_idx * cell + 22
        label = feature if len(feature) <= 28 else feature[:25] + "..."
        parts.append(f"<text x='20' y='{y}' font-family='Arial' font-size='10' fill='#1f2933'>{html.escape(label)}</text>")
    for y_idx, y_feature in enumerate(y_features):
        for x_idx, x_feature in enumerate(x_features):
            row = by_pair.get((x_feature, y_feature), {})
            x = left + x_idx * cell
            y = top + y_idx * cell
            fill = _cooccurrence_cell_color(row)
            phi = _float_or_none(row.get("phi_correlation", "")) if row else None
            label = "" if phi is None or row.get("significance_label") not in {"significant_positive", "significant_negative"} else f"{phi:.2f}"
            tip = (
                f"{x_feature} vs {y_feature}; phi={row.get('phi_correlation', '')}; q={row.get('q_value', '')}; "
                f"both={row.get('n_both_present', '')}/{row.get('n_total', '')}; {row.get('significance_label', '')}"
            ) if row else f"{x_feature} vs {y_feature}: no pair"
            parts.append(f"<rect x='{x}' y='{y}' width='{cell - 2}' height='{cell - 2}' fill='{fill}' stroke='#cbd5e1'><title>{html.escape(tip)}</title></rect>")
            if label:
                parts.append(f"<text x='{x + 4}' y='{y + 20}' font-family='Arial' font-size='10' fill='#111827'>{html.escape(label)}</text>")
    parts.append("</svg>\n")
    path.write_text("".join(parts), encoding="utf-8")


def _write_cooccurrence_heatmap_png(path: Path, rows: list[dict[str, str]]) -> None:
    x_features = list(dict.fromkeys(row.get("feature_a_id", "") for row in rows if row.get("feature_a_id")))
    y_features = list(dict.fromkeys(row.get("feature_b_id", "") for row in rows if row.get("feature_b_id")))
    cell = 34
    left = 230
    top = 170
    width = max(820, left + cell * max(len(x_features), 1) + 70)
    height = max(320, top + cell * max(len(y_features), 1) + 50)
    pixels = [bytearray([248, 250, 252] * width) for _ in range(height)]
    by_pair = {(row.get("feature_a_id", ""), row.get("feature_b_id", "")): row for row in rows}

    def rect(x0: int, y0: int, x1: int, y1: int, color: tuple[int, int, int]) -> None:
        for y in range(max(0, y0), min(height, y1)):
            row_pixels = pixels[y]
            for x in range(max(0, x0), min(width, x1)):
                idx = x * 3
                row_pixels[idx:idx + 3] = bytes(color)

    def parse_color(row: dict[str, str]) -> tuple[int, int, int]:
        label = row.get("significance_label", "")
        phi = abs(_float_or_none(row.get("phi_correlation", "")) or 0.0)
        intensity = min(phi, 1.0)
        if label == "significant_positive":
            return (254, int(226 - 140 * intensity), int(226 - 140 * intensity))
        if label == "significant_negative":
            return (int(219 - 150 * intensity), int(234 - 120 * intensity), 254)
        return (241, 245, 249)

    for y_idx, y_feature in enumerate(y_features):
        for x_idx, x_feature in enumerate(x_features):
            x = left + x_idx * cell
            y = top + y_idx * cell
            rect(x, y, x + cell - 2, y + cell - 2, parse_color(by_pair.get((x_feature, y_feature), {})))
    _write_png(path, width, height, pixels)


def _write_cooccurrence_network_svg(path: Path, nodes: list[dict[str, str]], edges: list[dict[str, str]], title: str) -> None:
    width, height = 920, 620
    cx, cy = width / 2, height / 2 + 20
    radius = 230
    positions: dict[str, tuple[float, float]] = {}
    node_ids = [row.get("node_id", "") for row in nodes if row.get("node_id")]
    for idx, node_id in enumerate(node_ids):
        angle = 2 * math.pi * idx / max(len(node_ids), 1)
        positions[node_id] = (cx + math.cos(angle) * radius, cy + math.sin(angle) * radius)
    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='#f8fafc'/>",
        f"<text x='20' y='32' font-family='Arial' font-size='20' font-weight='700' fill='#102a43'>{html.escape(title)}</text>",
        "<text x='20' y='56' font-family='Arial' font-size='12' fill='#52606d'>Edges show significant sample-level co-occurrence unless stronger evidence is listed in the edge table.</text>",
    ]
    for edge in edges:
        source = edge.get("source_feature", "")
        target = edge.get("target_feature", "")
        if source not in positions or target not in positions:
            continue
        x1, y1 = positions[source]
        x2, y2 = positions[target]
        phi = abs(_float_or_none(edge.get("phi_correlation", "")) or 0.0)
        color = "#dc2626" if edge.get("direction") == "positive" else "#2563eb"
        width_px = 1.0 + 5.0 * min(phi, 1.0)
        dash = " stroke-dasharray='5 4'" if edge.get("evidence_level") == "same_genome" else ""
        tip = f"{source} - {target}; phi={edge.get('phi_correlation', '')}; q={edge.get('q_value', '')}; both={edge.get('n_both_present', '')}"
        parts.append(f"<line x1='{x1:.1f}' y1='{y1:.1f}' x2='{x2:.1f}' y2='{y2:.1f}' stroke='{color}' stroke-width='{width_px:.1f}' opacity='0.72'{dash}><title>{html.escape(tip)}</title></line>")
    for node in nodes:
        node_id = node.get("node_id", "")
        if node_id not in positions:
            continue
        x, y = positions[node_id]
        prevalence = _float_or_none(node.get("prevalence", "")) or 0.0
        node_radius = max(6, min(24, 6 + prevalence * 22))
        color = _database_color(node.get("database", ""))
        label = node.get("node_label", node_id)
        label = label if len(label) <= 22 else label[:19] + "..."
        parts.append(f"<circle cx='{x:.1f}' cy='{y:.1f}' r='{node_radius:.1f}' fill='{color}' fill-opacity='0.85' stroke='#1f2933'><title>{html.escape(node_id)} prevalence={prevalence:.3f}</title></circle>")
        parts.append(f"<text x='{x + node_radius + 3:.1f}' y='{y + 4:.1f}' font-family='Arial' font-size='10' fill='#1f2933'>{html.escape(label)}</text>")
    parts.append("</svg>\n")
    path.write_text("".join(parts), encoding="utf-8")


def _write_cooccurrence_network_png(path: Path, nodes: list[dict[str, str]], edges: list[dict[str, str]]) -> None:
    width, height = 920, 620
    cx, cy = width / 2, height / 2 + 20
    radius = 230
    pixels = [bytearray([248, 250, 252] * width) for _ in range(height)]
    node_ids = [row.get("node_id", "") for row in nodes if row.get("node_id")]
    positions: dict[str, tuple[float, float]] = {}
    for idx, node_id in enumerate(node_ids):
        angle = 2 * math.pi * idx / max(len(node_ids), 1)
        positions[node_id] = (cx + math.cos(angle) * radius, cy + math.sin(angle) * radius)

    def set_pixel(x: int, y: int, color: tuple[int, int, int]) -> None:
        if 0 <= x < width and 0 <= y < height:
            idx = x * 3
            pixels[y][idx:idx + 3] = bytes(color)

    def line(x0: float, y0: float, x1: float, y1: float, color: tuple[int, int, int]) -> None:
        steps = int(max(abs(x1 - x0), abs(y1 - y0), 1))
        for i in range(steps + 1):
            t = i / steps
            set_pixel(int(x0 + (x1 - x0) * t), int(y0 + (y1 - y0) * t), color)

    def circle(cx0: float, cy0: float, r: float, color: tuple[int, int, int]) -> None:
        r2 = r * r
        for y in range(int(cy0 - r) - 1, int(cy0 + r) + 2):
            for x in range(int(cx0 - r) - 1, int(cx0 + r) + 2):
                if (x - cx0) ** 2 + (y - cy0) ** 2 <= r2:
                    set_pixel(x, y, color)

    for edge in edges:
        source = edge.get("source_feature", "")
        target = edge.get("target_feature", "")
        if source in positions and target in positions:
            line(*positions[source], *positions[target], (220, 38, 38) if edge.get("direction") == "positive" else (37, 99, 235))
    for node in nodes:
        node_id = node.get("node_id", "")
        if node_id in positions:
            circle(*positions[node_id], 9, (15, 118, 110))
    _write_png(path, width, height, pixels)


def _write_context_ladder_svg(path: Path, rows: list[dict[str, str]], title: str) -> None:
    ladder = [
        ("same_genome", "Same genome"),
        ("same_contig", "Same contig"),
        ("within_50kb", "Within 50 kb"),
        ("within_10kb", "Within 10 kb"),
        ("overlap_or_adjacent", "Overlap / adjacent"),
    ]
    values = {row.get("evidence_level", ""): _float_or_none(row.get("count", "")) or 0.0 for row in rows}
    plot_rows = [{"label": label, "count": f"{values.get(key, 0):.0f}"} for key, label in ladder]
    _write_bar_svg(path, plot_rows, title, "label", "count", "Feature-pair evidence count")


def _write_context_ladder_png(path: Path, rows: list[dict[str, str]]) -> None:
    _write_bar_png(path, rows, "count")


def _write_contig_neighborhood_svg(path: Path, rows: list[dict[str, str]], title: str) -> None:
    width = 1040
    height = max(220, 110 + 38 * max(len(rows), 1))
    intervals = []
    for row in rows:
        start = _as_int(row.get("start", ""))
        end = _as_int(row.get("end", ""))
        if start is None or end is None:
            continue
        if start > end:
            start, end = end, start
        intervals.append((start, end, row))
    if not intervals:
        Path(path).write_text(
            f"""<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='180' viewBox='0 0 {width} 180'>
<rect width='100%' height='100%' fill='#f8fafc'/>
<text x='20' y='32' font-family='Arial' font-size='20' font-weight='700' fill='#102a43'>{html.escape(title)}</text>
<text x='20' y='76' font-family='Arial' font-size='13' fill='#52606d'>No coordinate-complete feature rows were available for a neighborhood diagram.</text>
</svg>
""",
            encoding="utf-8",
        )
        return
    min_pos = min(start for start, _, _ in intervals)
    max_pos = max(end for _, end, _ in intervals)
    if min_pos == max_pos:
        max_pos = min_pos + 1
    left, right = 110, 950
    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='#f8fafc'/>",
        f"<text x='20' y='32' font-family='Arial' font-size='20' font-weight='700' fill='#102a43'>{html.escape(title)}</text>",
        f"<line x1='{left}' y1='78' x2='{right}' y2='78' stroke='#94a3b8' stroke-width='3'/>",
        f"<text x='{left}' y='100' font-family='Arial' font-size='11' fill='#52606d'>{min_pos} bp</text>",
        f"<text x='{right - 70}' y='100' font-family='Arial' font-size='11' fill='#52606d'>{max_pos} bp</text>",
    ]
    for idx, (start, end, row) in enumerate(sorted(intervals, key=lambda item: (item[0], item[1]))):
        x1 = left + (start - min_pos) / (max_pos - min_pos) * (right - left)
        x2 = left + (end - min_pos) / (max_pos - min_pos) * (right - left)
        if x2 - x1 < 8:
            x2 = x1 + 8
        y = 120 + idx * 34
        color = _database_color(row.get("database", ""))
        label = f"{row.get('database', '')}:{row.get('feature_id', '')}"
        short = label if len(label) <= 44 else label[:41] + "..."
        tip = f"{label}; {start}-{end}; identity={row.get('identity', '')}; coverage={row.get('coverage', '')}"
        parts.append(f"<rect x='{x1:.1f}' y='{y}' width='{x2 - x1:.1f}' height='18' rx='2' fill='{color}'><title>{html.escape(tip)}</title></rect>")
        parts.append(f"<text x='20' y='{y + 14}' font-family='Arial' font-size='11' fill='#1f2933'>{html.escape(short)}</text>")
    parts.append("</svg>\n")
    path.write_text("".join(parts), encoding="utf-8")


def _write_contig_neighborhood_png(path: Path, rows: list[dict[str, str]]) -> None:
    width = 1040
    height = max(220, 110 + 38 * max(len(rows), 1))
    pixels = [bytearray([248, 250, 252] * width) for _ in range(height)]
    intervals = []
    for row in rows:
        start = _as_int(row.get("start", ""))
        end = _as_int(row.get("end", ""))
        if start is None or end is None:
            continue
        if start > end:
            start, end = end, start
        intervals.append((start, end, row))
    if not intervals:
        _write_png(path, width, height, pixels)
        return
    min_pos = min(start for start, _, _ in intervals)
    max_pos = max(end for _, end, _ in intervals)
    if min_pos == max_pos:
        max_pos = min_pos + 1
    left, right = 110, 950

    def rect(x0: int, y0: int, x1: int, y1: int, color: tuple[int, int, int]) -> None:
        for y in range(max(0, y0), min(height, y1)):
            row_pixels = pixels[y]
            for x in range(max(0, x0), min(width, x1)):
                idx = x * 3
                row_pixels[idx:idx + 3] = bytes(color)

    color_map = {
        "amr": (220, 38, 38),
        "amrfinderplus": (239, 68, 68),
        "vfdb": (22, 163, 74),
        "plasmidfinder": (124, 58, 237),
        "integronfinder": (234, 88, 12),
        "isfinder": (202, 138, 4),
        "mobileelementfinder": (202, 138, 4),
    }
    rect(left, 76, right, 80, (148, 163, 184))
    for idx, (start, end, row) in enumerate(sorted(intervals, key=lambda item: (item[0], item[1]))):
        x1 = int(left + (start - min_pos) / (max_pos - min_pos) * (right - left))
        x2 = int(left + (end - min_pos) / (max_pos - min_pos) * (right - left))
        if x2 - x1 < 8:
            x2 = x1 + 8
        y = 120 + idx * 34
        rect(x1, y, x2, y + 18, color_map.get(row.get("database", ""), (100, 116, 139)))
    _write_png(path, width, height, pixels)

def _extract_year(value: str) -> int | None:
    match = re.search(r"(19|20)\d{2}", str(value or ""))
    if not match:
        return None
    year = int(match.group(0))
    if 1900 <= year <= 2100:
        return year
    return None


def _pearson_correlation(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 2 or len(xs) != len(ys):
        return None
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    numerator = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    denom_x = math.sqrt(sum((x - mean_x) ** 2 for x in xs))
    denom_y = math.sqrt(sum((y - mean_y) ** 2 for y in ys))
    if denom_x == 0 or denom_y == 0:
        return 0.0
    return numerator / (denom_x * denom_y)


def _temporal_trend_label(correlation: float | None, change_percent: float, years_observed: int) -> str:
    if years_observed < 3 or correlation is None:
        return "insufficient_temporal_support"
    if correlation >= 0.6 and change_percent >= 20:
        return "strong_increasing"
    if correlation >= 0.3 and change_percent > 0:
        return "increasing"
    if correlation <= -0.6 and change_percent <= -20:
        return "strong_decreasing"
    if correlation <= -0.3 and change_percent < 0:
        return "decreasing"
    return "stable"


def _temporal_support_label(years_observed: int, min_year_genomes: int, total_positive: int) -> str:
    if years_observed < 3 or min_year_genomes < 3 or total_positive < 3:
        return "low_support"
    if years_observed >= 5 and min_year_genomes >= 5 and total_positive >= 10:
        return "high_support"
    return "moderate_support"


def _temporal_pattern_label(prevalence_values: list[tuple[int, float, int, int]]) -> str:
    if not prevalence_values:
        return "insufficient_data"
    positives = [positive for _, _, _, positive in prevalence_values]
    prevalence = [value for _, value, _, _ in prevalence_values]
    if len(prevalence_values) >= 3:
        midpoint = max(1, len(prevalence_values) // 2)
        early_positive = sum(positives[:midpoint])
        late_positive = sum(positives[midpoint:])
        early_prevalence = max(prevalence[:midpoint]) if prevalence[:midpoint] else 0.0
        late_prevalence = max(prevalence[midpoint:]) if prevalence[midpoint:] else 0.0
        if early_positive == 0 and late_positive > 0:
            return "newly_detected"
        if early_positive > 0 and late_positive == 0:
            return "disappearing"
        if min(positives) > 0:
            return "persistent"
        if late_prevalence >= early_prevalence + 20:
            return "newly_detected"
        if early_prevalence >= late_prevalence + 20:
            return "disappearing"
    if sum(1 for positive in positives if positive > 0) <= max(1, len(positives) // 3):
        return "sporadic"
    return "persistent"


def write_important_qc_outputs(sample_dir: Path, out_dir: Path, important_dir: Path) -> dict[str, str]:
    key_tables = important_dir / "key_tables"
    figures = important_dir / "figures"
    key_tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    dataset_rows = read_table(sample_dir / "basic" / "enriched_genome_dataset.csv")
    total = len(dataset_rows)
    qc_pass = sum(1 for row in dataset_rows if row.get("qc_pass") == "true")
    qc_fail = sum(1 for row in dataset_rows if row.get("qc_pass") == "false" and row.get("qc_status"))
    qc_unknown = total - qc_pass - qc_fail

    def enabled_from_column(column: str) -> bool:
        return any(row.get(column) for row in dataset_rows)

    step_rows = [
        {"step_order": "1", "qc_step": "Genome download", "tool": "FetchM2", "enabled": "yes", "input_genomes": str(total), "pass": str(total), "warning": "0", "fail": "0", "skipped": "0", "output_genomes": str(total), "main_threshold": "downloaded FASTA present", "status": "PASS" if total else "WARNING_EMPTY", "notes": "Post-download metadata rows present in enriched dataset."},
        {"step_order": "2", "qc_step": "Sequence basic QC", "tool": "internal", "enabled": "yes", "input_genomes": str(total), "pass": str(qc_pass), "warning": str(qc_unknown), "fail": str(qc_fail), "skipped": "0", "output_genomes": str(qc_pass), "main_threshold": "FASTA statistics available", "status": "PASS" if qc_fail == 0 else "WARNING", "notes": "Combined QC status summarized per genome."},
        {"step_order": "3", "qc_step": "Assembly metrics", "tool": "QUAST", "enabled": "yes" if enabled_from_column("quast_status") else "no", "input_genomes": str(total), "pass": str(total if enabled_from_column("quast_status") else 0), "warning": "0", "fail": "0", "skipped": "0" if enabled_from_column("quast_status") else str(total), "output_genomes": str(total), "main_threshold": "assembly structure metrics", "status": "PASS" if enabled_from_column("quast_status") else "SKIPPED", "notes": "Enabled when QUAST output exists."},
        {"step_order": "4", "qc_step": "Completeness/contamination", "tool": "CheckM2", "enabled": "yes" if enabled_from_column("checkm2_completeness") else "no", "input_genomes": str(total), "pass": str(total if enabled_from_column("checkm2_completeness") else 0), "warning": "0", "fail": "0", "skipped": "0" if enabled_from_column("checkm2_completeness") else str(total), "output_genomes": str(total), "main_threshold": "completeness/contamination", "status": "PASS" if enabled_from_column("checkm2_completeness") else "SKIPPED", "notes": "Enabled when CheckM2 output exists."},
        {"step_order": "5", "qc_step": "ANI relatedness", "tool": "FastANI/skani", "enabled": "yes" if enabled_from_column("ani_status") else "no", "input_genomes": str(total), "pass": str(total if enabled_from_column("ani_status") else 0), "warning": "0", "fail": "0", "skipped": "0" if enabled_from_column("ani_status") else str(total), "output_genomes": str(total), "main_threshold": "species/cluster context", "status": "PASS" if enabled_from_column("ani_status") else "SKIPPED", "notes": "Skipped for one-genome or disabled runs."},
        {"step_order": "6", "qc_step": "Mash relatedness", "tool": "Mash", "enabled": "yes" if enabled_from_column("mash_status") else "no", "input_genomes": str(total), "pass": str(total if enabled_from_column("mash_status") else 0), "warning": "0", "fail": "0", "skipped": "0" if enabled_from_column("mash_status") else str(total), "output_genomes": str(total), "main_threshold": "sketch distance context", "status": "PASS" if enabled_from_column("mash_status") else "SKIPPED", "notes": "Skipped for one-genome or disabled runs."},
        {"step_order": "7", "qc_step": "Combined QC decision", "tool": "PanResistome", "enabled": "yes", "input_genomes": str(total), "pass": str(qc_pass), "warning": str(qc_unknown), "fail": str(qc_fail), "skipped": "0", "output_genomes": str(qc_pass), "main_threshold": "combined rules", "status": "PASS" if qc_fail == 0 else "WARNING", "notes": "Genomes passing combined QC are sent to annotation when filtering is enabled."},
    ]
    step_fields = ["step_order", "qc_step", "tool", "enabled", "input_genomes", "pass", "warning", "fail", "skipped", "output_genomes", "main_threshold", "status", "notes"]
    step_path = key_tables / "qc_step_summary.tsv"
    write_rows(step_path, step_rows, step_fields)

    qc_by_genome_fields = [
        "assembly_accession", "sample_id", "organism_name", "qc_status", "qc_pass", "qc_fail_reasons",
        "quast_status", "checkm2_completeness", "checkm2_contamination", "ani_status", "mash_status",
        "genome_size", "contig_count", "n50", "gc_percent", "ani_cluster", "mash_cluster",
    ]
    qc_by_genome = [{field: row.get(field, "") for field in qc_by_genome_fields} for row in dataset_rows]
    qc_by_genome_path = key_tables / "qc_by_genome.tsv"
    write_rows(qc_by_genome_path, qc_by_genome, qc_by_genome_fields)

    funnel_rows = [
        {"step": row["qc_step"], "genomes": row["output_genomes"], "status": row["status"]}
        for row in step_rows
    ]
    funnel_path = figures / "qc_funnel.data.tsv"
    write_rows(funnel_path, funnel_rows, ["step", "genomes", "status"])
    _write_bar_svg(figures / "qc_funnel.svg", funnel_rows, "QC Funnel", "step", "genomes", "Genomes")
    _write_bar_png(figures / "qc_funnel.png", funnel_rows, "genomes")

    status_rows = []
    for row in step_rows:
        for status_field in ["pass", "warning", "fail", "skipped"]:
            status_rows.append({"qc_step": row["qc_step"], "status": status_field.upper(), "count": row[status_field]})
    status_path = figures / "qc_status_overview.data.tsv"
    write_rows(status_path, status_rows, ["qc_step", "status", "count"])
    compact_status = [{"label": f"{row['qc_step']} {row['status']}", "count": row["count"]} for row in status_rows if row["count"] != "0"]
    _write_bar_svg(figures / "qc_status_overview.svg", compact_status, "QC Status Overview", "label", "count", "Genomes")
    _write_bar_png(figures / "qc_status_overview.png", compact_status, "count")

    outputs = {
        "important_qc_step_summary": str(step_path),
        "important_qc_by_genome": str(qc_by_genome_path),
        "important_qc_funnel_svg": str(figures / "qc_funnel.svg"),
        "important_qc_funnel_png": str(figures / "qc_funnel.png"),
        "important_qc_funnel_data": str(funnel_path),
        "important_qc_status_svg": str(figures / "qc_status_overview.svg"),
        "important_qc_status_png": str(figures / "qc_status_overview.png"),
        "important_qc_status_data": str(status_path),
    }

    checkm2_rows = [
        {
            "assembly_accession": row.get("assembly_accession", ""),
            "completeness": row.get("checkm2_completeness", ""),
            "contamination": row.get("checkm2_contamination", ""),
        }
        for row in dataset_rows
        if row.get("checkm2_completeness") or row.get("checkm2_contamination")
    ]
    if checkm2_rows:
        checkm2_path = figures / "checkm2_completeness_contamination.data.tsv"
        write_rows(checkm2_path, checkm2_rows, ["assembly_accession", "completeness", "contamination"])
        outputs["important_checkm2_data"] = str(checkm2_path)
    return outputs


def write_important_prevalence_outputs(sample_dir: Path, out_dir: Path, important_dir: Path, top_n: int = 20) -> dict[str, str]:
    key_tables = important_dir / "key_tables"
    figures = important_dir / "figures"
    key_tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    metadata_rows = read_table(sample_dir / "basic" / "enriched_genome_dataset.csv")
    sample_count = len(metadata_rows)
    features = [row for row in read_table(out_dir / "features" / "all_features.tsv") if row.get("presence", "1") != "0"]
    by_feature: dict[tuple[str, str], set[str]] = defaultdict(set)
    rows_by_feature: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in features:
        key = (row.get("database", ""), row.get("feature_id", ""))
        sample = row.get("assembly_accession", "") or row.get("sample_id", "")
        if key[0] and key[1] and sample:
            by_feature[key].add(sample)
            rows_by_feature[key].append(row)

    summary_rows = []
    for (database, feature_id), samples in sorted(by_feature.items()):
        feature_rows = rows_by_feature[(database, feature_id)]
        category = first_value(feature_rows[0], ["feature_category"], "")
        prevalence = len(samples) / sample_count if sample_count else 0.0
        summary_rows.append({
            "database": database,
            "feature_id": feature_id,
            "feature_category": category,
            "feature_rows": str(len(feature_rows)),
            "positive_genomes": str(len(samples)),
            "sample_count": str(sample_count),
            "prevalence": f"{prevalence:.4f}",
            "prevalence_percent": f"{prevalence * 100:.1f}",
        })
    summary_rows.sort(key=lambda row: (row["database"], -int(row["positive_genomes"]), row["feature_id"]))
    summary_path = key_tables / "feature_prevalence_summary.tsv"
    summary_fields = ["database", "feature_id", "feature_category", "feature_rows", "positive_genomes", "sample_count", "prevalence", "prevalence_percent"]
    write_rows(summary_path, summary_rows, summary_fields)

    outputs = {"important_feature_prevalence_summary": str(summary_path)}
    for database in sorted({row["database"] for row in summary_rows}):
        db_rows = [row for row in summary_rows if row["database"] == database][:top_n]
        if not db_rows:
            continue
        data_path = figures / f"prevalence_{database}_top20.data.tsv"
        svg_path = figures / f"prevalence_{database}_top20.svg"
        png_path = figures / f"prevalence_{database}_top20.png"
        write_rows(data_path, db_rows, summary_fields)
        _write_bar_svg(svg_path, db_rows, f"{database} Prevalence Top {top_n}", "feature_id", "prevalence_percent", "Prevalence (%)")
        _write_bar_png(png_path, db_rows, "prevalence_percent")
        outputs[f"important_prevalence_{database}_data"] = str(data_path)
        outputs[f"important_prevalence_{database}_svg"] = str(svg_path)
        outputs[f"important_prevalence_{database}_png"] = str(png_path)
    return outputs


def _summary_stats(values: list[float]) -> dict[str, str]:
    if not values:
        return {"median": "", "min": "", "max": "", "iqr": ""}
    values = sorted(values)
    mid = len(values) // 2
    median = values[mid] if len(values) % 2 else (values[mid - 1] + values[mid]) / 2
    q1 = values[len(values) // 4]
    q3 = values[(len(values) * 3) // 4]
    return {"median": f"{median:.2f}", "min": f"{values[0]:.2f}", "max": f"{values[-1]:.2f}", "iqr": f"{q3 - q1:.2f}"}


def write_important_variation_outputs(sample_dir: Path, out_dir: Path, important_dir: Path, top_n: int = 20) -> dict[str, str]:
    key_tables = important_dir / "key_tables"
    figures = important_dir / "figures"
    key_tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    features = [row for row in read_table(out_dir / "features" / "all_features.tsv") if row.get("presence", "1") != "0"]
    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    hit_rows = []
    for row in features:
        database = row.get("database", "")
        feature_id = row.get("feature_id", "")
        if not database or not feature_id:
            continue
        grouped[(database, feature_id)].append(row)
        hit_rows.append({
            "database": database,
            "feature_id": feature_id,
            "assembly_accession": row.get("assembly_accession", ""),
            "sample_id": row.get("sample_id", ""),
            "identity": row.get("identity", ""),
            "coverage": row.get("coverage", ""),
            "contig": row.get("contig", ""),
            "start": row.get("start", ""),
            "end": row.get("end", ""),
            "tool": row.get("tool", ""),
            "source_file": row.get("source_file", ""),
        })

    summary_rows = []
    for (database, feature_id), rows in sorted(grouped.items()):
        identities = [_float_or_none(row.get("identity", "")) for row in rows]
        coverages = [_float_or_none(row.get("coverage", "")) for row in rows]
        identities = [value for value in identities if value is not None]
        coverages = [value for value in coverages if value is not None]
        identity_stats = _summary_stats(identities)
        coverage_stats = _summary_stats(coverages)
        samples = {row.get("assembly_accession", "") or row.get("sample_id", "") for row in rows}
        low_identity = sum(1 for value in identities if value < 90)
        low_coverage = sum(1 for value in coverages if value < 80)
        warnings = []
        if low_identity:
            warnings.append("low_identity")
        if low_coverage:
            warnings.append("low_coverage")
        iqr_identity = _float_or_none(identity_stats["iqr"]) or 0.0
        label = "high_variation" if iqr_identity >= 10 else ("moderate_variation" if iqr_identity >= 3 else "low_variation")
        summary_rows.append({
            "database": database,
            "feature_id": feature_id,
            "total_hits": str(len(rows)),
            "positive_genomes": str(len(samples)),
            "median_identity": identity_stats["median"],
            "min_identity": identity_stats["min"],
            "max_identity": identity_stats["max"],
            "iqr_identity": identity_stats["iqr"],
            "median_coverage": coverage_stats["median"],
            "min_coverage": coverage_stats["min"],
            "max_coverage": coverage_stats["max"],
            "iqr_coverage": coverage_stats["iqr"],
            "low_identity_hits": str(low_identity),
            "low_coverage_hits": str(low_coverage),
            "variation_label": label,
            "warning_flags": ";".join(warnings),
        })
    summary_rows.sort(key=lambda row: (row["database"], -(_float_or_none(row["iqr_identity"]) or 0.0), row["feature_id"]))
    summary_fields = [
        "database", "feature_id", "total_hits", "positive_genomes",
        "median_identity", "min_identity", "max_identity", "iqr_identity",
        "median_coverage", "min_coverage", "max_coverage", "iqr_coverage",
        "low_identity_hits", "low_coverage_hits", "variation_label", "warning_flags",
    ]
    summary_path = key_tables / "feature_variation_summary.tsv"
    hits_path = key_tables / "feature_variation_hits.tsv"
    write_rows(summary_path, summary_rows, summary_fields)
    write_rows(hits_path, hit_rows, ["database", "feature_id", "assembly_accession", "sample_id", "identity", "coverage", "contig", "start", "end", "tool", "source_file"])

    outputs = {
        "important_feature_variation_summary": str(summary_path),
        "important_feature_variation_hits": str(hits_path),
    }
    for database in sorted({row["database"] for row in summary_rows}):
        db_rows = [row for row in summary_rows if row["database"] == database][:top_n]
        if not db_rows:
            continue
        data_path = figures / f"variation_identity_{database}_top20.data.tsv"
        svg_path = figures / f"variation_identity_{database}_top20.svg"
        png_path = figures / f"variation_identity_{database}_top20.png"
        write_rows(data_path, db_rows, summary_fields)
        _write_bar_svg(svg_path, db_rows, f"{database} Identity Variation Top {top_n}", "feature_id", "iqr_identity", "Identity IQR")
        _write_bar_png(png_path, db_rows, "iqr_identity")
        outputs[f"important_variation_{database}_data"] = str(data_path)
        outputs[f"important_variation_{database}_svg"] = str(svg_path)
        outputs[f"important_variation_{database}_png"] = str(png_path)
    return outputs


def write_important_temporal_outputs(sample_dir: Path, out_dir: Path, important_dir: Path, top_n: int = 20) -> dict[str, str]:
    key_tables = important_dir / "key_tables"
    figures = important_dir / "figures"
    key_tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    metadata_rows = read_table(sample_dir / "basic" / "enriched_genome_dataset.csv")
    features = [row for row in read_table(out_dir / "features" / "all_features.tsv") if row.get("presence", "1") != "0"]
    metadata_by_sample = {row.get("assembly_accession", ""): row for row in metadata_rows if row.get("assembly_accession")}
    sample_year: dict[str, int] = {}
    missing_year = 0
    for row in metadata_rows:
        sample = row.get("assembly_accession", "")
        year = _extract_year(row.get("collection_year", ""))
        if sample and year is not None:
            sample_year[sample] = year
        elif sample:
            missing_year += 1

    samples_by_year: dict[int, set[str]] = defaultdict(set)
    for sample, year in sample_year.items():
        samples_by_year[year].add(sample)
    years = sorted(samples_by_year)
    missing_year_fraction = missing_year / len(metadata_rows) if metadata_rows else 0.0

    feature_samples: dict[tuple[str, str], set[str]] = defaultdict(set)
    sample_database_features: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in features:
        database = row.get("database", "")
        feature_id = row.get("feature_id", "")
        sample = row.get("assembly_accession", "") or row.get("sample_id", "")
        if not database or not feature_id or not sample:
            continue
        feature_samples[(database, feature_id)].add(sample)
        sample_database_features[(sample, database)].add(feature_id)

    databases = sorted({database for database, _ in feature_samples})
    burden_rows = []
    for database in databases:
        for year in years:
            year_samples = samples_by_year[year]
            counts = [len(sample_database_features.get((sample, database), set())) for sample in year_samples]
            positive = sum(1 for count in counts if count > 0)
            total_features = sum(counts)
            mean_burden = total_features / len(year_samples) if year_samples else 0.0
            warnings = []
            if len(year_samples) < 3:
                warnings.append("small_year_group")
            if missing_year_fraction > 0:
                warnings.append("missing_collection_year")
            burden_rows.append({
                "database": database,
                "collection_year": str(year),
                "total_genomes": str(len(year_samples)),
                "positive_genomes": str(positive),
                "total_feature_count": str(total_features),
                "mean_feature_count_per_genome": f"{mean_burden:.4f}",
                "warning_flags": ";".join(warnings),
            })

    prevalence_rows = []
    trend_rows = []
    for (database, feature_id), present_samples in sorted(feature_samples.items()):
        prevalence_values = []
        total_positive = 0
        min_year_genomes = min((len(samples_by_year[year]) for year in years), default=0)
        max_year_genomes = max((len(samples_by_year[year]) for year in years), default=0)
        for year in years:
            year_samples = samples_by_year[year]
            positive = len(present_samples & year_samples)
            total = len(year_samples)
            total_positive += positive
            prevalence = positive / total if total else 0.0
            warnings = []
            if total < 3:
                warnings.append("small_year_group")
            if missing_year_fraction > 0:
                warnings.append("missing_collection_year")
            if positive < 3:
                warnings.append("few_positive_genomes")
            prevalence_rows.append({
                "database": database,
                "feature_id": feature_id,
                "collection_year": str(year),
                "total_genomes": str(total),
                "positive_genomes": str(positive),
                "prevalence": f"{prevalence:.4f}",
                "prevalence_percent": f"{prevalence * 100:.1f}",
                "warning_flags": ";".join(warnings),
            })
            prevalence_values.append((year, prevalence * 100, total, positive))

        year_numbers = [float(item[0]) for item in prevalence_values]
        prevalence_percents = [item[1] for item in prevalence_values]
        correlation = _pearson_correlation(year_numbers, prevalence_percents)
        first_year = prevalence_values[0] if prevalence_values else (None, 0.0, 0, 0)
        last_year = prevalence_values[-1] if prevalence_values else (None, 0.0, 0, 0)
        change = last_year[1] - first_year[1]
        years_observed = len(prevalence_values)
        trend_label = _temporal_trend_label(correlation, change, years_observed)
        support_label = _temporal_support_label(years_observed, min_year_genomes, total_positive)
        pattern_label = _temporal_pattern_label(prevalence_values)

        support_samples = [sample for sample in present_samples if sample in sample_year]
        support_bioprojects = Counter(metadata_by_sample.get(sample, {}).get("bioproject", "") for sample in support_samples)
        support_bioprojects.pop("", None)
        largest_bioproject, largest_bioproject_count = support_bioprojects.most_common(1)[0] if support_bioprojects else ("", 0)
        largest_bioproject_fraction = largest_bioproject_count / len(support_samples) if support_samples else 0.0

        warnings = ["exploratory_only"]
        if years_observed < 3:
            warnings.append("insufficient_temporal_support")
        if min_year_genomes < 3:
            warnings.append("small_year_group")
        if total_positive < 3:
            warnings.append("few_positive_genomes")
        if missing_year_fraction > 0:
            warnings.append("missing_collection_year")
        if max_year_genomes and min_year_genomes and max_year_genomes / max(min_year_genomes, 1) >= 5:
            warnings.append("uneven_sampling_by_year")
        if largest_bioproject_fraction >= 0.8 and len(support_samples) >= 3:
            warnings.append("single_bioproject_dominance")
        country_year_dominated = False
        lineage_year_dominated = False
        for year in years:
            if len(samples_by_year[year]) < 3:
                continue
            country_counts = Counter(metadata_by_sample.get(sample, {}).get("country", "") for sample in samples_by_year[year])
            country_counts.pop("", None)
            lineage_counts = Counter(metadata_by_sample.get(sample, {}).get("mlst_ST", "") for sample in samples_by_year[year])
            lineage_counts.pop("", None)
            if country_counts and country_counts.most_common(1)[0][1] / len(samples_by_year[year]) >= 0.8:
                country_year_dominated = True
            if lineage_counts and lineage_counts.most_common(1)[0][1] / len(samples_by_year[year]) >= 0.8:
                lineage_year_dominated = True
        if country_year_dominated:
            warnings.append("country_year_confounding")
        if lineage_year_dominated:
            warnings.append("lineage_year_confounding")

        trend_rows.append({
            "database": database,
            "feature_id": feature_id,
            "years_observed": str(years_observed),
            "first_year": str(first_year[0] or ""),
            "last_year": str(last_year[0] or ""),
            "first_year_prevalence_percent": f"{first_year[1]:.1f}",
            "last_year_prevalence_percent": f"{last_year[1]:.1f}",
            "change_percent_points": f"{change:.1f}",
            "correlation": "" if correlation is None else f"{correlation:.4f}",
            "trend_label": trend_label,
            "support_label": support_label,
            "temporal_pattern_label": pattern_label,
            "total_positive_genomes": str(total_positive),
            "min_year_genomes": str(min_year_genomes),
            "max_year_genomes": str(max_year_genomes),
            "largest_bioproject": largest_bioproject,
            "largest_bioproject_fraction": f"{largest_bioproject_fraction:.4f}",
            "warning_flags": ";".join(dict.fromkeys(warnings)),
        })

    burden_fields = ["database", "collection_year", "total_genomes", "positive_genomes", "total_feature_count", "mean_feature_count_per_genome", "warning_flags"]
    prevalence_fields = ["database", "feature_id", "collection_year", "total_genomes", "positive_genomes", "prevalence", "prevalence_percent", "warning_flags"]
    trend_fields = [
        "database", "feature_id", "years_observed", "first_year", "last_year",
        "first_year_prevalence_percent", "last_year_prevalence_percent", "change_percent_points",
        "correlation", "trend_label", "support_label", "temporal_pattern_label", "total_positive_genomes",
        "min_year_genomes", "max_year_genomes", "largest_bioproject",
        "largest_bioproject_fraction", "warning_flags",
    ]
    burden_path = key_tables / "temporal_database_burden.tsv"
    prevalence_path = key_tables / "temporal_feature_prevalence.tsv"
    trend_path = key_tables / "temporal_trend_summary.tsv"
    increasing_path = key_tables / "temporal_increasing_features.tsv"
    decreasing_path = key_tables / "temporal_decreasing_features.tsv"
    write_rows(burden_path, burden_rows, burden_fields)
    write_rows(prevalence_path, prevalence_rows, prevalence_fields)
    write_rows(trend_path, trend_rows, trend_fields)

    increasing = [
        row for row in trend_rows
        if row["trend_label"] in {"increasing", "strong_increasing"}
    ]
    increasing.sort(key=lambda row: (-(_float_or_none(row["change_percent_points"]) or 0.0), row["database"], row["feature_id"]))
    decreasing = [
        row for row in trend_rows
        if row["trend_label"] in {"decreasing", "strong_decreasing"}
    ]
    decreasing.sort(key=lambda row: (_float_or_none(row["change_percent_points"]) or 0.0, row["database"], row["feature_id"]))
    write_rows(increasing_path, increasing, trend_fields)
    write_rows(decreasing_path, decreasing, trend_fields)

    outputs = {
        "important_temporal_database_burden": str(burden_path),
        "important_temporal_feature_prevalence": str(prevalence_path),
        "important_temporal_trend_summary": str(trend_path),
        "important_temporal_increasing_features": str(increasing_path),
        "important_temporal_decreasing_features": str(decreasing_path),
    }

    latest_burden = []
    for database in databases:
        db_rows = [row for row in burden_rows if row["database"] == database]
        if not db_rows:
            continue
        latest = max(db_rows, key=lambda row: int(row["collection_year"]))
        latest_burden.append({"database": database, "mean_feature_count_per_genome": latest["mean_feature_count_per_genome"]})
    latest_burden.sort(key=lambda row: -(_float_or_none(row["mean_feature_count_per_genome"]) or 0.0))
    burden_data = figures / "temporal_database_burden_top20.data.tsv"
    burden_svg = figures / "temporal_database_burden_top20.svg"
    burden_png = figures / "temporal_database_burden_top20.png"
    write_rows(burden_data, latest_burden[:top_n], ["database", "mean_feature_count_per_genome"])
    _write_bar_svg(burden_svg, latest_burden[:top_n], "Temporal Database Burden", "database", "mean_feature_count_per_genome", "Mean features/genome in latest year")
    _write_bar_png(burden_png, latest_burden[:top_n], "mean_feature_count_per_genome")
    outputs.update({
        "important_temporal_database_burden_data": str(burden_data),
        "important_temporal_database_burden_svg": str(burden_svg),
        "important_temporal_database_burden_png": str(burden_png),
    })

    for label, rows, filename in [
        ("Top Increasing Temporal Features", increasing[:top_n], "temporal_top_increasing_features"),
        ("Top Decreasing Temporal Features", decreasing[:top_n], "temporal_top_decreasing_features"),
    ]:
        data_path = figures / f"{filename}.data.tsv"
        svg_path = figures / f"{filename}.svg"
        png_path = figures / f"{filename}.png"
        figure_rows = [
            {
                **row,
                "feature_label": f"{row['database']}:{row['feature_id']}",
                "display_change_percent_points": f"{abs(_float_or_none(row['change_percent_points']) or 0.0):.1f}",
            }
            for row in rows
        ]
        write_rows(data_path, figure_rows, trend_fields + ["feature_label", "display_change_percent_points"])
        _write_bar_svg(svg_path, figure_rows, label, "feature_label", "display_change_percent_points", "Absolute change in prevalence percentage points")
        _write_bar_png(png_path, figure_rows, "display_change_percent_points")
        outputs[f"important_{filename}_data"] = str(data_path)
        outputs[f"important_{filename}_svg"] = str(svg_path)
        outputs[f"important_{filename}_png"] = str(png_path)

    heatmap_features = increasing[:top_n] + decreasing[:top_n]
    heatmap_keys = {(row["database"], row["feature_id"]) for row in heatmap_features[: top_n * 2]}
    heatmap_rows = [
        {**row, "feature_label": f"{row['database']}:{row['feature_id']}"}
        for row in prevalence_rows
        if (row["database"], row["feature_id"]) in heatmap_keys
    ]
    heatmap_data = figures / "temporal_feature_heatmap_top40.data.tsv"
    heatmap_svg = figures / "temporal_feature_heatmap_top40.svg"
    heatmap_png = figures / "temporal_feature_heatmap_top40.png"
    write_rows(heatmap_data, heatmap_rows, prevalence_fields + ["feature_label"])
    _write_heatmap_svg(heatmap_svg, heatmap_rows, "Temporal Feature Prevalence Heatmap", "feature_label", "collection_year", "prevalence_percent")
    _write_heatmap_png(heatmap_png, heatmap_rows, "feature_label", "collection_year", "prevalence_percent")
    outputs.update({
        "important_temporal_heatmap_data": str(heatmap_data),
        "important_temporal_heatmap_svg": str(heatmap_svg),
        "important_temporal_heatmap_png": str(heatmap_png),
    })

    selected_feature = (increasing or decreasing or sorted(trend_rows, key=lambda row: -abs(_float_or_none(row.get("change_percent_points", "")) or 0.0)))[:1]
    selected_rows = []
    selected_title = "Selected Feature Temporal Prevalence"
    if selected_feature:
        selected = selected_feature[0]
        selected_title = f"{selected['database']}:{selected['feature_id']} Temporal Prevalence"
        selected_rows = [
            row for row in prevalence_rows
            if row["database"] == selected["database"] and row["feature_id"] == selected["feature_id"]
        ]
    selected_data = figures / "temporal_selected_feature_prevalence.data.tsv"
    selected_svg = figures / "temporal_selected_feature_prevalence.svg"
    selected_png = figures / "temporal_selected_feature_prevalence.png"
    write_rows(selected_data, selected_rows, prevalence_fields)
    _write_temporal_series_svg(selected_svg, selected_rows, selected_title)
    _write_temporal_series_png(selected_png, selected_rows)

    slope_rows = (increasing[:top_n] + decreasing[:top_n])[: top_n * 2]
    if not slope_rows:
        slope_rows = sorted(trend_rows, key=lambda row: -abs(_float_or_none(row.get("change_percent_points", "")) or 0.0))[:top_n]
    slope_data = figures / "temporal_slope_top40.data.tsv"
    slope_svg = figures / "temporal_slope_top40.svg"
    slope_png = figures / "temporal_slope_top40.png"
    write_rows(slope_data, slope_rows, trend_fields)
    _write_temporal_slope_svg(slope_svg, slope_rows, "Temporal Slope Plot")
    _write_temporal_slope_png(slope_png, slope_rows)

    temporal_html = figures / "temporal_trends.html"
    temporal_html.write_text(
        f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Temporal Trends</title>
<style>
body {{ font-family: Arial, sans-serif; color: #1f2933; margin: 1.5rem; }}
label {{ font-weight: 700; margin-right: 0.35rem; }}
select {{ margin: 0 1rem 0.75rem 0; padding: 0.35rem; max-width: 22rem; }}
.warning {{ background: #fff7ed; border-left: 4px solid #c2410c; padding: 0.75rem; margin: 0.75rem 0; }}
.grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(420px, 1fr)); gap: 1rem; }}
.panel {{ border: 1px solid #d9e2ec; border-radius: 6px; padding: 0.75rem; background: #f8fafc; }}
svg {{ max-width: 100%; height: auto; }}
table {{ border-collapse: collapse; width: 100%; font-size: 0.88rem; margin-top: 0.75rem; }}
th, td {{ border: 1px solid #d9e2ec; padding: 0.35rem; text-align: left; }}
th {{ background: #f0f4f8; }}
</style></head><body>
<h1>Temporal Trends</h1>
<p>Use the controls to inspect database burden, feature prevalence, and trend support over collection years.</p>
<div class="warning">Temporal trends are exploratory. Sampling year, BioProject composition, geography, lineage, and missing collection-year metadata can change the apparent pattern.</div>
<label for="database">Database</label><select id="database"></select>
<label for="trend">Trend</label><select id="trend"><option value="all">All</option><option value="increasing">Increasing</option><option value="decreasing">Decreasing</option><option value="stable">Stable</option></select>
<label for="support">Support</label><select id="support"><option value="all">All</option><option value="high_support">High</option><option value="moderate_support">Moderate</option><option value="low_support">Low</option></select>
<label for="feature">Feature</label><select id="feature"></select>
<div id="summary"></div>
<div class="grid"><div class="panel"><h2>Selected Feature Prevalence</h2><div id="linePlot"></div></div><div class="panel"><h2>First-to-Last Year Slope</h2><div id="slopePlot"></div></div></div>
<div class="panel"><h2>Filtered Trend Table</h2><div id="trendTable"></div></div>
<p><a href="temporal_selected_feature_prevalence.png">Download default line PNG</a> | <a href="temporal_selected_feature_prevalence.svg">Download default line SVG</a> | <a href="temporal_slope_top40.png">Download slope PNG</a> | <a href="temporal_slope_top40.svg">Download slope SVG</a> | <a href="temporal_trend_summary.data.tsv">Download trend data TSV</a></p>
<script>
const trends = {json.dumps(trend_rows)};
const prevalence = {json.dumps(prevalence_rows)};
const burden = {json.dumps(burden_rows)};
const dbSelect = document.getElementById('database');
const trendSelect = document.getElementById('trend');
const supportSelect = document.getElementById('support');
const featureSelect = document.getElementById('feature');
function num(value) {{ const n = Number(value); return Number.isFinite(n) ? n : 0; }}
function label(row) {{ return row.database + ':' + row.feature_id; }}
function trendGroup(value) {{
  if (value.includes('increasing')) return 'increasing';
  if (value.includes('decreasing')) return 'decreasing';
  if (value === 'stable') return 'stable';
  return 'all';
}}
function filteredTrends() {{
  return trends.filter(row => {{
    if (dbSelect.value !== 'all' && row.database !== dbSelect.value) return false;
    if (trendSelect.value !== 'all' && trendGroup(row.trend_label) !== trendSelect.value) return false;
    if (supportSelect.value !== 'all' && row.support_label !== supportSelect.value) return false;
    return true;
  }}).sort((a, b) => Math.abs(num(b.change_percent_points)) - Math.abs(num(a.change_percent_points)));
}}
function populateDatabases() {{
  const dbs = ['all', ...Array.from(new Set(trends.map(row => row.database))).sort()];
  dbSelect.innerHTML = dbs.map(db => `<option value="${{db}}">${{db === 'all' ? 'All' : db}}</option>`).join('');
  const amr = dbs.find(db => db.toLowerCase() === 'amr');
  if (amr) dbSelect.value = amr;
  trendSelect.value = 'increasing';
  supportSelect.value = 'all';
}}
function populateFeatures() {{
  const rows = filteredTrends();
  featureSelect.innerHTML = rows.map(row => `<option value="${{label(row)}}">${{label(row)}} (${{row.trend_label}}, Δ${{row.change_percent_points}})</option>`).join('');
  if (!rows.length) featureSelect.innerHTML = '<option value="">No matching features</option>';
}}
function renderLine() {{
  const selected = featureSelect.value;
  const rows = prevalence.filter(row => label(row) === selected).sort((a, b) => num(a.collection_year) - num(b.collection_year));
  const width = 760, height = 320, left = 58, top = 30, plotW = 650, plotH = 220;
  const years = rows.map(row => num(row.collection_year));
  const values = rows.map(row => num(row.prevalence_percent));
  const minYear = years.length ? Math.min(...years) : 0;
  const maxYear = years.length && Math.max(...years) !== minYear ? Math.max(...years) : minYear + 1;
  let svg = `<svg width="${{width}}" height="${{height}}" viewBox="0 0 ${{width}} ${{height}}" xmlns="http://www.w3.org/2000/svg"><rect width="100%" height="100%" fill="#f8fafc"/>`;
  svg += `<line x1="${{left}}" y1="${{top + plotH}}" x2="${{left + plotW}}" y2="${{top + plotH}}" stroke="#9fb3c8"/><line x1="${{left}}" y1="${{top}}" x2="${{left}}" y2="${{top + plotH}}" stroke="#9fb3c8"/>`;
  const points = rows.map(row => {{
    const x = left + (num(row.collection_year) - minYear) / (maxYear - minYear) * plotW;
    const y = top + plotH - num(row.prevalence_percent) / 100 * plotH;
    return [x, y, row];
  }});
  svg += `<polyline points="${{points.map(p => p[0].toFixed(1) + ',' + p[1].toFixed(1)).join(' ')}}" fill="none" stroke="#0f766e" stroke-width="3"/>`;
  for (const [x, y, row] of points) {{
    const tip = `${{row.collection_year}}: ${{row.prevalence_percent}}% (${{row.positive_genomes}}/${{row.total_genomes}})`;
    svg += `<circle cx="${{x.toFixed(1)}}" cy="${{y.toFixed(1)}}" r="5" fill="#0f766e"><title>${{tip}}</title></circle><text x="${{(x - 14).toFixed(1)}}" y="${{top + plotH + 20}}" font-size="11" fill="#52606d">${{row.collection_year}}</text>`;
  }}
  svg += `<text x="8" y="${{top + 15}}" font-size="12" fill="#52606d">Prevalence %</text></svg>`;
  document.getElementById('linePlot').innerHTML = svg;
  if (rows.length) {{
    const first = rows[0], last = rows[rows.length - 1];
    document.getElementById('summary').innerHTML = `<p><strong>${{selected}}</strong> changed from ${{first.prevalence_percent}}% (${{first.positive_genomes}}/${{first.total_genomes}}) in ${{first.collection_year}} to ${{last.prevalence_percent}}% (${{last.positive_genomes}}/${{last.total_genomes}}) in ${{last.collection_year}}.</p>`;
  }} else {{
    document.getElementById('summary').innerHTML = '<p>No matching temporal data.</p>';
  }}
}}
function renderSlope() {{
  const rows = filteredTrends().slice(0, 30);
  const width = 760, rowH = 24, top = 40, leftX = 250, rightX = 560;
  const height = Math.max(180, top + rows.length * rowH + 30);
  let svg = `<svg width="${{width}}" height="${{height}}" viewBox="0 0 ${{width}} ${{height}}" xmlns="http://www.w3.org/2000/svg"><rect width="100%" height="100%" fill="#f8fafc"/>`;
  svg += `<text x="${{leftX - 30}}" y="24" font-size="12" fill="#52606d">First year</text><text x="${{rightX - 20}}" y="24" font-size="12" fill="#52606d">Last year</text>`;
  rows.forEach((row, idx) => {{
    const yBase = top + idx * rowH + 12;
    const y1 = yBase - num(row.first_year_prevalence_percent) / 100 * 8;
    const y2 = yBase - num(row.last_year_prevalence_percent) / 100 * 8;
    const color = row.trend_label.includes('increasing') ? '#0f766e' : (row.trend_label.includes('decreasing') ? '#b91c1c' : '#64748b');
    const short = label(row).length > 30 ? label(row).slice(0, 27) + '...' : label(row);
    svg += `<text x="10" y="${{yBase + 4}}" font-size="11" fill="#1f2933">${{short}}</text><line x1="${{leftX}}" y1="${{y1.toFixed(1)}}" x2="${{rightX}}" y2="${{y2.toFixed(1)}}" stroke="${{color}}" stroke-width="2.4"><title>${{label(row)}}: ${{row.first_year_prevalence_percent}}% to ${{row.last_year_prevalence_percent}}%</title></line><circle cx="${{leftX}}" cy="${{y1.toFixed(1)}}" r="3.5" fill="${{color}}"/><circle cx="${{rightX}}" cy="${{y2.toFixed(1)}}" r="3.5" fill="${{color}}"/>`;
  }});
  svg += '</svg>';
  document.getElementById('slopePlot').innerHTML = svg;
}}
function renderTable() {{
  const rows = filteredTrends().slice(0, 50);
  const cols = ['database','feature_id','first_year','last_year','change_percent_points','correlation','trend_label','support_label','temporal_pattern_label','warning_flags'];
  let html = '<table><thead><tr>' + cols.map(c => `<th>${{c}}</th>`).join('') + '</tr></thead><tbody>';
  for (const row of rows) html += '<tr>' + cols.map(c => `<td>${{row[c] || ''}}</td>`).join('') + '</tr>';
  html += '</tbody></table>';
  document.getElementById('trendTable').innerHTML = html;
}}
function renderAll() {{ populateFeatures(); renderLine(); renderSlope(); renderTable(); }}
populateDatabases(); renderAll();
dbSelect.addEventListener('change', renderAll);
trendSelect.addEventListener('change', renderAll);
supportSelect.addEventListener('change', renderAll);
featureSelect.addEventListener('change', () => {{ renderLine(); }});
</script></body></html>
""",
        encoding="utf-8",
    )
    (figures / "temporal_trend_summary.data.tsv").write_text(trend_path.read_text(encoding="utf-8"), encoding="utf-8")
    outputs.update({
        "important_temporal_selected_feature_data": str(selected_data),
        "important_temporal_selected_feature_svg": str(selected_svg),
        "important_temporal_selected_feature_png": str(selected_png),
        "important_temporal_slope_data": str(slope_data),
        "important_temporal_slope_svg": str(slope_svg),
        "important_temporal_slope_png": str(slope_png),
        "important_temporal_interactive_html": str(temporal_html),
    })
    return outputs


def write_important_cooccurrence_context_outputs(sample_dir: Path, out_dir: Path, important_dir: Path, top_n: int = 20, min_support: int = 30) -> dict[str, str]:
    tables = important_dir / "tables"
    figures = important_dir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    features = [row for row in read_table(out_dir / "features" / "all_features.tsv") if row.get("presence", "1") != "0"]
    presence = feature_presence(features)
    sample_ids = sorted({row.get("assembly_accession", "") or row.get("sample_id", "") for row in features if row.get("assembly_accession", "") or row.get("sample_id", "")})
    sample_count = len(sample_ids)
    cooccurrence_source = read_table(out_dir / "cross_database" / "feature_cooccurrence.tsv")
    proximity_source = read_table(out_dir / "cross_database" / "feature_proximity_all.tsv") or read_table(out_dir / "cross_database" / "feature_proximity.tsv")

    pair_rows = []
    feature_category: dict[tuple[str, str], str] = {}
    feature_prevalence: dict[tuple[str, str], tuple[int, float]] = {}
    for row in features:
        key = (row.get("database", ""), row.get("feature_id", ""))
        if key[0] and key[1] and key not in feature_category:
            feature_category[key] = row.get("feature_category", "")
    for key, samples in presence.items():
        prevalence_value = len(samples) / sample_count if sample_count else 0.0
        feature_prevalence[key] = (len(samples), prevalence_value)

    for row in cooccurrence_source:
        n_total = int(_float_or_none(row.get("n_total", "")) or sample_count)
        n_both = int(_float_or_none(row.get("n_both_present", "")) or 0)
        n_a_only = int(_float_or_none(row.get("n_a_only", "")) or 0)
        n_b_only = int(_float_or_none(row.get("n_b_only", "")) or 0)
        n_neither = int(_float_or_none(row.get("n_neither", "")) or 0)
        phi = _float_or_none(row.get("phi", row.get("phi_correlation", ""))) or 0.0
        p_value = _float_or_none(row.get("p_value", ""))
        q_value = _float_or_none(row.get("q_value", ""))
        feature_a = (row.get("feature_a_database", ""), row.get("feature_a_id", ""))
        feature_b = (row.get("feature_b_database", ""), row.get("feature_b_id", ""))
        a_count = len(presence.get(feature_a, set()))
        b_count = len(presence.get(feature_b, set()))
        prevalence_a = a_count / sample_count if sample_count else 0.0
        prevalence_b = b_count / sample_count if sample_count else 0.0
        direction, significance = _cooccurrence_direction(phi, q_value)
        warnings = ["multiple_testing", "same_genome_only", "exploratory_only"]
        support_label = "supported"
        if n_total < min_support:
            support_label = "low_sample_support"
            significance = "insufficient_support"
            warnings.append("low_sample_support")
        if prevalence_a < 0.01 or prevalence_b < 0.01:
            warnings.append("low_feature_prevalence")
        pair_rows.append({
            "feature_a_database": feature_a[0],
            "feature_a_id": feature_a[1],
            "feature_b_database": feature_b[0],
            "feature_b_id": feature_b[1],
            "n_total": str(n_total),
            "n_both_present": str(n_both),
            "n_a_only": str(n_a_only),
            "n_b_only": str(n_b_only),
            "n_neither": str(n_neither),
            "prevalence_a": f"{prevalence_a:.4f}",
            "prevalence_b": f"{prevalence_b:.4f}",
            "cooccurrence_prevalence": f"{(n_both / n_total if n_total else 0.0):.4f}",
            "phi_correlation": f"{phi:.4f}",
            "odds_ratio": row.get("odds_ratio", ""),
            "p_value": "" if p_value is None else f"{p_value:.6g}",
            "q_value": "" if q_value is None else f"{q_value:.6g}",
            "direction": direction,
            "significance_label": significance,
            "support_label": support_label,
            "evidence_level": "same_genome",
            "warning_flags": ";".join(dict.fromkeys(warnings)),
        })

    pair_fields = [
        "feature_a_database", "feature_a_id", "feature_b_database", "feature_b_id",
        "n_total", "n_both_present", "n_a_only", "n_b_only", "n_neither",
        "prevalence_a", "prevalence_b", "cooccurrence_prevalence",
        "phi_correlation", "odds_ratio", "p_value", "q_value", "direction",
        "significance_label", "support_label", "evidence_level", "warning_flags",
    ]
    pair_summary_path = tables / "cooccurrence_pair_summary.tsv"
    write_rows(pair_summary_path, pair_rows, pair_fields)

    databases = sorted({key[0] for key in presence})
    x_database = "amr" if "amr" in databases else (databases[0] if databases else "all")
    y_database = "plasmidfinder" if "plasmidfinder" in databases else (next((db for db in databases if db != x_database), x_database) if databases else "all")
    if not any(row["feature_a_database"] == x_database and row["feature_b_database"] == y_database for row in pair_rows):
        if any(row["feature_b_database"] == x_database and row["feature_a_database"] == y_database for row in pair_rows):
            x_database, y_database = y_database, x_database
        elif pair_rows:
            x_database = pair_rows[0]["feature_a_database"]
            y_database = pair_rows[0]["feature_b_database"]

    x_features = [
        feature for feature, (_count, _prev) in sorted(
            feature_prevalence.items(),
            key=lambda item: (-item[1][0], item[0][0], item[0][1]),
        )
        if feature[0] == x_database
    ][:top_n]
    y_features = [
        feature for feature, (_count, _prev) in sorted(
            feature_prevalence.items(),
            key=lambda item: (-item[1][0], item[0][0], item[0][1]),
        )
        if feature[0] == y_database
    ][:top_n]
    by_pair = {}
    for row in pair_rows:
        key = ((row["feature_a_database"], row["feature_a_id"]), (row["feature_b_database"], row["feature_b_id"]))
        by_pair[key] = row
        reverse = ((row["feature_b_database"], row["feature_b_id"]), (row["feature_a_database"], row["feature_a_id"]))
        by_pair[reverse] = {
            **row,
            "feature_a_database": row["feature_b_database"],
            "feature_a_id": row["feature_b_id"],
            "feature_b_database": row["feature_a_database"],
            "feature_b_id": row["feature_a_id"],
            "n_a_only": row["n_b_only"],
            "n_b_only": row["n_a_only"],
            "prevalence_a": row["prevalence_b"],
            "prevalence_b": row["prevalence_a"],
        }

    heatmap_rows = []
    for feature_b in y_features:
        for feature_a in x_features:
            row = by_pair.get((feature_a, feature_b))
            if row:
                heatmap_rows.append(row)
            else:
                heatmap_rows.append({
                    "feature_a_database": feature_a[0],
                    "feature_a_id": feature_a[1],
                    "feature_b_database": feature_b[0],
                    "feature_b_id": feature_b[1],
                    "n_total": str(sample_count),
                    "n_both_present": "0",
                    "n_a_only": "0",
                    "n_b_only": "0",
                    "n_neither": str(sample_count),
                    "prevalence_a": f"{feature_prevalence.get(feature_a, (0, 0.0))[1]:.4f}",
                    "prevalence_b": f"{feature_prevalence.get(feature_b, (0, 0.0))[1]:.4f}",
                    "cooccurrence_prevalence": "0.0000",
                    "phi_correlation": "0.0000",
                    "odds_ratio": "",
                    "p_value": "",
                    "q_value": "",
                    "direction": "not_significant",
                    "significance_label": "insufficient_support" if sample_count < min_support else "not_significant",
                    "support_label": "low_sample_support" if sample_count < min_support else "supported",
                    "evidence_level": "same_genome",
                    "warning_flags": "low_sample_support;same_genome_only;exploratory_only" if sample_count < min_support else "same_genome_only;exploratory_only",
                })
    heatmap_fields = pair_fields + ["heatmap_color_rule"]
    for row in heatmap_rows:
        row["heatmap_color_rule"] = "colored" if row["significance_label"] in {"significant_positive", "significant_negative"} else "uncolored"
    heatmap_matrix_path = tables / "cooccurrence_heatmap_matrix.tsv"
    write_rows(heatmap_matrix_path, heatmap_rows, heatmap_fields)

    heatmap_base = f"cooccurrence_heatmap_{_safe_filename(x_database)}_vs_{_safe_filename(y_database)}"
    heatmap_data = figures / f"{heatmap_base}.data.tsv"
    heatmap_svg = figures / f"{heatmap_base}.svg"
    heatmap_png = figures / f"{heatmap_base}.png"
    heatmap_pdf = figures / f"{heatmap_base}.pdf"
    write_rows(heatmap_data, heatmap_rows, heatmap_fields)
    _write_cooccurrence_heatmap_svg(heatmap_svg, heatmap_rows, f"{x_database} vs {y_database} Co-occurrence")
    _write_cooccurrence_heatmap_png(heatmap_png, heatmap_rows)
    _write_simple_pdf(
        heatmap_pdf,
        f"{x_database} vs {y_database} Co-occurrence",
        [
            "PDF companion for the SVG/PNG co-occurrence heatmap.",
            "Red cells indicate significant positive association; blue cells indicate significant negative association.",
            "Gray cells are not significant or below support/effect thresholds.",
            "Use the SVG for publication-quality vector rendering and the data TSV for exact values.",
        ],
    )

    network_edges = [
        {
            "source_feature": f"{row['feature_a_database']}:{row['feature_a_id']}",
            "source_database": row["feature_a_database"],
            "target_feature": f"{row['feature_b_database']}:{row['feature_b_id']}",
            "target_database": row["feature_b_database"],
            "edge_weight": f"{abs(_float_or_none(row['phi_correlation']) or 0.0):.4f}",
            "phi_correlation": row["phi_correlation"],
            "q_value": row["q_value"],
            "direction": row["direction"],
            "evidence_level": row["evidence_level"],
            "n_both_present": row["n_both_present"],
        }
        for row in sorted(pair_rows, key=lambda item: -abs(_float_or_none(item.get("phi_correlation", "")) or 0.0))
        if row["significance_label"] in {"significant_positive", "significant_negative"}
    ][:100]
    node_ids = set()
    for edge in network_edges:
        node_ids.add(edge["source_feature"])
        node_ids.add(edge["target_feature"])
    node_rows = []
    for node_id in sorted(node_ids):
        database, feature_id = node_id.split(":", 1)
        positive_genomes, prevalence_value = feature_prevalence.get((database, feature_id), (0, 0.0))
        node_rows.append({
            "node_id": node_id,
            "feature_id": feature_id,
            "database": database,
            "feature_category": feature_category.get((database, feature_id), ""),
            "prevalence": f"{prevalence_value:.4f}",
            "positive_genomes": str(positive_genomes),
            "node_size": f"{max(4, 4 + prevalence_value * 24):.1f}",
            "node_label": feature_id,
        })
    edge_path = tables / "cooccurrence_network_edges.tsv"
    node_path = tables / "cooccurrence_network_nodes.tsv"
    edge_fields = ["source_feature", "source_database", "target_feature", "target_database", "edge_weight", "phi_correlation", "q_value", "direction", "evidence_level", "n_both_present"]
    node_fields = ["node_id", "feature_id", "database", "feature_category", "prevalence", "positive_genomes", "node_size", "node_label"]
    write_rows(edge_path, network_edges, edge_fields)
    write_rows(node_path, node_rows, node_fields)
    network_base = f"cooccurrence_network_{_safe_filename(x_database)}_vs_{_safe_filename(y_database)}"
    network_data = figures / f"{network_base}.data.tsv"
    network_svg = figures / f"{network_base}.svg"
    network_png = figures / f"{network_base}.png"
    network_pdf = figures / f"{network_base}.pdf"
    write_rows(network_data, network_edges, edge_fields)
    _write_cooccurrence_network_svg(network_svg, node_rows[:50], network_edges[:100], f"{x_database} vs {y_database} Co-occurrence Network")
    _write_cooccurrence_network_png(network_png, node_rows[:50], network_edges[:100])
    _write_simple_pdf(
        network_pdf,
        f"{x_database} vs {y_database} Co-occurrence Network",
        [
            "PDF companion for the co-occurrence network.",
            "Node color indicates database, node size indicates prevalence, and edge width indicates absolute phi correlation.",
            "Network edges are exploratory unless supported by same-contig/proximity evidence.",
        ],
    )

    feature_lookup = {}
    feature_lookup_fallback = {}
    for row in features:
        sample = row.get("assembly_accession", "") or row.get("sample_id", "")
        key = (sample, row.get("contig", ""), row.get("database", ""), row.get("feature_id", ""), row.get("start", ""), row.get("end", ""))
        feature_lookup[key] = row
        feature_lookup_fallback[(sample, row.get("database", ""), row.get("feature_id", ""))] = row

    context_rows = []
    for row in proximity_source:
        sample = row.get("assembly_accession", "")
        a_key = (sample, row.get("contig", ""), row.get("feature_a_database", ""), row.get("feature_a_id", ""), row.get("feature_a_start", ""), row.get("feature_a_end", ""))
        b_key = (sample, row.get("contig", ""), row.get("feature_b_database", ""), row.get("feature_b_id", ""), row.get("feature_b_start", ""), row.get("feature_b_end", ""))
        a_feature = feature_lookup.get(a_key) or feature_lookup_fallback.get((sample, row.get("feature_a_database", ""), row.get("feature_a_id", "")), {})
        b_feature = feature_lookup.get(b_key) or feature_lookup_fallback.get((sample, row.get("feature_b_database", ""), row.get("feature_b_id", "")), {})
        evidence = _context_evidence_label(row.get("evidence_level", "") or row.get("interpretation_level", ""))
        warnings = ["exploratory_only"]
        if evidence == "same_contig":
            warnings.append("same_genome_only" if not row.get("distance_bp") else "same_contig_context")
        if not row.get("feature_a_start") or not row.get("feature_b_start"):
            warnings.append("missing_coordinates")
        context_rows.append({
            "selected_database": row.get("feature_a_database", ""),
            "selected_feature": row.get("feature_a_id", ""),
            "context_database": row.get("feature_b_database", ""),
            "context_feature": row.get("feature_b_id", ""),
            "sample_id": sample,
            "assembly_accession": sample,
            "contig": row.get("contig", ""),
            "selected_start": row.get("feature_a_start", ""),
            "selected_end": row.get("feature_a_end", ""),
            "context_start": row.get("feature_b_start", ""),
            "context_end": row.get("feature_b_end", ""),
            "distance_bp": row.get("distance_bp", ""),
            "evidence_level": evidence,
            "identity_selected": a_feature.get("identity", ""),
            "coverage_selected": a_feature.get("coverage", ""),
            "identity_context": b_feature.get("identity", ""),
            "coverage_context": b_feature.get("coverage", ""),
            "interpretation_warning": row.get("interpretation_warning", "Same-contig/proximity evidence does not prove transfer, expression, phenotype, or plasmid localization."),
            "warning_flags": ";".join(dict.fromkeys(warnings)),
        })
    context_fields = [
        "selected_database", "selected_feature", "context_database", "context_feature",
        "sample_id", "assembly_accession", "contig", "selected_start", "selected_end",
        "context_start", "context_end", "distance_bp", "evidence_level",
        "identity_selected", "coverage_selected", "identity_context", "coverage_context",
        "interpretation_warning", "warning_flags",
    ]
    context_path = tables / "genomic_context_evidence.tsv"
    write_rows(context_path, context_rows, context_fields)

    selected_feature_key = None
    for row in context_rows:
        if row["selected_database"] in {"amr", "amrfinderplus"}:
            selected_feature_key = (row["selected_database"], row["selected_feature"])
            break
    if selected_feature_key is None and pair_rows:
        selected_feature_key = (pair_rows[0]["feature_a_database"], pair_rows[0]["feature_a_id"])
    selected_database, selected_feature = selected_feature_key if selected_feature_key else ("", "")

    same_genome_count = sum(
        1 for row in pair_rows
        if row["n_both_present"] != "0"
        and ((row["feature_a_database"], row["feature_a_id"]) == selected_feature_key or (row["feature_b_database"], row["feature_b_id"]) == selected_feature_key)
    )
    ladder_counts = Counter(row["evidence_level"] for row in context_rows if (row["selected_database"], row["selected_feature"]) == selected_feature_key)
    ladder_rows = [
        {"evidence_level": "same_genome", "count": str(same_genome_count)},
        {"evidence_level": "same_contig", "count": str(ladder_counts.get("same_contig", 0) + ladder_counts.get("within_10kb", 0) + ladder_counts.get("overlap_or_adjacent", 0))},
        {"evidence_level": "within_50kb", "count": str(ladder_counts.get("within_50kb", 0) + ladder_counts.get("within_10kb", 0) + ladder_counts.get("overlap_or_adjacent", 0))},
        {"evidence_level": "within_10kb", "count": str(ladder_counts.get("within_10kb", 0) + ladder_counts.get("overlap_or_adjacent", 0))},
        {"evidence_level": "overlap_or_adjacent", "count": str(ladder_counts.get("overlap_or_adjacent", 0))},
    ]
    ladder_base = f"genomic_context_evidence_ladder_{_safe_filename(selected_database)}_{_safe_filename(selected_feature)}"
    ladder_data = figures / f"{ladder_base}.data.tsv"
    ladder_svg = figures / f"{ladder_base}.svg"
    ladder_png = figures / f"{ladder_base}.png"
    ladder_pdf = figures / f"{ladder_base}.pdf"
    write_rows(ladder_data, ladder_rows, ["evidence_level", "count"])
    _write_context_ladder_svg(ladder_svg, ladder_rows, f"{selected_database}:{selected_feature} Context Evidence")
    _write_context_ladder_png(ladder_png, ladder_rows)
    _write_simple_pdf(
        ladder_pdf,
        f"{selected_database}:{selected_feature} Context Evidence",
        [f"{row['evidence_level']}: {row['count']}" for row in ladder_rows],
    )

    context_counter = Counter()
    for row in context_rows:
        if (row["selected_database"], row["selected_feature"]) == selected_feature_key:
            context_counter[(row["context_database"], row["context_feature"], row["evidence_level"])] += 1
    top_context_rows = [
        {
            "context_database": database,
            "context_feature": feature_id,
            "evidence_level": evidence,
            "count": str(count),
            "feature_label": f"{database}:{feature_id}",
        }
        for (database, feature_id, evidence), count in context_counter.most_common(top_n)
    ]
    top_context_base = f"top_context_features_{_safe_filename(selected_database)}_{_safe_filename(selected_feature)}"
    top_context_data = figures / f"{top_context_base}.data.tsv"
    top_context_svg = figures / f"{top_context_base}.svg"
    top_context_png = figures / f"{top_context_base}.png"
    top_context_pdf = figures / f"{top_context_base}.pdf"
    write_rows(top_context_data, top_context_rows, ["context_database", "context_feature", "evidence_level", "count", "feature_label"])
    _write_bar_svg(top_context_svg, top_context_rows, f"Top Context Features For {selected_feature}", "feature_label", "count", "Context evidence count")
    _write_bar_png(top_context_png, top_context_rows, "count")
    _write_simple_pdf(
        top_context_pdf,
        f"Top Context Features For {selected_feature}",
        [f"{row.get('feature_label', '')}: {row.get('count', '')} ({row.get('evidence_level', '')})" for row in top_context_rows],
    )

    neighborhood_rows = []
    neighborhood_sample = ""
    neighborhood_contig = ""
    selected_context = next((row for row in context_rows if row.get("contig") and row.get("selected_start") and row.get("context_start")), None)
    if selected_context:
        neighborhood_sample = selected_context["assembly_accession"]
        neighborhood_contig = selected_context["contig"]
        for row in features:
            sample = row.get("assembly_accession", "") or row.get("sample_id", "")
            if sample == neighborhood_sample and row.get("contig", "") == neighborhood_contig:
                neighborhood_rows.append({
                    "assembly_accession": sample,
                    "contig": row.get("contig", ""),
                    "database": row.get("database", ""),
                    "feature_id": row.get("feature_id", ""),
                    "start": row.get("start", ""),
                    "end": row.get("end", ""),
                    "strand": row.get("strand", ""),
                    "identity": row.get("identity", ""),
                    "coverage": row.get("coverage", ""),
                })
    neighborhood_path = tables / "contig_neighborhoods.tsv"
    neighborhood_fields = ["assembly_accession", "contig", "database", "feature_id", "start", "end", "strand", "identity", "coverage"]
    write_rows(neighborhood_path, neighborhood_rows, neighborhood_fields)
    neighborhood_base = f"contig_neighborhood_{_safe_filename(neighborhood_sample)}_{_safe_filename(neighborhood_contig)}" if neighborhood_rows else "contig_neighborhood_unavailable"
    neighborhood_data = figures / f"{neighborhood_base}.data.tsv"
    neighborhood_svg = figures / f"{neighborhood_base}.svg"
    neighborhood_png = figures / f"{neighborhood_base}.png"
    neighborhood_pdf = figures / f"{neighborhood_base}.pdf"
    write_rows(neighborhood_data, neighborhood_rows, neighborhood_fields)
    _write_contig_neighborhood_svg(neighborhood_svg, neighborhood_rows, f"{neighborhood_sample} {neighborhood_contig} Neighborhood")
    _write_contig_neighborhood_png(neighborhood_png, neighborhood_rows)
    _write_simple_pdf(
        neighborhood_pdf,
        f"{neighborhood_sample} {neighborhood_contig} Neighborhood",
        [f"{row.get('database', '')}:{row.get('feature_id', '')} {row.get('start', '')}-{row.get('end', '')}" for row in neighborhood_rows]
        or ["No coordinate-compatible contig neighborhood was available."],
    )

    same_contig_count = sum(1 for row in context_rows if row["evidence_level"] in {"same_contig", "within_10kb", "overlap_or_adjacent"})
    within_10kb_count = sum(1 for row in context_rows if row["evidence_level"] in {"within_10kb", "overlap_or_adjacent"})
    overlap_count = sum(1 for row in context_rows if row["evidence_level"] == "overlap_or_adjacent")
    significant_positive = sum(1 for row in pair_rows if row["significance_label"] == "significant_positive")
    significant_negative = sum(1 for row in pair_rows if row["significance_label"] == "significant_negative")
    summary_rows = [{
        "tested_pairs": str(len(pair_rows)),
        "significant_positive_pairs": str(significant_positive),
        "significant_negative_pairs": str(significant_negative),
        "same_contig_context_pairs": str(same_contig_count),
        "within_10kb_context_pairs": str(within_10kb_count),
        "overlap_or_adjacent_context_pairs": str(overlap_count),
        "selected_default_x_database": x_database,
        "selected_default_y_database": y_database,
        "message": "Co-occurrence is sample-level exploratory evidence. Same-contig/proximity evidence is stronger but does not prove transfer, expression, phenotype, or plasmid localization.",
    }]
    summary_path = tables / "cooccurrence_context_summary.tsv"
    write_rows(summary_path, summary_rows, ["tested_pairs", "significant_positive_pairs", "significant_negative_pairs", "same_contig_context_pairs", "within_10kb_context_pairs", "overlap_or_adjacent_context_pairs", "selected_default_x_database", "selected_default_y_database", "message"])
    cooccurrence_tables_zip = important_dir / "cooccurrence_tables.zip"
    cooccurrence_figures_zip = important_dir / "cooccurrence_figures.zip"
    _write_zip_bundle(
        cooccurrence_tables_zip,
        [pair_summary_path, heatmap_matrix_path, edge_path, node_path, context_path, neighborhood_path, summary_path],
        important_dir,
    )
    _write_zip_bundle(
        cooccurrence_figures_zip,
        [
            heatmap_svg, heatmap_png, heatmap_pdf, heatmap_data,
            network_svg, network_png, network_pdf, network_data,
            ladder_svg, ladder_png, ladder_pdf, ladder_data,
            top_context_svg, top_context_png, top_context_pdf, top_context_data,
            neighborhood_svg, neighborhood_png, neighborhood_pdf, neighborhood_data,
        ],
        important_dir,
    )

    interactive_html = figures / "cooccurrence_context.html"
    interactive_html.write_text(
        f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Co-occurrence / Genomic Context</title>
<style>
body {{ font-family: Arial, sans-serif; color: #1f2933; margin: 1.5rem; }}
label {{ font-weight: 700; margin-right: 0.35rem; }}
select {{ margin: 0 1rem 0.75rem 0; padding: 0.35rem; }}
.warning {{ background: #fff7ed; border-left: 4px solid #c2410c; padding: 0.75rem; margin: 0.75rem 0; }}
.heatmap-box {{ max-width: 100%; max-height: 700px; overflow: auto; border: 1px solid #d9e2ec; background: white; }}
table {{ border-collapse: collapse; width: 100%; font-size: 0.88rem; }}
th, td {{ border: 1px solid #d9e2ec; padding: 0.35rem; text-align: left; }}
th {{ background: #f0f4f8; }}
.panel {{ border: 1px solid #d9e2ec; background: #f8fafc; padding: 0.75rem; margin: 0.75rem 0; }}
.figure-row {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 1rem; }}
.figure-row img {{ max-width: 100%; border: 1px solid #d9e2ec; background: white; }}
</style></head><body>
<h1>Co-occurrence / Genomic Context</h1>
<div class="warning">Sample-level co-occurrence does not prove physical linkage. Same-contig/proximity evidence is stronger, but still does not prove transfer, expression, phenotype, or plasmid localization.</div>
<label for="analysisMode">Analysis mode</label><select id="analysisMode"><option value="heatmap" selected>Co-occurrence heatmap</option><option value="network">Co-occurrence network</option><option value="same_contig">Same-contig context</option><option value="proximity">Proximity context</option><option value="selected_feature">Selected feature report</option></select>
<label for="xDatabase">X database</label><select id="xDatabase"></select>
<label for="yDatabase">Y database</label><select id="yDatabase"></select>
<label for="featureSet">Feature set</label><select id="featureSet"><option value="10">Top 10</option><option value="20" selected>Top 20</option><option value="50">Top 50</option><option value="99999">Complete</option></select>
<label for="support">Minimum sample support</label><select id="support"><option value="0">All</option><option value="10">n >= 10</option><option value="20">n >= 20</option><option value="30" selected>n >= 30</option></select>
<label for="prevalence">Minimum feature prevalence</label><select id="prevalence"><option value="0">Any</option><option value="0.01" selected>&gt;= 1%</option><option value="0.05">&gt;= 5%</option><option value="0.10">&gt;= 10%</option></select>
<label for="significance">Significance</label><select id="significance"><option value="q0.05" selected>FDR q &lt; 0.05</option><option value="q0.10">FDR q &lt; 0.10</option><option value="p0.05">p &lt; 0.05</option><option value="all">show all</option></select>
<label for="effect">Effect size</label><select id="effect"><option value="0">Any</option><option value="0.2" selected>|phi| >= 0.2</option><option value="0.4">|phi| >= 0.4</option><option value="0.6">|phi| >= 0.6</option></select>
<label for="evidence">Evidence level</label><select id="evidence"><option value="same_genome" selected>same genome</option><option value="same_contig">same contig</option><option value="within_10kb">within 10 kb</option><option value="overlap_or_adjacent">overlap / adjacent</option></select>
<div id="summary"></div>
<div id="heatmapPanel" class="panel"><h2>Co-occurrence heatmap</h2><div class="heatmap-box" id="heatmap"></div></div>
<div id="networkPanel" class="panel"><h2>Co-occurrence network</h2><p>Network edges represent statistical co-occurrence unless supported by same-contig/proximity evidence.</p><img src="{html.escape(network_base)}.svg" alt="Co-occurrence network"><p><a href="{html.escape(network_base)}.png">PNG</a> | <a href="{html.escape(network_base)}.svg">SVG</a> | <a href="{html.escape(network_base)}.pdf">PDF</a> | <a href="{html.escape(network_base)}.data.tsv">Data TSV</a></p></div>
<div id="contextPanel" class="panel"><h2>Same-contig / proximity context</h2><p>Physical context evidence is shown separately from sample-level co-occurrence.</p><div class="figure-row"><div><h3>Evidence ladder</h3><img src="{html.escape(ladder_base)}.svg" alt="Context evidence ladder"><p><a href="{html.escape(ladder_base)}.png">PNG</a> | <a href="{html.escape(ladder_base)}.svg">SVG</a> | <a href="{html.escape(ladder_base)}.pdf">PDF</a> | <a href="{html.escape(ladder_base)}.data.tsv">Data TSV</a></p></div><div><h3>Top context features</h3><img src="{html.escape(top_context_base)}.svg" alt="Top context features"><p><a href="{html.escape(top_context_base)}.png">PNG</a> | <a href="{html.escape(top_context_base)}.svg">SVG</a> | <a href="{html.escape(top_context_base)}.pdf">PDF</a> | <a href="{html.escape(top_context_base)}.data.tsv">Data TSV</a></p></div></div></div>
<div id="selectedFeaturePanel" class="panel"><h2>Selected feature report</h2><p>Default selected feature: {html.escape(selected_database)}:{html.escape(selected_feature)}.</p><img src="{html.escape(neighborhood_base)}.svg" alt="Contig neighborhood"><p><a href="{html.escape(neighborhood_base)}.png">PNG</a> | <a href="{html.escape(neighborhood_base)}.svg">SVG</a> | <a href="{html.escape(neighborhood_base)}.pdf">PDF</a> | <a href="{html.escape(neighborhood_base)}.data.tsv">Data TSV</a></p></div>
<h2>Top Pair Table</h2><div id="pairTable"></div>
<p><a href="{html.escape(heatmap_base)}.png">Download default heatmap PNG</a> | <a href="{html.escape(heatmap_base)}.svg">Download default heatmap SVG</a> | <a href="{html.escape(heatmap_base)}.pdf">Download default heatmap PDF</a> | <a href="{html.escape(heatmap_base)}.data.tsv">Download default heatmap data</a> | <a href="../tables/cooccurrence_pair_summary.tsv">Download full pair table</a> | <a href="../cooccurrence_tables.zip">Download all co-occurrence tables ZIP</a> | <a href="../cooccurrence_figures.zip">Download all co-occurrence figures ZIP</a></p>
<script>
const pairs = {json.dumps(pair_rows)};
function num(value) {{ const n = Number(value); return Number.isFinite(n) ? n : 0; }}
const analysisMode = document.getElementById('analysisMode'), xSelect = document.getElementById('xDatabase'), ySelect = document.getElementById('yDatabase'), featureSet = document.getElementById('featureSet'), support = document.getElementById('support'), prevalence = document.getElementById('prevalence'), significance = document.getElementById('significance'), effect = document.getElementById('effect'), evidence = document.getElementById('evidence');
const databases = Array.from(new Set(pairs.flatMap(r => [r.feature_a_database, r.feature_b_database]))).filter(Boolean).sort();
function fillSelect(select, preferred) {{ select.innerHTML = databases.map(db => `<option value="${{db}}">${{db}}</option>`).join(''); if (databases.includes(preferred)) select.value = preferred; }}
fillSelect(xSelect, '{html.escape(x_database)}'); fillSelect(ySelect, '{html.escape(y_database)}');
function label(row, side) {{ return side === 'a' ? row.feature_a_database + ':' + row.feature_a_id : row.feature_b_database + ':' + row.feature_b_id; }}
function passesSignificance(r) {{
  if (significance.value === 'all') return true;
  if (significance.value === 'p0.05') return num(r.p_value) > 0 && num(r.p_value) < 0.05;
  if (significance.value === 'q0.10') return num(r.q_value) > 0 && num(r.q_value) < 0.10;
  return num(r.q_value) > 0 && num(r.q_value) < 0.05;
}}
function updatePanels() {{
  const mode = analysisMode.value;
  document.getElementById('heatmapPanel').style.display = mode === 'heatmap' ? 'block' : 'none';
  document.getElementById('networkPanel').style.display = mode === 'network' ? 'block' : 'none';
  document.getElementById('contextPanel').style.display = (mode === 'same_contig' || mode === 'proximity') ? 'block' : 'none';
  document.getElementById('selectedFeaturePanel').style.display = mode === 'selected_feature' ? 'block' : 'none';
}}
function render() {{
  updatePanels();
  const minSupport = Number(support.value), minEffect = Number(effect.value), minPrevalence = Number(prevalence.value), limit = Number(featureSet.value);
  let active = pairs.filter(r => r.feature_a_database === xSelect.value && r.feature_b_database === ySelect.value);
  if (!active.length) active = pairs.filter(r => r.feature_b_database === xSelect.value && r.feature_a_database === ySelect.value).map(r => Object.assign({{}}, r, {{feature_a_database: r.feature_b_database, feature_a_id: r.feature_b_id, feature_b_database: r.feature_a_database, feature_b_id: r.feature_a_id}}));
  active = active.filter(r => num(r.prevalence_a) >= minPrevalence && num(r.prevalence_b) >= minPrevalence);
  const xFeatures = Array.from(new Set(active.map(r => r.feature_a_id))).slice(0, limit);
  const yFeatures = Array.from(new Set(active.map(r => r.feature_b_id))).slice(0, limit);
  const byPair = new Map(active.map(r => [r.feature_a_id + '||' + r.feature_b_id, r]));
  let html = '<table><thead><tr><th></th>' + xFeatures.map(f => `<th>${{f}}</th>`).join('') + '</tr></thead><tbody>';
  for (const y of yFeatures) {{
    html += `<tr><th>${{y}}</th>`;
    for (const x of xFeatures) {{
      const r = byPair.get(x + '||' + y);
      let color = '#f1f5f9', text = '';
      if (r && num(r.n_total) >= minSupport && Math.abs(num(r.phi_correlation)) >= minEffect && passesSignificance(r) && (r.significance_label === 'significant_positive' || r.significance_label === 'significant_negative')) {{
        const intensity = Math.min(Math.abs(num(r.phi_correlation)), 1);
        color = r.significance_label === 'significant_positive' ? `rgb(254,${{Math.round(226 - 140 * intensity)}},${{Math.round(226 - 140 * intensity)}})` : `rgb(${{Math.round(219 - 150 * intensity)}},${{Math.round(234 - 120 * intensity)}},254)`;
        text = num(r.phi_correlation).toFixed(2);
      }}
      const tip = r ? `phi=${{r.phi_correlation}}; p=${{r.p_value}}; q=${{r.q_value}}; both=${{r.n_both_present}}/${{r.n_total}}; odds=${{r.odds_ratio}}; evidence=${{evidence.value}}; ${{r.significance_label}}` : 'No pair';
      html += `<td style="background:${{color}}" title="${{tip}}">${{text}}</td>`;
    }}
    html += '</tr>';
  }}
  html += '</tbody></table>';
  document.getElementById('heatmap').innerHTML = html;
  const significant = active.filter(r => r.significance_label === 'significant_positive' || r.significance_label === 'significant_negative');
  document.getElementById('summary').innerHTML = `<p>${{active.length}} tested pairs for ${{xSelect.value}} vs ${{ySelect.value}}; ${{significant.length}} significant pairs before current display filters. Gray cells are not significant or below support/effect/prevalence thresholds. Evidence level selector controls interpretation context; sample-level co-occurrence remains separate from same-contig/proximity evidence.</p>`;
  const top = active.slice().sort((a,b) => Math.abs(num(b.phi_correlation)) - Math.abs(num(a.phi_correlation))).slice(0, 50);
  const cols = ['feature_a_database','feature_a_id','feature_b_database','feature_b_id','n_total','n_both_present','phi_correlation','q_value','significance_label','warning_flags'];
  document.getElementById('pairTable').innerHTML = '<table><thead><tr>' + cols.map(c => `<th>${{c}}</th>`).join('') + '</tr></thead><tbody>' + top.map(r => '<tr>' + cols.map(c => `<td>${{r[c] || ''}}</td>`).join('') + '</tr>').join('') + '</tbody></table>';
}}
[analysisMode, xSelect, ySelect, featureSet, support, prevalence, significance, effect, evidence].forEach(el => el.addEventListener('change', render));
render();
</script></body></html>
""",
        encoding="utf-8",
    )

    outputs = {
        "important_cooccurrence_pair_summary": str(pair_summary_path),
        "important_cooccurrence_heatmap_matrix": str(heatmap_matrix_path),
        "important_cooccurrence_network_edges": str(edge_path),
        "important_cooccurrence_network_nodes": str(node_path),
        "important_genomic_context_evidence": str(context_path),
        "important_contig_neighborhoods": str(neighborhood_path),
        "important_cooccurrence_context_summary": str(summary_path),
        "important_cooccurrence_heatmap_svg": str(heatmap_svg),
        "important_cooccurrence_heatmap_png": str(heatmap_png),
        "important_cooccurrence_heatmap_pdf": str(heatmap_pdf),
        "important_cooccurrence_heatmap_data": str(heatmap_data),
        "important_cooccurrence_network_svg": str(network_svg),
        "important_cooccurrence_network_png": str(network_png),
        "important_cooccurrence_network_pdf": str(network_pdf),
        "important_cooccurrence_network_data": str(network_data),
        "important_context_ladder_svg": str(ladder_svg),
        "important_context_ladder_png": str(ladder_png),
        "important_context_ladder_pdf": str(ladder_pdf),
        "important_context_ladder_data": str(ladder_data),
        "important_top_context_features_svg": str(top_context_svg),
        "important_top_context_features_png": str(top_context_png),
        "important_top_context_features_pdf": str(top_context_pdf),
        "important_top_context_features_data": str(top_context_data),
        "important_contig_neighborhood_svg": str(neighborhood_svg),
        "important_contig_neighborhood_png": str(neighborhood_png),
        "important_contig_neighborhood_pdf": str(neighborhood_pdf),
        "important_contig_neighborhood_data": str(neighborhood_data),
        "important_cooccurrence_context_html": str(interactive_html),
        "important_cooccurrence_tables_zip": str(cooccurrence_tables_zip),
        "important_cooccurrence_figures_zip": str(cooccurrence_figures_zip),
    }
    return outputs


def write_important_results_report(
    sample_dir: Path,
    out_dir: Path,
    important_dir: Path,
    geographic_outputs: dict[str, str],
    qc_outputs: dict[str, str],
    prevalence_outputs: dict[str, str],
    variation_outputs: dict[str, str],
    temporal_outputs: dict[str, str],
    cooccurrence_outputs: dict[str, str],
) -> dict[str, str]:
    important_dir.mkdir(parents=True, exist_ok=True)
    basic_csv = sample_dir / "basic" / "enriched_genome_dataset.csv"
    dataset_rows = read_table(basic_csv)
    all_features = read_table(out_dir / "features" / "all_features.tsv")
    schema_summary = (out_dir / "manifest" / "schema_validation_summary.txt").read_text(encoding="utf-8") if (out_dir / "manifest" / "schema_validation_summary.txt").exists() else ""
    total_features = len(all_features)
    databases = sorted({row.get("database", "") for row in all_features if row.get("database")})
    qc_pass = sum(1 for row in dataset_rows if row.get("qc_pass") == "true")
    warning_count = sum(1 for row in dataset_rows if row.get("modules_warning"))
    cards = [
        ("Genomes", str(len(dataset_rows))),
        ("QC PASS", f"{qc_pass}/{len(dataset_rows)}" if dataset_rows else "0"),
        ("Feature rows", str(total_features)),
        ("Databases", str(len(databases))),
        ("Schema", "PASS" if "unmatched_feature_rows=0" in schema_summary and "invalid_feature_rows=0" in schema_summary else "CHECK"),
        ("Warnings", str(warning_count)),
    ]
    card_html = "".join(f"<div class='card'><span>{html.escape(label)}</span><strong>{html.escape(value)}</strong></div>" for label, value in cards)
    db_badges = " ".join(f"<span class='badge'>{html.escape(db)}</span>" for db in databases)
    qc_steps = read_table(important_dir / "key_tables" / "qc_step_summary.tsv")
    prevalence_rows = read_table(important_dir / "key_tables" / "feature_prevalence_summary.tsv")
    variation_rows = read_table(important_dir / "key_tables" / "feature_variation_summary.tsv")
    temporal_rows = read_table(important_dir / "key_tables" / "temporal_trend_summary.tsv")
    cooccurrence_rows = read_table(important_dir / "tables" / "cooccurrence_pair_summary.tsv")
    cooccurrence_summary_rows = read_table(important_dir / "tables" / "cooccurrence_context_summary.tsv")
    context_rows = read_table(important_dir / "tables" / "genomic_context_evidence.tsv")
    top_prevalence = sorted(prevalence_rows, key=lambda row: (row.get("database", ""), -(_float_or_none(row.get("prevalence_percent", "")) or 0.0), row.get("feature_id", "")))[:20]
    top_variation = sorted(variation_rows, key=lambda row: (-(_float_or_none(row.get("iqr_identity", "")) or 0.0), row.get("database", ""), row.get("feature_id", "")))[:20]
    top_temporal = sorted(temporal_rows, key=lambda row: (-abs(_float_or_none(row.get("change_percent_points", "")) or 0.0), row.get("database", ""), row.get("feature_id", "")))[:20]
    top_cooccurrence = sorted(cooccurrence_rows, key=lambda row: (-abs(_float_or_none(row.get("phi_correlation", "")) or 0.0), row.get("feature_a_database", ""), row.get("feature_a_id", "")))[:20]
    top_context = context_rows[:20]
    qc_table_html = _html_table(qc_steps, ["step_order", "qc_step", "tool", "enabled", "pass", "warning", "fail", "skipped", "status", "notes"], max_rows=20)
    prevalence_table_html = _html_table(top_prevalence, ["database", "feature_id", "positive_genomes", "sample_count", "prevalence_percent", "feature_rows"], max_rows=20)
    variation_table_html = _html_table(top_variation, ["database", "feature_id", "total_hits", "positive_genomes", "median_identity", "iqr_identity", "median_coverage", "iqr_coverage", "variation_label", "warning_flags"], max_rows=20)
    temporal_table_html = _html_table(top_temporal, ["database", "feature_id", "first_year", "last_year", "first_year_prevalence_percent", "last_year_prevalence_percent", "change_percent_points", "correlation", "trend_label", "support_label", "warning_flags"], max_rows=20)
    cooccurrence_table_html = _html_table(top_cooccurrence, ["feature_a_database", "feature_a_id", "feature_b_database", "feature_b_id", "n_total", "n_both_present", "phi_correlation", "q_value", "direction", "significance_label", "warning_flags"], max_rows=20)
    context_table_html = _html_table(top_context, ["selected_database", "selected_feature", "context_database", "context_feature", "assembly_accession", "contig", "distance_bp", "evidence_level", "warning_flags"], max_rows=20)
    prevalence_figures = []
    for path in sorted((important_dir / "figures").glob("prevalence_*_top20.svg")):
        database = path.name.replace("prevalence_", "").replace("_top20.svg", "")
        prevalence_figures.append(
            f"<div><h3>{html.escape(database)} prevalence</h3><img src='figures/{html.escape(path.name)}' alt='{html.escape(database)} prevalence'>"
            f"<p><a href='figures/{html.escape(path.with_suffix('.png').name)}'>PNG</a> | <a href='figures/{html.escape(path.name)}'>SVG</a> | <a href='figures/{html.escape(path.name.replace('.svg', '.data.tsv'))}'>Data TSV</a></p></div>"
        )
    prevalence_figures_html = "<div class='figure-row'>" + "".join(prevalence_figures[:6]) + "</div>" if prevalence_figures else "<p>No prevalence figures were generated because no feature rows were available.</p>"
    variation_figures = []
    for path in sorted((important_dir / "figures").glob("variation_identity_*_top20.svg")):
        database = path.name.replace("variation_identity_", "").replace("_top20.svg", "")
        variation_figures.append(
            f"<div><h3>{html.escape(database)} identity variation</h3><img src='figures/{html.escape(path.name)}' alt='{html.escape(database)} identity variation'>"
            f"<p><a href='figures/{html.escape(path.with_suffix('.png').name)}'>PNG</a> | <a href='figures/{html.escape(path.name)}'>SVG</a> | <a href='figures/{html.escape(path.name.replace('.svg', '.data.tsv'))}'>Data TSV</a></p></div>"
        )
    variation_figures_html = "<div class='figure-row'>" + "".join(variation_figures[:6]) + "</div>" if variation_figures else "<p>No variation figures were generated because no identity/coverage feature rows were available.</p>"
    temporal_figure_items = []
    for figure_name, title in [
        ("temporal_selected_feature_prevalence", "Selected feature prevalence"),
        ("temporal_slope_top40", "First-to-last-year slope plot"),
        ("temporal_database_burden_top20", "Database burden over time"),
        ("temporal_top_increasing_features", "Top increasing features"),
        ("temporal_top_decreasing_features", "Top decreasing features"),
        ("temporal_feature_heatmap_top40", "Temporal feature heatmap"),
    ]:
        svg_path = important_dir / "figures" / f"{figure_name}.svg"
        if not svg_path.exists():
            continue
        temporal_figure_items.append(
            f"<div><h3>{html.escape(title)}</h3><img src='figures/{html.escape(svg_path.name)}' alt='{html.escape(title)}'>"
            f"<p><a href='figures/{html.escape(figure_name)}.png'>PNG</a> | <a href='figures/{html.escape(figure_name)}.svg'>SVG</a> | <a href='figures/{html.escape(figure_name)}.data.tsv'>Data TSV</a></p></div>"
        )
    temporal_figures_html = "<div class='figure-row'>" + "".join(temporal_figure_items) + "</div>" if temporal_figure_items else "<p>No temporal figures were generated because collection-year metadata or feature rows were unavailable.</p>"
    cooccurrence_summary = cooccurrence_summary_rows[0] if cooccurrence_summary_rows else {}
    cooccurrence_summary_html = (
        "<p>"
        f"Co-occurrence analysis tested {html.escape(cooccurrence_summary.get('tested_pairs', '0'))} feature pairs, "
        f"including {html.escape(cooccurrence_summary.get('significant_positive_pairs', '0'))} significant positive pairs and "
        f"{html.escape(cooccurrence_summary.get('significant_negative_pairs', '0'))} significant negative pairs. "
        f"Same-contig/proximity evidence rows: {html.escape(cooccurrence_summary.get('same_contig_context_pairs', '0'))}; "
        f"within 10 kb: {html.escape(cooccurrence_summary.get('within_10kb_context_pairs', '0'))}."
        "</p>"
    )
    cooccurrence_figure_items = []
    for key, title in [
        ("important_cooccurrence_heatmap_svg", "Co-occurrence heatmap"),
        ("important_cooccurrence_network_svg", "Co-occurrence network"),
        ("important_context_ladder_svg", "Context evidence ladder"),
        ("important_top_context_features_svg", "Top context features"),
        ("important_contig_neighborhood_svg", "Contig neighborhood"),
    ]:
        figure_path = Path(cooccurrence_outputs.get(key, ""))
        if not figure_path.exists():
            continue
        stem = figure_path.stem
        png_name = figure_path.with_suffix(".png").name
        pdf_name = figure_path.with_suffix(".pdf").name
        data_name = f"{stem}.data.tsv"
        cooccurrence_figure_items.append(
            f"<div><h3>{html.escape(title)}</h3><img src='figures/{html.escape(figure_path.name)}' alt='{html.escape(title)}'>"
            f"<p><a href='figures/{html.escape(png_name)}'>PNG</a> | <a href='figures/{html.escape(figure_path.name)}'>SVG</a> | <a href='figures/{html.escape(pdf_name)}'>PDF</a> | <a href='figures/{html.escape(data_name)}'>Data TSV</a></p></div>"
        )
    cooccurrence_figures_html = "<div class='figure-row'>" + "".join(cooccurrence_figure_items) + "</div>" if cooccurrence_figure_items else "<p>No co-occurrence/context figures were generated because feature-pair data were unavailable.</p>"
    report_path = important_dir / "results.html"
    report_path.write_text(
        f"""<!doctype html>
<html><head><meta charset="utf-8"><title>PanResistome Important Results</title>
<style>
body {{ font-family: Arial, sans-serif; margin: 0; color: #1f2933; background: #f8fafc; }}
nav {{ position: fixed; left: 0; top: 0; bottom: 0; width: 210px; background: #102a43; color: white; padding: 1rem; }}
nav a {{ display: block; color: white; text-decoration: none; margin: 0.8rem 0; }}
main {{ margin-left: 240px; padding: 1.5rem 2rem; }}
section {{ background: white; border: 1px solid #d9e2ec; border-radius: 6px; padding: 1rem; margin-bottom: 1rem; }}
.cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(140px, 1fr)); gap: 0.75rem; }}
.card {{ border: 1px solid #d9e2ec; border-radius: 6px; padding: 0.75rem; background: #f8fafc; }}
.card span {{ display: block; font-size: 0.8rem; color: #52606d; }}
.card strong {{ display: block; font-size: 1.4rem; margin-top: 0.3rem; }}
.badge {{ display: inline-block; background: #e0f2fe; border: 1px solid #7dd3fc; padding: 0.2rem 0.45rem; border-radius: 999px; margin: 0.15rem; }}
.downloads a {{ display: inline-block; margin: 0.25rem 0.5rem 0.25rem 0; padding: 0.45rem 0.7rem; background: #0f766e; color: white; text-decoration: none; border-radius: 4px; }}
.warning {{ background: #fff7ed; border-left: 4px solid #c2410c; padding: 0.75rem; }}
iframe {{ width: 100%; height: 680px; border: 1px solid #d9e2ec; }}
table {{ border-collapse: collapse; width: 100%; margin: 0.75rem 0; font-size: 0.9rem; }}
th, td {{ border: 1px solid #d9e2ec; padding: 0.4rem; text-align: left; vertical-align: top; }}
th {{ background: #f0f4f8; }}
.figure-row {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(320px, 1fr)); gap: 1rem; }}
.figure-row img {{ max-width: 100%; border: 1px solid #d9e2ec; background: white; }}
</style></head>
<body>
<nav>
<h2>Results</h2>
<a href="#featured">Featured Results</a>
<a href="#overview">Run Overview</a>
<a href="#qc">QC Summary</a>
<a href="#prevalence">Prevalence</a>
<a href="#geography">Geographic Distribution</a>
<a href="#variations">Variations</a>
<a href="#temporal">Temporal Trends</a>
<a href="#cooccurrence">Co-occurrence / Genomic Context</a>
<a href="#files">Important Files</a>
<a href="#warnings">Warnings</a>
</nav>
<main>
<section id="featured"><h1>Featured Results</h1><div class="cards">{card_html}</div><p>{db_badges}</p></section>
<section id="overview"><h2>Run Overview</h2><p>This curated report summarizes the key outputs while preserving complete advanced outputs in the full PanResistome bundle.</p>
<div class="downloads"><a href="../basic/enriched_genome_dataset.csv">Download enriched dataset CSV</a><a href="../basic/enriched_genome_dataset.tsv">Download enriched dataset TSV</a><a href="../panr2_inputs/report/panr2_handoff_index.html">Open complete PanR2 handoff report</a></div></section>
<section id="qc"><h2>QC Summary</h2><p>This section shows which QC steps were enabled, skipped, or passed before annotation.</p>
<div class="figure-row"><div><h3>QC Funnel</h3><img src="figures/qc_funnel.svg" alt="QC funnel"></div><div><h3>QC Status</h3><img src="figures/qc_status_overview.svg" alt="QC status overview"></div></div>
{qc_table_html}
<div class="downloads"><a href="key_tables/qc_step_summary.tsv">Download QC step summary</a><a href="key_tables/qc_by_genome.tsv">Download per-genome QC table</a><a href="figures/qc_funnel.png">Download funnel PNG</a><a href="figures/qc_funnel.svg">Download funnel SVG</a><a href="figures/qc_funnel.data.tsv">Download funnel data</a></div></section>
<section id="prevalence"><h2>Prevalence</h2><p>Feature prevalence is summarized by database. Top plots are capped for readability; the complete prevalence table is downloadable.</p>
{prevalence_figures_html}
{prevalence_table_html}
<div class="downloads"><a href="key_tables/feature_prevalence_summary.tsv">Download complete prevalence table</a></div></section>
<section id="geography"><h2>Geographic Distribution</h2><div class="warning">Geographic patterns reflect the analyzed dataset only. They are not global prevalence estimates and can be affected by BioProject, lineage, country, and year sampling bias.</div>
<iframe src="figures/geographic_distribution_map.html" title="Geographic distribution map"></iframe>
<div class="downloads"><a href="figures/geographic_distribution_map.html">Open map</a><a href="figures/geographic_distribution_map.png">Download initial PNG</a><a href="figures/geographic_distribution_map.svg">Download initial SVG</a><a href="figures/geographic_distribution.data.tsv">Download data TSV</a></div></section>
<section id="variations"><h2>Variations</h2><p>Variation summaries use identity and coverage values when available. Low identity, low coverage, and high variation are review flags, not automatic failures.</p>
{variation_figures_html}
{variation_table_html}
<div class="downloads"><a href="key_tables/feature_variation_summary.tsv">Download variation summary</a><a href="key_tables/feature_variation_hits.tsv">Download hit-level variation table</a></div></section>
<section id="temporal"><h2>Temporal Trends</h2><div class="warning">Temporal trends reflect the analyzed dataset only. They can be affected by sampling year, BioProject, country, lineage, and missing collection-year metadata.</div>
<p>Prevalence trends use yearly percentages with genome-count denominators. Database burden is summarized as mean detected features per genome by collection year.</p>
<iframe src="figures/temporal_trends.html" title="Temporal trends interactive report"></iframe>
{temporal_figures_html}
{temporal_table_html}
<div class="downloads"><a href="figures/temporal_trends.html">Open interactive temporal report</a><a href="key_tables/temporal_database_burden.tsv">Download database burden by year</a><a href="key_tables/temporal_feature_prevalence.tsv">Download yearly feature prevalence</a><a href="key_tables/temporal_trend_summary.tsv">Download temporal trend summary</a><a href="key_tables/temporal_increasing_features.tsv">Download increasing features</a><a href="key_tables/temporal_decreasing_features.tsv">Download decreasing features</a></div></section>
<section id="cooccurrence"><h2>Co-occurrence / Genomic Context</h2><div class="warning">Sample-level co-occurrence does not prove physical linkage. Same-contig and proximity evidence are stronger context signals, but do not prove transfer, expression, phenotype, or plasmid localization.</div>
{cooccurrence_summary_html}
<iframe src="figures/cooccurrence_context.html" title="Co-occurrence and genomic context interactive report"></iframe>
{cooccurrence_figures_html}
<h3>Top co-occurrence pairs</h3>
{cooccurrence_table_html}
<h3>Genomic context evidence</h3>
{context_table_html}
<div class="downloads"><a href="figures/cooccurrence_context.html">Open interactive co-occurrence report</a><a href="cooccurrence_tables.zip">Download all co-occurrence tables ZIP</a><a href="cooccurrence_figures.zip">Download all co-occurrence figures ZIP</a><a href="tables/cooccurrence_pair_summary.tsv">Download pair summary</a><a href="tables/cooccurrence_heatmap_matrix.tsv">Download heatmap matrix</a><a href="tables/cooccurrence_network_edges.tsv">Download network edges</a><a href="tables/cooccurrence_network_nodes.tsv">Download network nodes</a><a href="tables/genomic_context_evidence.tsv">Download genomic context evidence</a><a href="tables/contig_neighborhoods.tsv">Download contig neighborhoods</a></div></section>
<section id="files"><h2>Important Files</h2><ul>
<li><a href="../basic/enriched_genome_dataset.csv">Enriched genome dataset CSV</a></li>
<li><a href="key_tables/qc_step_summary.tsv">QC step summary</a></li>
<li><a href="key_tables/feature_prevalence_summary.tsv">Feature prevalence summary</a></li>
<li><a href="key_tables/geographic_distribution.tsv">Geographic distribution table</a></li>
<li><a href="key_tables/feature_variation_summary.tsv">Feature variation summary</a></li>
<li><a href="key_tables/temporal_trend_summary.tsv">Temporal trend summary</a></li>
<li><a href="tables/cooccurrence_pair_summary.tsv">Co-occurrence pair summary</a></li>
<li><a href="tables/genomic_context_evidence.tsv">Genomic context evidence</a></li>
<li><a href="../panr2_inputs/features/all_features.tsv">Complete standardized feature table</a></li>
<li><a href="../panr2_inputs/manifest/schema_validation_summary.txt">Feature-contract validation summary</a></li>
</ul></section>
<section id="warnings"><h2>Warnings And Limitations</h2><p>Association, geography, and co-occurrence summaries are exploratory. Confirm important findings with denominator checks, lineage context, BioProject balance, and independent datasets.</p></section>
</main></body></html>
""",
        encoding="utf-8",
    )
    return {"important_results_html": str(report_path), **geographic_outputs, **qc_outputs, **prevalence_outputs, **variation_outputs, **temporal_outputs, **cooccurrence_outputs}


def write_user_output_bundles(
    sample_dir: Path,
    out_dir: Path,
    output_mode: str = "all",
    figure_formats: str = "png,svg,tsv",
    publication_figures: bool = False,
    pipeline_version: str = "",
) -> dict[str, str]:
    mode = (output_mode or "all").strip().lower()
    if mode not in {"basic", "important", "all"}:
        raise ValueError(f"Unsupported output_mode: {output_mode}")
    outputs = write_basic_enriched_dataset(sample_dir, out_dir, sample_dir / "basic", pipeline_version=pipeline_version)
    if mode in {"important", "all"}:
        important_dir = sample_dir / "important"
        manifest_dir = important_dir / "manifest"
        manifest_dir.mkdir(parents=True, exist_ok=True)
        qc_outputs = write_important_qc_outputs(sample_dir, out_dir, important_dir)
        prevalence_outputs = write_important_prevalence_outputs(sample_dir, out_dir, important_dir)
        geographic_outputs = write_important_geographic_outputs(sample_dir, out_dir, important_dir)
        variation_outputs = write_important_variation_outputs(sample_dir, out_dir, important_dir)
        temporal_outputs = write_important_temporal_outputs(sample_dir, out_dir, important_dir)
        cooccurrence_outputs = write_important_cooccurrence_context_outputs(sample_dir, out_dir, important_dir)
        outputs.update(write_important_results_report(sample_dir, out_dir, important_dir, geographic_outputs, qc_outputs, prevalence_outputs, variation_outputs, temporal_outputs, cooccurrence_outputs))
        write_rows(
            manifest_dir / "important_output_manifest.tsv",
            [
                {"setting": "output_mode", "value": mode, "message": "User-facing output bundle mode."},
                {"setting": "figure_formats_requested", "value": figure_formats, "message": "Requested figure formats. Current important report writes portable HTML, PNG, SVG, and TSV outputs without new plotting dependencies."},
                {"setting": "publication_figures", "value": str(publication_figures).lower(), "message": "Reserved for PDF/static publication figure expansion."},
            ],
            ["setting", "value", "message"],
        )
    return outputs


def write_interpretation_reports(
    out_dir: Path,
    report_mode: str = "publication",
    max_metadata_columns: int = 80,
    skip_heavy_interactive_plots: bool = False,
) -> dict[str, str]:
    report_dir = out_dir / "report"
    report_dir.mkdir(parents=True, exist_ok=True)
    analysis_dir = out_dir / "metadata_feature_analysis"
    cross_dir = out_dir / "cross_database"
    diversity_dir = out_dir / "diversity"
    manifest_dir = out_dir / "manifest"
    row_limit = 40 if report_mode == "compact" else (200 if report_mode == "exploratory" else 80)
    metadata_row_limit = max_metadata_columns if report_mode == "compact" else max(max_metadata_columns, row_limit)

    styles = """
body { font-family: Arial, sans-serif; margin: 2rem; color: #1f2933; }
h1, h2 { color: #102a43; }
table { border-collapse: collapse; width: 100%; margin: 1rem 0; font-size: 0.9rem; }
th, td { border: 1px solid #d9e2ec; padding: 0.45rem; text-align: left; vertical-align: top; }
th { background: #f0f4f8; }
.warning { background: #fff7ed; border-left: 4px solid #c2410c; padding: 0.75rem; }
.grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 0.75rem; }
.card { border: 1px solid #d9e2ec; padding: 0.75rem; border-radius: 6px; }
"""

    def page(title: str, body: str) -> str:
        return (
            "<!doctype html><html><head><meta charset='utf-8'>"
            f"<title>{html.escape(title)}</title><style>{styles}</style></head><body>"
            f"<h1>{html.escape(title)}</h1>{body}</body></html>\n"
        )

    top_findings_rows = read_table(analysis_dir / "top_findings.tsv")
    top_feature_rows = read_table(analysis_dir / "top_features_by_database.tsv")
    report_control_rows = read_table(manifest_dir / "report_controls.tsv")
    top_findings_path = report_dir / "top_findings.html"
    top_findings_body = """
<div class="warning">These are screening summaries. They indicate metadata-feature patterns in this dataset and do not prove causality, transmission, plasmid localization, or physical linkage.</div>
<p>Inspect metadata completeness, group sizes, and BioProject/study balance before interpreting any association.</p>
"""
    top_findings_body += _html_table(
        top_findings_rows,
        [
            "finding_type",
            "database",
            "feature_id",
            "metadata_column",
            "metadata_value",
            "effect_size",
            "supporting_samples",
            "largest_bioproject_fraction",
            "dominant_ST",
            "dominant_ST_fraction",
            "dominant_ani_cluster",
            "dominant_ani_cluster_fraction",
            "lineage_warning_flags",
            "warning_flags",
            "interpretation_label",
            "message",
        ],
        max_rows=row_limit,
    )
    top_findings_body += "<h2>Top Features By Database</h2>"
    top_findings_body += _html_table(
        top_feature_rows,
        ["database", "rank", "feature_id", "feature_category", "present_count", "prevalence"],
        max_rows=row_limit,
    )
    top_findings_path.write_text(page("Top Metadata-Feature Findings", top_findings_body), encoding="utf-8")

    metadata_path = report_dir / "metadata_quality_and_bias.html"
    metadata_rows = read_table(analysis_dir / "metadata_column_eligibility.tsv")
    audit_rows = read_table(analysis_dir / "fetchm2_metadata_audit.tsv")
    usability_rows = read_table(analysis_dir / "metadata_usability_summary.tsv")
    bioproject_rows = read_table(analysis_dir / "bioproject_bias_report.tsv")
    metadata_body = """
<p>This page summarizes which FetchM2 metadata columns are usable for comparative analysis and which are sparse, dominated by one value, or likely identifiers.</p>
"""
    metadata_body += "<h2>Metadata Usability Summary</h2>"
    metadata_body += _html_table(
        usability_rows,
        ["metadata_column", "non_missing_count", "missing_fraction", "unique_values", "largest_group_fraction", "recommended_for_analysis", "reason"],
        max_rows=metadata_row_limit,
    )
    metadata_body += "<h2>Metadata Column Eligibility</h2>"
    metadata_body += _html_table(
        metadata_rows,
        ["metadata_column", "data_type", "non_missing_count", "missing_fraction", "unique_values", "largest_group_fraction", "eligible", "reason"],
        max_rows=metadata_row_limit,
    )
    metadata_body += "<h2>FetchM2 Metadata Audit</h2>"
    metadata_body += _html_table(
        audit_rows,
        ["column", "standardized_name", "data_type", "non_missing_count", "missing_fraction", "unique_values", "top_value", "top_value_fraction", "recommended_for_analysis", "reason"],
        max_rows=metadata_row_limit,
    )
    metadata_body += "<h2>BioProject / Study Bias</h2>"
    metadata_body += "<p>Dominance by one BioProject can make metadata-feature associations reflect study design rather than biology.</p>"
    metadata_body += _html_table(
        bioproject_rows,
        ["row_type", "database", "feature_id", "metadata_column", "metadata_value", "samples_evaluated", "n_bioprojects", "largest_bioproject", "largest_bioproject_fraction", "status", "warning"],
        max_rows=row_limit,
    )
    metadata_path.write_text(page("Metadata Quality And Bias", metadata_body), encoding="utf-8")

    bioproject_path = report_dir / "bioproject_bias.html"
    bioproject_body = """
<div class="warning">Public genome sets often reflect study design. Treat associations dominated by one BioProject as exploratory until validated with independent sampling.</div>
"""
    bioproject_body += _html_table(
        bioproject_rows,
        ["row_type", "database", "feature_id", "metadata_column", "metadata_value", "samples_evaluated", "n_bioprojects", "largest_bioproject", "largest_bioproject_count", "largest_bioproject_fraction", "status", "warning"],
        max_rows=row_limit,
    )
    bioproject_path.write_text(page("BioProject Bias", bioproject_body), encoding="utf-8")

    burden_path = report_dir / "database_burden_by_metadata.html"
    burden_rows = read_table(analysis_dir / "database_burden_by_sample.tsv")
    burden_body = """
<p>Database burden is the number of unique features observed per sample for each database. Use this with metadata columns to identify groups worth deeper analysis.</p>
"""
    burden_body += _html_table(
        burden_rows,
        ["assembly_accession", "database", "unique_feature_count", "category_count", "country", "host", "sample_type", "isolation_source", "environment_medium", "collection_year"],
        max_rows=row_limit,
    )
    burden_path.write_text(page("Database Burden By Metadata", burden_body), encoding="utf-8")

    lineage_path = report_dir / "lineage_context.html"
    lineage_rows = read_table(analysis_dir / "lineage_summary.tsv")
    lineage_burden_rows = read_table(analysis_dir / "lineage_feature_burden.tsv")
    lineage_confounding_rows = read_table(analysis_dir / "lineage_metadata_confounding.tsv")
    lineage_warning_rows = read_table(analysis_dir / "lineage_adjusted_warnings.tsv")
    lineage_body = """
<div class="warning">Lineage context helps identify when metadata-feature associations may reflect ST, ANI cluster, BioProject, or sampling structure rather than independent biology.</div>
<p>These warnings do not invalidate findings; they mark results that need lineage-aware interpretation or independent validation.</p>
"""
    lineage_body += "<h2>Sample Lineage Summary</h2>"
    lineage_body += _html_table(lineage_rows, ["assembly_accession", "mlst_ST", "ani_cluster", "bioproject", "lineage_data_status"], max_rows=row_limit)
    lineage_body += "<h2>Lineage Feature Burden</h2>"
    lineage_body += _html_table(lineage_burden_rows, ["lineage_type", "lineage_value", "database", "feature_id", "lineage_sample_count", "present_count", "prevalence"], max_rows=row_limit)
    lineage_body += "<h2>Metadata-Lineage Confounding</h2>"
    lineage_body += _html_table(lineage_confounding_rows, ["metadata_column", "lineage_type", "lineage_value", "samples_evaluated", "dominant_metadata_value", "dominant_metadata_fraction", "status", "warning"], max_rows=row_limit)
    lineage_body += "<h2>Top-Finding Lineage Warnings</h2>"
    lineage_body += _html_table(lineage_warning_rows, ["finding_type", "database", "feature_id", "metadata_column", "metadata_value", "dominant_ST", "dominant_ST_fraction", "dominant_ani_cluster", "dominant_ani_cluster_fraction", "lineage_warning_flags", "interpretation"], max_rows=row_limit)
    lineage_path.write_text(page("Lineage Context", lineage_body), encoding="utf-8")

    diversity_path = report_dir / "diversity_summary.html"
    richness_rows = read_table(diversity_dir / "feature_richness_by_sample.tsv")
    database_diversity_rows = read_table(diversity_dir / "database_diversity_by_sample.tsv")
    core_rows = read_table(diversity_dir / "core_accessory_rare_features.tsv")
    accumulation_rows = read_table(diversity_dir / "pan_feature_accumulation.tsv")
    diversity_body = """
<p>Diversity summaries are computed from standardized feature tables. They support broad comparison of feature burden, database richness, core/accessory/rare features, and sample-level Jaccard distances.</p>
"""
    diversity_body += "<h2>Feature Richness By Sample</h2>"
    diversity_body += _html_table(richness_rows, ["assembly_accession", "total_feature_richness", "database_count", "shannon_database_diversity", "simpson_database_diversity", "country", "host", "isolation_source"], max_rows=row_limit)
    diversity_body += "<h2>Database Richness By Sample</h2>"
    diversity_body += _html_table(database_diversity_rows, ["assembly_accession", "database", "feature_richness"], max_rows=row_limit)
    diversity_body += "<h2>Core / Accessory / Rare Features</h2>"
    diversity_body += _html_table(core_rows, ["database", "feature_id", "present_count", "sample_count", "prevalence", "feature_class"], max_rows=row_limit)
    diversity_body += "<h2>Pan-Feature Accumulation</h2>"
    diversity_body += _html_table(accumulation_rows, ["sample_order", "assembly_accession", "cumulative_unique_features", "new_features_added"], max_rows=row_limit)
    diversity_path.write_text(page("Feature Diversity Summary", diversity_body), encoding="utf-8")

    statistical_path = report_dir / "statistical_summary.html"
    statistical_rows = read_table(analysis_dir / "statistical_summary.tsv")
    statistical_body = """
<div class="warning">Statistical summaries are exploratory screening summaries. Use warning flags and validation datasets before drawing biological conclusions.</div>
"""
    statistical_body += _html_table(statistical_rows, ["metric", "value", "message"], max_rows=80)
    statistical_path.write_text(page("Statistical Summary", statistical_body), encoding="utf-8")

    cross_path = report_dir / "cross_database_interpretation.html"
    proximity_rows = read_table(cross_dir / "feature_proximity.tsv")
    cooccurrence_rows = read_table(cross_dir / "database_cooccurrence_summary.tsv")
    concordance_rows = read_table(cross_dir / "amrfinder_abricate_concordance.tsv")
    cross_body = """
<div class="warning">Genome-level co-occurrence means features were detected in the same sample/genome. Same-contig and proximity rows provide stronger context, but still do not prove transfer, expression, phenotype, or plasmid localization.</div>
<p>Evidence levels are ordered from weaker sample-level context to stronger same-contig coordinate context. Assembly fragmentation and naming differences can still limit interpretation.</p>
<div class="grid">
"""
    if skip_heavy_interactive_plots:
        cross_body += "<div class='card'><strong>Large-dataset mode</strong><br>Heavy interactive plots are skipped or deprioritized when supported; use TSV summaries first.</div>"
    for label, path in [
        ("Feature co-occurrence", cross_dir / "feature_cooccurrence.tsv"),
        ("Database co-occurrence", cross_dir / "database_cooccurrence_summary.tsv"),
        ("AMR-MGE same-contig evidence", cross_dir / "amr_mge_same_contig.tsv"),
        ("AMR-plasmid same-contig evidence", cross_dir / "amr_plasmid_same_contig.tsv"),
        ("AMR-integron same-contig evidence", cross_dir / "amr_integron_same_contig.tsv"),
        ("Report-capped feature proximity evidence", cross_dir / "feature_proximity.tsv"),
        ("Complete feature proximity evidence", cross_dir / "feature_proximity_all.tsv"),
        ("AMRFinderPlus vs ABRicate concordance", cross_dir / "amrfinder_abricate_concordance.tsv"),
    ]:
        cross_body += f"<div class='card'><strong>{html.escape(label)}</strong><br><a href='../{html.escape(_relative_link(path, out_dir))}'>{html.escape(_relative_link(path, out_dir))}</a></div>"
    cross_body += "</div><h2>Database Co-occurrence</h2>"
    cross_body += _html_table(cooccurrence_rows, ["database_a", "database_b", "n_total", "n_both_present", "jaccard"], max_rows=row_limit)
    cross_body += "<h2>Same-Contig And Proximity Evidence</h2>"
    cross_body += _html_table(
        proximity_rows,
        ["assembly_accession", "contig", "context", "feature_a_database", "feature_a_id", "feature_b_database", "feature_b_id", "distance_bp", "evidence_level", "interpretation_warning"],
        max_rows=row_limit,
    )
    cross_body += "<h2>AMRFinderPlus vs ABRicate Concordance</h2>"
    cross_body += _html_table(
        concordance_rows,
        ["feature_id", "sample_id", "abricate_feature_ids", "amrfinderplus_feature_ids", "shared_sample_count", "status", "possible_match_basis", "interpretation_note"],
        max_rows=row_limit,
    )
    cross_path.write_text(page("Cross-Database Interpretation", cross_body), encoding="utf-8")

    concordance_path = report_dir / "amrfinder_abricate_concordance.html"
    concordance_body = """
<div class="warning">ABRicate and AMRFinderPlus use different databases, algorithms, and naming conventions. Discordance should be inspected in raw outputs before drawing biological conclusions.</div>
"""
    concordance_body += _html_table(
        concordance_rows,
        ["feature_id", "sample_id", "normalized_feature_id", "abricate_feature_ids", "amrfinderplus_feature_ids", "abricate_sample_count", "amrfinderplus_sample_count", "shared_sample_count", "status", "possible_match_basis"],
        max_rows=row_limit,
    )
    concordance_path.write_text(page("AMRFinderPlus vs ABRicate Concordance", concordance_body), encoding="utf-8")

    setup_path = report_dir / "database_setup_and_contract.html"
    setup_rows = read_table(manifest_dir / "database_setup_status.tsv")
    abricate_setup_rows = read_table(manifest_dir / "abricate_database_setup_status.tsv")
    mobsuite_setup_rows = read_table(manifest_dir / "mobsuite_database_setup_status.tsv")
    genomad_setup_rows = read_table(manifest_dir / "genomad_database_setup_status.tsv")
    audit_rows = read_table(manifest_dir / "feature_completeness_audit.tsv")
    setup_body = "<h2>Database And Tool Setup</h2>"
    setup_body += _html_table(
        setup_rows,
        ["database_or_tool", "required_for_profile", "checked", "status", "setup_action", "version_or_path", "message"],
        max_rows=80,
    )
    if abricate_setup_rows:
        setup_body += "<h2>ABRicate Database Setup Actions</h2>"
        setup_body += _html_table(
            abricate_setup_rows,
            ["database", "present_before", "setup_requested", "update_requested", "setup_status", "update_status", "present_after", "status", "message"],
            max_rows=80,
        )
    if mobsuite_setup_rows:
        setup_body += "<h2>MOB-suite Database Setup Actions</h2>"
        setup_body += _html_table(
            mobsuite_setup_rows,
            ["database_dir", "auto_init_requested", "auto_init_taxa_requested", "mob_init_status", "taxa_init_status", "core_status", "taxa_status", "status", "message"],
            max_rows=80,
        )
    if genomad_setup_rows:
        setup_body += "<h2>geNomad Database Setup Actions</h2>"
        setup_body += _html_table(
            genomad_setup_rows,
            ["requested_database_dir", "resolved_database_dir", "auto_download_requested", "genomad_available", "download_status", "status", "message"],
            max_rows=80,
        )
    setup_body += "<h2>Feature Completeness Audit</h2>"
    setup_body += _html_table(
        audit_rows,
        ["database", "expected_from_profile", "module_enabled", "feature_table_found", "feature_rows", "unique_features", "samples_with_features", "status", "message"],
        max_rows=row_limit,
    )
    setup_path.write_text(page("Database Setup And Feature Contract", setup_body), encoding="utf-8")

    controls_path = report_dir / "report_controls.html"
    controls_body = """
<p>This page records report density, large-dataset safeguards, and output-size limits applied while building the PanR2 handoff bundle.</p>
"""
    controls_body += _html_table(report_control_rows, ["setting", "value", "message"], max_rows=80)
    controls_path.write_text(page("Report Controls", controls_body), encoding="utf-8")

    index_path = report_dir / "panr2_handoff_index.html"
    index_body = "<p>PanResistome-generated PanR2 handoff interpretation pages.</p><ul>"
    for label, path in [
        ("Top findings", top_findings_path),
        ("Metadata quality and bias", metadata_path),
        ("BioProject bias", bioproject_path),
        ("Database burden by metadata", burden_path),
        ("Lineage context", lineage_path),
        ("Feature diversity summary", diversity_path),
        ("Statistical summary", statistical_path),
        ("Cross-database interpretation", cross_path),
        ("AMRFinderPlus vs ABRicate concordance", concordance_path),
        ("Database setup and feature contract", setup_path),
        ("Report controls", controls_path),
    ]:
        index_body += f"<li><a href='{html.escape(path.name)}'>{html.escape(label)}</a></li>"
    index_body += "</ul>"
    index_path.write_text(page("PanR2 Handoff Report Index", index_body), encoding="utf-8")

    return {
        "handoff_report_index": str(index_path),
        "top_findings_html": str(top_findings_path),
        "metadata_quality_html": str(metadata_path),
        "bioproject_bias_html": str(bioproject_path),
        "database_burden_html": str(burden_path),
        "lineage_context_html": str(lineage_path),
        "diversity_summary_html": str(diversity_path),
        "statistical_summary_html": str(statistical_path),
        "cross_database_interpretation_html": str(cross_path),
        "amrfinder_abricate_concordance_html": str(concordance_path),
        "database_setup_contract_html": str(setup_path),
        "report_controls_html": str(controls_path),
    }


def export_contract(
    sample_dir: Path,
    out_dir: Path,
    large_dataset: bool = False,
    report_mode: str = "publication",
    max_features_heatmap: int = 300,
    max_features_network: int = 300,
    max_metadata_columns: int = 80,
    top_n_features_per_database: int = 25,
    skip_heavy_interactive_plots: bool = False,
    core_feature_threshold: float = 0.95,
    rare_feature_threshold: float = 0.05,
    output_mode: str = "all",
    figure_formats: str = "png,svg,tsv",
    publication_figures: bool = False,
    pipeline_version: str = "",
) -> dict[str, str]:
    written = write_feature_tables(sample_dir, out_dir)
    validation = validate_feature_tables(sample_dir, out_dir)
    all_features = read_table(out_dir / "features" / "all_features.tsv")
    metadata_rows = load_metadata_rows(sample_dir)
    metadata_outputs = write_metadata_analysis(sample_dir, out_dir)
    feature_outputs = write_feature_eligibility_and_prevalence(
        all_features,
        metadata_rows,
        out_dir,
        top_n_features_per_database=top_n_features_per_database,
    )
    matrix_outputs = write_feature_matrices(
        all_features,
        metadata_rows,
        out_dir,
        max_features=max_features_heatmap,
    )
    diversity_outputs = write_diversity_analysis(
        all_features,
        metadata_rows,
        out_dir,
        core_feature_threshold=core_feature_threshold,
        rare_feature_threshold=rare_feature_threshold,
    )
    cross_outputs = write_cross_database_outputs(
        all_features,
        metadata_rows,
        out_dir,
        max_features=max_features_network,
    )
    audit_outputs = write_feature_completeness_audit(sample_dir, out_dir)
    contract_outputs = write_feature_contract_manifest(out_dir)
    control_outputs = write_report_controls(
        out_dir,
        all_features,
        metadata_rows,
        large_dataset=large_dataset,
        report_mode=report_mode,
        max_features_heatmap=max_features_heatmap,
        max_features_network=max_features_network,
        max_metadata_columns=max_metadata_columns,
        top_n_features_per_database=top_n_features_per_database,
        skip_heavy_interactive_plots=skip_heavy_interactive_plots,
    )
    report_outputs = write_interpretation_reports(
        out_dir,
        report_mode=report_mode,
        max_metadata_columns=max_metadata_columns,
        skip_heavy_interactive_plots=skip_heavy_interactive_plots,
    )
    user_outputs = write_user_output_bundles(
        sample_dir,
        out_dir,
        output_mode=output_mode,
        figure_formats=figure_formats,
        publication_figures=publication_figures,
        pipeline_version=pipeline_version,
    )
    return {
        **written,
        **validation,
        **metadata_outputs,
        **feature_outputs,
        **matrix_outputs,
        **diversity_outputs,
        **cross_outputs,
        **audit_outputs,
        **contract_outputs,
        **control_outputs,
        **report_outputs,
        **user_outputs,
    }
