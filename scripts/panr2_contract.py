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
        {"setting": "important_lineage_default_top_n", "value": "20", "message": "Default Top-N lineages/features shown in important lineage figures."},
        {"setting": "important_lineage_feature_cap_per_database", "value": "200", "message": "Report-facing feature cap per database for important lineage enrichment; complete lineage TSVs are preserved."},
        {"setting": "important_lineage_complete_tsvs_preserved", "value": "true", "message": "Complete important lineage distribution, burden, enrichment, and selected-feature TSVs are written even when figures are capped."},
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


COUNTRY_REGIONS = {
    "argentina": ("South America", "South America"),
    "australia": ("Oceania", "Australia and New Zealand"),
    "austria": ("Europe", "Western Europe"),
    "bangladesh": ("Asia", "South Asia"),
    "belgium": ("Europe", "Western Europe"),
    "brazil": ("South America", "South America"),
    "canada": ("North America", "Northern America"),
    "chile": ("South America", "South America"),
    "china": ("Asia", "East Asia"),
    "colombia": ("South America", "South America"),
    "denmark": ("Europe", "Northern Europe"),
    "egypt": ("Africa", "Northern Africa"),
    "finland": ("Europe", "Northern Europe"),
    "france": ("Europe", "Western Europe"),
    "germany": ("Europe", "Western Europe"),
    "ghana": ("Africa", "Western Africa"),
    "greece": ("Europe", "Southern Europe"),
    "india": ("Asia", "South Asia"),
    "indonesia": ("Asia", "Southeast Asia"),
    "iran": ("Asia", "Western Asia"),
    "iraq": ("Asia", "Western Asia"),
    "ireland": ("Europe", "Northern Europe"),
    "israel": ("Asia", "Western Asia"),
    "italy": ("Europe", "Southern Europe"),
    "japan": ("Asia", "East Asia"),
    "kenya": ("Africa", "Eastern Africa"),
    "malaysia": ("Asia", "Southeast Asia"),
    "mexico": ("North America", "Central America"),
    "nepal": ("Asia", "South Asia"),
    "netherlands": ("Europe", "Western Europe"),
    "new zealand": ("Oceania", "Australia and New Zealand"),
    "nigeria": ("Africa", "Western Africa"),
    "norway": ("Europe", "Northern Europe"),
    "pakistan": ("Asia", "South Asia"),
    "peru": ("South America", "South America"),
    "philippines": ("Asia", "Southeast Asia"),
    "poland": ("Europe", "Eastern Europe"),
    "portugal": ("Europe", "Southern Europe"),
    "russia": ("Europe", "Eastern Europe"),
    "saudi arabia": ("Asia", "Western Asia"),
    "singapore": ("Asia", "Southeast Asia"),
    "south africa": ("Africa", "Southern Africa"),
    "south korea": ("Asia", "East Asia"),
    "spain": ("Europe", "Southern Europe"),
    "sri lanka": ("Asia", "South Asia"),
    "sweden": ("Europe", "Northern Europe"),
    "switzerland": ("Europe", "Western Europe"),
    "taiwan": ("Asia", "East Asia"),
    "thailand": ("Asia", "Southeast Asia"),
    "turkey": ("Asia", "Western Asia"),
    "united kingdom": ("Europe", "Northern Europe"),
    "uk": ("Europe", "Northern Europe"),
    "united states": ("North America", "Northern America"),
    "usa": ("North America", "Northern America"),
    "vietnam": ("Asia", "Southeast Asia"),
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


def _country_region(country: str) -> tuple[str, str]:
    return COUNTRY_REGIONS.get(_clean_country(country).lower(), ("", ""))


def _geo_group_size_label(total: int) -> str:
    if total < 3:
        return "very_small_group"
    if total < 5:
        return "small_group"
    if total < 10:
        return "limited_group"
    if total < 30:
        return "exploratory_group"
    return "standard_group"


def _geo_percent_display(percent: float, positive: int, total: int) -> str:
    return f"{percent:.1f}% ({positive}/{total})"


def _geo_metric_value(row: dict[str, str], metric: str) -> float:
    field = {
        "prevalence_percent": "prevalence_percent",
        "positive_genomes": "positive_genomes",
        "total_genomes": "total_genomes",
        "mean_feature_burden_per_genome": "mean_feature_burden_per_genome",
        "median_feature_burden_per_genome": "median_feature_burden_per_genome",
        "feature_rows": "feature_rows",
        "total_feature_rows": "total_feature_rows",
    }.get(metric, "prevalence_percent")
    return _float_or_none(row.get(field, "")) or 0.0


def _svg_geographic_map(rows: list[dict[str, str]], title: str) -> str:
    width, height = 960, 480
    points = []
    for row in rows:
        xy = _country_xy(row.get("country", ""), width, height)
        if not xy:
            continue
        prevalence = (_float_or_none(row.get("prevalence_percent", "")) or 0.0) / 100.0
        if row.get("prevalence_percent", "") == "":
            prevalence = _float_or_none(row.get("prevalence", "")) or 0.0
        total = int(_float_or_none(row.get("total_genomes", "")) or 0)
        radius = max(5, min(28, 4 + math.sqrt(max(total, 1)) * 4))
        red = int(230 * prevalence)
        blue = int(200 * (1 - prevalence))
        fill = "#cbd5e1" if "small_group_warning" in row.get("warning_flags", "") else f"rgb({red},80,{blue})"
        positive = row.get("positive_genomes", "") or row.get("positive_genomes_with_database", "0")
        label = f"{row.get('country', '')}: {positive}/{row.get('total_genomes', '0')} ({prevalence * 100:.1f}%); warnings={row.get('warning_flags', '')}"
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
        prevalence = (_float_or_none(row.get("prevalence_percent", "")) or 0.0) / 100.0
        if row.get("prevalence_percent", "") == "":
            prevalence = _float_or_none(row.get("prevalence", "")) or 0.0
        total = int(_float_or_none(row.get("total_genomes", "")) or 0)
        radius = max(5, min(28, 4 + math.sqrt(max(total, 1)) * 4))
        color = (203, 213, 225) if "small_group_warning" in row.get("warning_flags", "") else (int(230 * prevalence), 80, int(200 * (1 - prevalence)))
        draw_circle(xy[0], xy[1] + header, radius, color)
    _write_png(path, width, height, pixels)


GEOGRAPHIC_DATABASE_FIELDS = [
    "database",
    "geo_level",
    "group_name",
    "country",
    "continent",
    "subcontinent",
    "collection_year",
    "total_genomes",
    "positive_genomes_with_database",
    "positive_genomes",
    "prevalence_percent",
    "prevalence_display",
    "total_feature_rows",
    "feature_rows",
    "mean_feature_burden_per_genome",
    "median_feature_burden_per_genome",
    "min_collection_year",
    "max_collection_year",
    "top_bioproject",
    "largest_bioproject_fraction",
    "dominant_lineage",
    "dominant_lineage_fraction",
    "group_size_label",
    "warning_flags",
    "interpretation_label",
]


GEOGRAPHIC_FEATURE_FIELDS = [
    "database",
    "feature_id",
    "feature_name",
    "geo_level",
    "group_name",
    "country",
    "continent",
    "subcontinent",
    "collection_year",
    "total_genomes",
    "positive_genomes",
    "prevalence_percent",
    "prevalence_display",
    "feature_rows",
    "mean_hits_per_positive_genome",
    "min_collection_year",
    "max_collection_year",
    "top_bioproject",
    "largest_bioproject_fraction",
    "dominant_lineage",
    "dominant_lineage_fraction",
    "group_size_label",
    "warning_flags",
    "interpretation_label",
]


def _write_geographic_bar_figure(
    figures: Path,
    stem: str,
    rows: list[dict[str, str]],
    title: str,
    fields: list[str],
    value_field: str = "prevalence_percent",
) -> list[Path]:
    plot_rows = []
    for row in rows:
        plot = dict(row)
        plot["feature_id"] = row.get("group_name") or row.get("country") or row.get("continent") or row.get("subcontinent")
        plot["prevalence_display"] = row.get("prevalence_display") or row.get(value_field, "")
        plot_rows.append(plot)
    data_path = figures / f"{stem}.data.tsv"
    write_rows(data_path, rows, fields)
    svg_path = figures / f"{stem}.svg"
    png_path = figures / f"{stem}.png"
    pdf_path = figures / f"{stem}.pdf"
    _write_prevalence_bar_svg(svg_path, plot_rows, title, value_field)
    _write_bar_png(png_path, plot_rows, value_field)
    _write_simple_pdf(
        pdf_path,
        title,
        [
            f"{row.get('group_name', '')}: {row.get('prevalence_display', row.get(value_field, ''))}; warnings={row.get('warning_flags', '')}"
            for row in rows[:30]
        ],
    )
    return [data_path, svg_path, png_path, pdf_path]


def _write_geographic_map_figure(
    figures: Path,
    stem: str,
    rows: list[dict[str, str]],
    title: str,
    fields: list[str],
) -> list[Path]:
    data_path = figures / f"{stem}.data.tsv"
    write_rows(data_path, rows, fields)
    svg_path = figures / f"{stem}.svg"
    png_path = figures / f"{stem}.png"
    pdf_path = figures / f"{stem}.pdf"
    svg_path.write_text(_svg_geographic_map(rows, title), encoding="utf-8")
    _geographic_map_png(rows, png_path)
    _write_simple_pdf(
        pdf_path,
        title,
        [
            f"{row.get('country', '')}: {row.get('prevalence_display', '')}; warnings={row.get('warning_flags', '')}"
            for row in rows[:30]
        ],
    )
    return [data_path, svg_path, png_path, pdf_path]


def write_important_geographic_outputs(sample_dir: Path, out_dir: Path, important_dir: Path, top_n: int = 20) -> dict[str, str]:
    key_tables = important_dir / "key_tables"
    tables = important_dir / "tables"
    figures = important_dir / "figures"
    key_tables.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    metadata_rows = normalize_metadata_rows(load_metadata_rows(sample_dir))
    features = read_table(out_dir / "features" / "all_features.tsv")
    feature_samples = {
        row.get("assembly_accession", "") or row.get("sample_id", "")
        for row in features
        if row.get("assembly_accession", "") or row.get("sample_id", "")
    }
    metadata_samples = {row.get("assembly_accession", "") for row in metadata_rows if row.get("assembly_accession")}
    samples = sorted(feature_samples | metadata_samples)
    metadata_by_sample = {row["assembly_accession"]: row for row in metadata_rows if row.get("assembly_accession")}
    presence = feature_presence(features)
    database_presence: dict[str, set[str]] = defaultdict(set)
    database_counts_by_sample: dict[tuple[str, str], int] = defaultdict(int)
    feature_counts_by_sample: dict[tuple[str, str, str], int] = defaultdict(int)
    feature_metadata: dict[tuple[str, str], dict[str, str]] = {}
    lineage_by_sample: dict[str, str] = {}
    for feature in features:
        sample = feature.get("assembly_accession", "") or feature.get("sample_id", "")
        database = feature.get("database", "")
        feature_id = feature.get("feature_id", "")
        if not sample or not database or not feature_id or feature.get("presence", "1") == "0":
            continue
        database_presence[database].add(sample)
        database_counts_by_sample[(sample, database)] += 1
        feature_counts_by_sample[(sample, database, feature_id)] += 1
        feature_metadata.setdefault((database, feature_id), feature)
        if database == "mlst" and "ST-" not in feature_id:
            category = feature.get("feature_category", "").lower()
            if "sequence_type" in category or feature_id.lower().startswith("st") or ":st" in feature_id.lower():
                lineage_by_sample.setdefault(sample, feature_id)
    for sample, meta in metadata_by_sample.items():
        lineage = first_value(meta, ["dominant_lineage_label", "mlst_ST", "lineage", "ani_cluster", "mash_cluster"], "")
        if lineage:
            lineage_by_sample.setdefault(sample, lineage)

    feature_rank = sorted(presence.items(), key=lambda item: (-len(item[1]), item[0][0], item[0][1]))
    selected_features = feature_rank[:max(top_n * 2, top_n, 20)]
    selected_feature_keys = {key for key, _ in selected_features}

    def enriched_meta(sample: str) -> dict[str, str]:
        meta = dict(metadata_by_sample.get(sample, {}))
        country = _clean_country(meta.get("country", ""))
        derived_continent, derived_subcontinent = _country_region(country)
        meta["country"] = country
        meta["continent"] = meta.get("continent", "") or derived_continent
        meta["subcontinent"] = meta.get("subcontinent", "") or derived_subcontinent
        return meta

    enriched_metadata = {sample: enriched_meta(sample) for sample in samples}
    missing_country_count = sum(1 for sample in samples if not enriched_metadata.get(sample, {}).get("country"))
    missing_country_fraction = missing_country_count / len(samples) if samples else 0.0

    grouped_samples: dict[str, dict[str, set[str]]] = {
        "country": defaultdict(set),
        "continent": defaultdict(set),
        "subcontinent": defaultdict(set),
        "country_year": defaultdict(set),
    }
    for sample in samples:
        meta = enriched_metadata.get(sample, {})
        country = meta.get("country", "")
        continent = meta.get("continent", "")
        subcontinent = meta.get("subcontinent", "")
        year = meta.get("collection_year", "")
        grouped_samples["country"][country or "missing"].add(sample)
        grouped_samples["continent"][continent or "unknown"].add(sample)
        grouped_samples["subcontinent"][subcontinent or "unknown"].add(sample)
        if country and year:
            grouped_samples["country_year"][f"{country} ({year})"].add(sample)
        elif country:
            grouped_samples["country_year"][f"{country} (unknown)"].add(sample)
        else:
            grouped_samples["country_year"]["missing (unknown)"].add(sample)

    def group_context(geo_level: str, group_name: str, members: set[str]) -> dict[str, str]:
        total = len(members)
        metas = [enriched_metadata.get(sample, {}) for sample in members]
        countries = [meta.get("country", "") for meta in metas if meta.get("country", "")]
        continents = [meta.get("continent", "") for meta in metas if meta.get("continent", "")]
        subcontinents = [meta.get("subcontinent", "") for meta in metas if meta.get("subcontinent", "")]
        years = [meta.get("collection_year", "") for meta in metas if meta.get("collection_year", "")]
        bioprojects = [meta.get("bioproject", "") for meta in metas if meta.get("bioproject", "")]
        lineages = [lineage_by_sample.get(sample, "") for sample in members if lineage_by_sample.get(sample, "")]

        def dominant(values: list[str]) -> tuple[str, float]:
            if not values or not total:
                return "", 0.0
            label, count = Counter(values).most_common(1)[0]
            return label, count / total

        top_bioproject, top_bioproject_fraction = dominant(bioprojects)
        dominant_lineage, dominant_lineage_fraction = dominant(lineages)
        top_country, top_country_fraction = dominant(countries)
        top_year, top_year_fraction = dominant(years)
        warnings = []
        if total < 5:
            warnings.append("small_group_warning")
        if geo_level in {"country", "country_year"} and (not countries or group_name.startswith("missing")):
            warnings.append("missing_country_metadata")
        if geo_level in {"continent", "subcontinent"} and group_name == "unknown":
            warnings.append("missing_region_metadata")
        if top_bioproject_fraction >= 0.8 and total >= 5:
            warnings.append("bioproject_dominance")
        if dominant_lineage_fraction >= 0.8 and total >= 5:
            warnings.append("lineage_dominance")
        if top_year_fraction >= 0.8 and total >= 5:
            warnings.append("collection_year_bias")
        if geo_level in {"continent", "subcontinent"} and top_country_fraction >= 0.8 and total >= 5:
            warnings.append("single_country_dominance")
        warnings.append("exploratory_only")
        min_year = min(years) if years else ""
        max_year = max(years) if years else ""
        return {
            "group_name": group_name,
            "country": top_country if geo_level not in {"country", "country_year"} else (countries[0] if countries else "missing"),
            "continent": continents[0] if len(set(continents)) == 1 else (group_name if geo_level == "continent" else (continents[0] if continents else "unknown")),
            "subcontinent": subcontinents[0] if len(set(subcontinents)) == 1 else (group_name if geo_level == "subcontinent" else (subcontinents[0] if subcontinents else "unknown")),
            "collection_year": top_year if geo_level == "country_year" else "",
            "min_collection_year": min_year,
            "max_collection_year": max_year,
            "top_bioproject": top_bioproject,
            "largest_bioproject_fraction": f"{top_bioproject_fraction:.3f}" if top_bioproject else "",
            "dominant_lineage": dominant_lineage,
            "dominant_lineage_fraction": f"{dominant_lineage_fraction:.3f}" if dominant_lineage else "",
            "group_size_label": _geo_group_size_label(total),
            "warning_flags": ";".join(dict.fromkeys(warnings)),
            "interpretation_label": "exploratory",
        }

    database_rows = []
    for database, present_samples in sorted(database_presence.items()):
        for geo_level, groups in grouped_samples.items():
            for group_name, members in sorted(groups.items()):
                total = len(members)
                positive = len(present_samples & members)
                prevalence_percent = positive / total * 100 if total else 0.0
                counts = [database_counts_by_sample.get((sample, database), 0) for sample in members]
                total_rows = sum(counts)
                database_rows.append({
                    "database": database,
                    "geo_level": geo_level,
                    **group_context(geo_level, group_name, members),
                    "total_genomes": str(total),
                    "positive_genomes_with_database": str(positive),
                    "positive_genomes": str(positive),
                    "prevalence_percent": f"{prevalence_percent:.1f}",
                    "prevalence_display": _geo_percent_display(prevalence_percent, positive, total),
                    "total_feature_rows": str(total_rows),
                    "feature_rows": str(total_rows),
                    "mean_feature_burden_per_genome": f"{_mean([float(value) for value in counts]):.2f}" if counts else "0.00",
                    "median_feature_burden_per_genome": f"{_median([float(value) for value in counts]):.2f}" if counts else "0.00",
                })

    feature_rows = []
    for (database, feature_id), present_samples in sorted(presence.items(), key=lambda item: (item[0][0], item[0][1])):
        feature_meta = feature_metadata.get((database, feature_id), {})
        feature_name = feature_meta.get("feature_name", "") or feature_id
        for geo_level, groups in grouped_samples.items():
            for group_name, members in sorted(groups.items()):
                total = len(members)
                positive_samples = present_samples & members
                positive = len(positive_samples)
                prevalence_percent = positive / total * 100 if total else 0.0
                feature_rows_count = sum(feature_counts_by_sample.get((sample, database, feature_id), 0) for sample in members)
                feature_rows.append({
                    "database": database,
                    "feature_id": feature_id,
                    "feature_name": feature_name,
                    "geo_level": geo_level,
                    **group_context(geo_level, group_name, members),
                    "total_genomes": str(total),
                    "positive_genomes": str(positive),
                    "prevalence_percent": f"{prevalence_percent:.1f}",
                    "prevalence_display": _geo_percent_display(prevalence_percent, positive, total),
                    "feature_rows": str(feature_rows_count),
                    "mean_hits_per_positive_genome": f"{(feature_rows_count / positive):.2f}" if positive else "0.00",
                })

    warning_rows = []
    for geo_level, groups in grouped_samples.items():
        for group_name, members in sorted(groups.items()):
            context = group_context(geo_level, group_name, members)
            warning_rows.append({
                "geo_level": geo_level,
                "group_name": group_name,
                "total_genomes": str(len(members)),
                "missing_metadata": str(sum(1 for sample in members if not enriched_metadata.get(sample, {}).get("country"))),
                "largest_bioproject": context.get("top_bioproject", ""),
                "largest_bioproject_fraction": context.get("largest_bioproject_fraction", ""),
                "dominant_lineage": context.get("dominant_lineage", ""),
                "dominant_lineage_fraction": context.get("dominant_lineage_fraction", ""),
                "warning_flags": context.get("warning_flags", ""),
            })

    def summarize_rows(rows: list[dict[str, str]], mode: str, database: str, feature_id: str, feature_name: str, geo_level: str) -> dict[str, str]:
        candidates = [row for row in rows if row.get("geo_level") == geo_level]
        passing = [
            row for row in candidates
            if int(_float_or_none(row.get("total_genomes", "")) or 0) >= 5
            and row.get("group_name") not in {"missing", "unknown", "missing (unknown)"}
        ]
        ranked = sorted(
            passing or candidates,
            key=lambda row: (
                -(_float_or_none(row.get("prevalence_percent", "")) or 0.0),
                -(_float_or_none(row.get("positive_genomes", row.get("positive_genomes_with_database", ""))) or 0.0),
                row.get("group_name", ""),
            ),
        )
        top = ranked[0] if ranked else {}
        warning_flags = sorted({flag for row in candidates for flag in row.get("warning_flags", "").split(";") if flag})
        return {
            "database": database,
            "mode": mode,
            "feature_id": feature_id,
            "feature_name": feature_name,
            "geo_level": geo_level,
            "metric": "prevalence_percent",
            "total_geographic_groups": str(len(candidates)),
            "groups_passing_min_n": str(len(passing)),
            "missing_country_count": str(missing_country_count),
            "missing_country_fraction": f"{missing_country_fraction:.3f}",
            "top_group": top.get("group_name", ""),
            "top_group_prevalence_percent": top.get("prevalence_percent", ""),
            "top_group_positive_genomes": top.get("positive_genomes", top.get("positive_genomes_with_database", "")),
            "top_group_total_genomes": top.get("total_genomes", ""),
            "warning_flags": ";".join(warning_flags),
        }

    summary_fields = [
        "database", "mode", "feature_id", "feature_name", "geo_level", "metric",
        "total_geographic_groups", "groups_passing_min_n", "missing_country_count",
        "missing_country_fraction", "top_group", "top_group_prevalence_percent",
        "top_group_positive_genomes", "top_group_total_genomes", "warning_flags",
    ]
    summary_rows = []
    for database in sorted(database_presence):
        for geo_level in ["country", "continent", "subcontinent", "country_year"]:
            rows_for_key = [row for row in database_rows if row.get("database") == database and row.get("geo_level") == geo_level]
            summary_rows.append(summarize_rows(rows_for_key, "database_burden", database, "__any_feature__", "", geo_level))
    for (database, feature_id), _present_samples in selected_features:
        feature_name = feature_metadata.get((database, feature_id), {}).get("feature_name", "") or feature_id
        for geo_level in ["country", "continent", "subcontinent", "country_year"]:
            rows_for_key = [row for row in feature_rows if row.get("database") == database and row.get("feature_id") == feature_id and row.get("geo_level") == geo_level]
            summary_rows.append(summarize_rows(rows_for_key, "individual_feature", database, feature_id, feature_name, geo_level))

    summary_path = tables / "geographic_distribution_summary.tsv"
    feature_path = tables / "geographic_feature_distribution.tsv"
    burden_path = tables / "geographic_database_burden.tsv"
    warning_path = tables / "geographic_warning_summary.tsv"
    write_rows(summary_path, summary_rows, summary_fields)
    write_rows(feature_path, feature_rows, GEOGRAPHIC_FEATURE_FIELDS)
    write_rows(burden_path, database_rows, GEOGRAPHIC_DATABASE_FIELDS)
    write_rows(warning_path, warning_rows, ["geo_level", "group_name", "total_genomes", "missing_metadata", "largest_bioproject", "largest_bioproject_fraction", "dominant_lineage", "dominant_lineage_fraction", "warning_flags"])

    legacy_fields = ["mode", "database", "feature_id", "country", "continent", "subcontinent", "collection_year", "total_genomes", "positive_genomes", "prevalence", "prevalence_percent", "warning_flags"]
    legacy_rows = []
    for row in database_rows:
        legacy_rows.append({
            "mode": "database_burden",
            "database": row.get("database", ""),
            "feature_id": "__any_feature__",
            "country": row.get("country", ""),
            "continent": row.get("continent", ""),
            "subcontinent": row.get("subcontinent", ""),
            "collection_year": row.get("collection_year", "") or "all",
            "total_genomes": row.get("total_genomes", ""),
            "positive_genomes": row.get("positive_genomes", ""),
            "prevalence": f"{((_float_or_none(row.get('prevalence_percent', '')) or 0.0) / 100.0):.4f}",
            "prevalence_percent": row.get("prevalence_percent", ""),
            "warning_flags": row.get("warning_flags", ""),
        })
    for row in [row for row in feature_rows if (row.get("database", ""), row.get("feature_id", "")) in selected_feature_keys]:
        legacy_rows.append({
            "mode": "feature",
            "database": row.get("database", ""),
            "feature_id": row.get("feature_id", ""),
            "country": row.get("country", ""),
            "continent": row.get("continent", ""),
            "subcontinent": row.get("subcontinent", ""),
            "collection_year": row.get("collection_year", "") or "all",
            "total_genomes": row.get("total_genomes", ""),
            "positive_genomes": row.get("positive_genomes", ""),
            "prevalence": f"{((_float_or_none(row.get('prevalence_percent', '')) or 0.0) / 100.0):.4f}",
            "prevalence_percent": row.get("prevalence_percent", ""),
            "warning_flags": row.get("warning_flags", ""),
        })
    data_path = key_tables / "geographic_distribution.tsv"
    write_rows(data_path, legacy_rows, legacy_fields)
    (figures / "geographic_distribution.data.tsv").write_text(data_path.read_text(encoding="utf-8"), encoding="utf-8")

    def plot_ready(rows: list[dict[str, str]], level: str, limit: int = top_n) -> list[dict[str, str]]:
        active = [row for row in rows if row.get("geo_level") == level and row.get("group_name") not in {"missing", "unknown", "missing (unknown)"}]
        passing = [row for row in active if int(_float_or_none(row.get("total_genomes", "")) or 0) >= 5]
        ranked = sorted(
            passing or active,
            key=lambda row: (
                -_geo_metric_value(row, "prevalence_percent"),
                -(_float_or_none(row.get("positive_genomes", row.get("positive_genomes_with_database", ""))) or 0.0),
                row.get("group_name", ""),
            ),
        )
        return ranked[:limit]

    figure_files: list[Path] = []
    default_database = "amr" if "amr" in database_presence else (sorted(database_presence)[0] if database_presence else "")
    if default_database:
        default_rows = [row for row in database_rows if row.get("database") == default_database]
        country_rows = [row for row in default_rows if row.get("geo_level") == "country" and row.get("country") not in {"", "missing"}]
        figure_files += _write_geographic_map_figure(figures, f"geographic_map_{_safe_filename(default_database)}_burden", country_rows, f"{default_database} geographic distribution", GEOGRAPHIC_DATABASE_FIELDS)
        figure_files += _write_geographic_bar_figure(figures, f"geographic_country_bar_{_safe_filename(default_database)}_burden", plot_ready(default_rows, "country"), f"{default_database} by country", GEOGRAPHIC_DATABASE_FIELDS)
        figure_files += _write_geographic_bar_figure(figures, f"geographic_continent_bar_{_safe_filename(default_database)}_burden", plot_ready(default_rows, "continent"), f"{default_database} by continent", GEOGRAPHIC_DATABASE_FIELDS)
        figure_files += _write_geographic_bar_figure(figures, f"geographic_region_bar_{_safe_filename(default_database)}_burden", plot_ready(default_rows, "subcontinent"), f"{default_database} by region", GEOGRAPHIC_DATABASE_FIELDS)

    if selected_features:
        feature_database, feature_id = selected_features[0][0]
        feature_default_rows = [row for row in feature_rows if row.get("database") == feature_database and row.get("feature_id") == feature_id]
        feature_stem = f"{_safe_filename(feature_database)}_{_safe_filename(feature_id)}"
        figure_files += _write_geographic_map_figure(figures, f"geographic_map_{feature_stem}", [row for row in feature_default_rows if row.get("geo_level") == "country" and row.get("country") not in {"", "missing"}], f"{feature_database}:{feature_id} geographic distribution", GEOGRAPHIC_FEATURE_FIELDS)
        figure_files += _write_geographic_bar_figure(figures, f"geographic_country_bar_{feature_stem}", plot_ready(feature_default_rows, "country"), f"{feature_database}:{feature_id} by country", GEOGRAPHIC_FEATURE_FIELDS)
        figure_files += _write_geographic_bar_figure(figures, f"geographic_continent_bar_{feature_stem}", plot_ready(feature_default_rows, "continent"), f"{feature_database}:{feature_id} by continent", GEOGRAPHIC_FEATURE_FIELDS)
        figure_files += _write_geographic_bar_figure(figures, f"geographic_region_bar_{feature_stem}", plot_ready(feature_default_rows, "subcontinent"), f"{feature_database}:{feature_id} by region", GEOGRAPHIC_FEATURE_FIELDS)

    initial_rows = [row for row in database_rows if row.get("database") == default_database and row.get("geo_level") == "country" and row.get("country") not in {"", "missing"}] if default_database else []
    svg_path = figures / "geographic_distribution_map.svg"
    svg_path.write_text(_svg_geographic_map(initial_rows, "Geographic Distribution"), encoding="utf-8")
    png_path = figures / "geographic_distribution_map.png"
    _geographic_map_png(initial_rows, png_path)

    report_feature_rows = [row for row in feature_rows if (row.get("database", ""), row.get("feature_id", "")) in selected_feature_keys]
    geographic_html = """<!doctype html>
<html><head><meta charset="utf-8"><title>Geographic Distribution</title>
<style>
body { font-family: Arial, sans-serif; color: #1f2933; margin: 1.5rem; }
.controls { display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 0.75rem; align-items: end; }
label { font-weight: 700; display: block; margin-bottom: 0.25rem; }
select { width: 100%; padding: 0.35rem; box-sizing: border-box; }
#map svg { max-width: 100%; height: auto; border: 1px solid #d9e2ec; background: white; }
.warning { background: #fff7ed; border-left: 4px solid #c2410c; padding: 0.75rem; margin: 1rem 0; }
.cards { display: grid; grid-template-columns: repeat(auto-fit, minmax(145px, 1fr)); gap: 0.65rem; margin: 1rem 0; }
.card { border: 1px solid #d9e2ec; border-radius: 6px; padding: 0.7rem; background: #f8fafc; }
.card span { display: block; color: #52606d; font-size: 0.8rem; }
.card strong { font-size: 1.25rem; }
.bar { display: grid; grid-template-columns: minmax(150px, 260px) 1fr 130px; gap: 0.5rem; align-items: center; margin: 0.3rem 0; }
.bar-track { background: #e2e8f0; height: 18px; border-radius: 3px; overflow: hidden; }
.bar-fill { background: #0f766e; height: 100%; }
.bar.small .bar-fill { background: #94a3b8; }
table { border-collapse: collapse; width: 100%; margin-top: 1rem; font-size: 0.9rem; }
th, td { border: 1px solid #d9e2ec; padding: 0.35rem; text-align: left; }
th { background: #f0f4f8; }
a.button { display: inline-block; padding: 0.45rem 0.65rem; margin: 0.2rem 0.35rem 0.2rem 0; background: #0f766e; color: white; text-decoration: none; border-radius: 4px; }
</style></head><body>
<h1>Geographic Distribution</h1>
<p>This section summarizes where selected databases or features were detected in the analyzed dataset. Percentages always use genome-count denominators.</p>
<div class="warning">Geographic distribution reflects the analyzed dataset and may not represent true regional or global prevalence.</div>
<div class="controls">
<div><label for="database">Database</label><select id="database"></select></div>
<div><label for="mode">Mode</label><select id="mode"><option value="database_burden">Database burden / any feature</option><option value="individual_feature">Individual feature / gene</option></select></div>
<div><label for="feature">Feature</label><select id="feature"></select></div>
<div><label for="geo">Geographic level</label><select id="geo"><option value="country">Country</option><option value="continent">Continent</option><option value="subcontinent">Subcontinent / region</option><option value="country_year">Country + collection year</option></select></div>
<div><label for="metric">Metric</label><select id="metric"><option value="prevalence_percent">Prevalence %</option><option value="positive_genomes">Positive genome count</option><option value="total_genomes">Total genome count</option><option value="mean_feature_burden_per_genome">Mean feature burden per genome</option><option value="median_feature_burden_per_genome">Median feature burden per genome</option><option value="feature_rows">Feature row count</option></select></div>
<div><label for="minn">Minimum group size</label><select id="minn"><option value="5">n&gt;=5</option><option value="0">All</option><option value="3">n&gt;=3</option><option value="10">n&gt;=10</option></select></div>
<div><label for="display">Display</label><select id="display"><option value="20">Top 20</option><option value="10">Top 10</option><option value="50">Top 50</option><option value="999999">Complete</option></select></div>
<div><label for="warnings">Warning filter</label><select id="warnings"><option value="all">Show all</option><option value="hide_small">Hide small groups</option><option value="no_major">No major warnings</option></select></div>
</div>
<div id="summary"></div>
<div id="map"></div>
<h2>Ranked groups</h2>
<div id="bars"></div>
<h2>Preview table</h2>
<div id="table"></div>
<p><a class="button" href="../geographic_tables.zip">Download geographic tables ZIP</a><a class="button" href="../geographic_figures.zip">Download geographic figures ZIP</a><a class="button" href="../tables/geographic_database_burden.tsv">Database burden table</a><a class="button" href="../tables/geographic_feature_distribution.tsv">Feature distribution table</a><a class="button" href="../tables/geographic_warning_summary.tsv">Warning summary</a></p>
<script>
const databaseRows = __DATABASE_ROWS__;
const featureRows = __FEATURE_ROWS__;
const coords = __COUNTRY_COORDS__;
const width = 960, height = 480;
function cleanCountry(value) { return (value || '').split(':')[0].trim(); }
function xy(country) {
  const item = coords[cleanCountry(country).toLowerCase()];
  if (!item) return null;
  const lat = item[0], lon = item[1];
  return [((lon + 180) / 360) * width, ((90 - lat) / 180) * height];
}
const databaseSelect = document.getElementById('database');
const modeSelect = document.getElementById('mode');
const featureSelect = document.getElementById('feature');
const geoSelect = document.getElementById('geo');
const metricSelect = document.getElementById('metric');
const minSelect = document.getElementById('minn');
const displaySelect = document.getElementById('display');
const warningSelect = document.getElementById('warnings');
for (const value of [...new Set(databaseRows.map(r => r.database).concat(featureRows.map(r => r.database)))].sort()) {
  const opt = document.createElement('option'); opt.value = value; opt.textContent = value; databaseSelect.appendChild(opt);
}
if ([...databaseSelect.options].some(o => o.value === 'amr')) databaseSelect.value = 'amr';
function updateFeatures() {
  const db = databaseSelect.value;
  const current = featureSelect.value;
  featureSelect.innerHTML = '';
  const features = [...new Set(featureRows.filter(r => r.database === db).map(r => r.feature_id))].sort();
  for (const value of features) { const opt = document.createElement('option'); opt.value = value; opt.textContent = value; featureSelect.appendChild(opt); }
  if (features.includes(current)) featureSelect.value = current;
  featureSelect.disabled = modeSelect.value !== 'individual_feature' || features.length === 0;
}
function rowValue(row, metric) {
  const fallback = metric === 'feature_rows' ? row.total_feature_rows : 0;
  return Number(row[metric] || fallback || 0);
}
function activeRows() {
  const db = databaseSelect.value, mode = modeSelect.value, geo = geoSelect.value, minN = Number(minSelect.value || 0);
  let rows = mode === 'database_burden'
    ? databaseRows.filter(r => r.database === db && r.geo_level === geo)
    : featureRows.filter(r => r.database === db && r.feature_id === featureSelect.value && r.geo_level === geo);
  rows = rows.filter(r => Number(r.total_genomes || 0) >= minN || minN === 0);
  if (warningSelect.value === 'hide_small') rows = rows.filter(r => !(r.warning_flags || '').includes('small_group_warning'));
  if (warningSelect.value === 'no_major') rows = rows.filter(r => !/(small_group_warning|bioproject_dominance|lineage_dominance|single_country_dominance)/.test(r.warning_flags || ''));
  rows = rows.filter(r => !['missing','unknown','missing (unknown)'].includes(r.group_name || ''));
  const metric = metricSelect.value;
  rows.sort((a, b) => rowValue(b, metric) - rowValue(a, metric) || (a.group_name || '').localeCompare(b.group_name || ''));
  return rows.slice(0, Number(displaySelect.value || 20));
}
function renderMap(active) {
  const geo = geoSelect.value;
  if (geo !== 'country' && geo !== 'country_year') {
    document.getElementById('map').innerHTML = '<p>Map view is available for country-level rows. Use the bar plots for continent and region summaries.</p>';
    return;
  }
  let svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height + 45}" viewBox="0 0 ${width} ${height + 45}">`;
  svg += `<rect width="100%" height="100%" fill="#f8fafc"/><text x="20" y="28" font-size="20" font-family="Arial" font-weight="700" fill="#102a43">Geographic Distribution</text>`;
  svg += `<g transform="translate(0,45)"><rect x="0" y="0" width="${width}" height="${height}" fill="#eff6ff" stroke="#bcccdc"/>`;
  for (let lon = -120; lon <= 180; lon += 60) { const x = ((lon + 180) / 360) * width; svg += `<line x1="${x}" y1="0" x2="${x}" y2="${height}" stroke="#d9e2ec"/>`; }
  for (let lat = -60; lat <= 90; lat += 30) { const y = ((90 - lat) / 180) * height; svg += `<line x1="0" y1="${y}" x2="${width}" y2="${y}" stroke="#d9e2ec"/>`; }
  for (const row of active) {
    const point = xy(row.country); if (!point) continue;
    const total = Number(row.total_genomes || 0);
    const radius = Math.max(5, Math.min(28, 4 + Math.sqrt(Math.max(total, 1)) * 4));
    const p = Number(row.prevalence_percent || 0) / 100;
    const fill = (row.warning_flags || '').includes('small_group_warning') ? '#cbd5e1' : `rgb(${Math.round(230 * p)},80,${Math.round(200 * (1 - p))})`;
    const label = `${row.group_name}: ${row.prevalence_display}; warnings=${row.warning_flags || ''}`;
    svg += `<circle cx="${point[0].toFixed(1)}" cy="${point[1].toFixed(1)}" r="${radius.toFixed(1)}" fill="${fill}" fill-opacity="0.75" stroke="#1f2933"><title>${label}</title></circle>`;
    svg += `<text x="${(point[0] + radius + 3).toFixed(1)}" y="${(point[1] + 4).toFixed(1)}" font-size="11" fill="#1f2933">${row.country}</text>`;
  }
  svg += '</g></svg>';
  document.getElementById('map').innerHTML = svg;
}
function renderBars(active) {
  const metric = metricSelect.value;
  const maxValue = Math.max(...active.map(r => rowValue(r, metric)), 1);
  document.getElementById('bars').innerHTML = active.map(row => {
    const value = rowValue(row, metric);
    const width = Math.max(1, value / maxValue * 100);
    const cls = (row.warning_flags || '').includes('small_group_warning') ? 'bar small' : 'bar';
    const display = metric === 'prevalence_percent' ? row.prevalence_display : value.toFixed(metric.includes('burden') ? 2 : 0);
    return `<div class="${cls}"><div>${row.group_name}</div><div class="bar-track"><div class="bar-fill" style="width:${width}%"></div></div><div>${display}</div></div>`;
  }).join('') || '<p>No rows match the selected filters.</p>';
}
function renderTable(active) {
  const headers = ['database','feature_id','geo_level','group_name','total_genomes','positive_genomes','prevalence_display','mean_feature_burden_per_genome','median_feature_burden_per_genome','warning_flags'];
  let html = '<table><thead><tr>' + headers.map(h => `<th>${h}</th>`).join('') + '</tr></thead><tbody>';
  for (const row of active.slice(0, 50)) html += '<tr>' + headers.map(h => `<td>${row[h] || ''}</td>`).join('') + '</tr>';
  html += '</tbody></table>';
  document.getElementById('table').innerHTML = html;
}
function render() {
  updateFeatures();
  const active = activeRows();
  const db = databaseSelect.value, mode = modeSelect.value;
  const positives = active.reduce((a, r) => a + Number(r.positive_genomes || 0), 0);
  const total = active.reduce((a, r) => a + Number(r.total_genomes || 0), 0);
  const top = active[0] || {};
  document.getElementById('summary').innerHTML = `<div class="cards"><div class="card"><span>Mode</span><strong>${mode === 'database_burden' ? 'Database burden' : 'Feature'}</strong></div><div class="card"><span>Database</span><strong>${db}</strong></div><div class="card"><span>Groups shown</span><strong>${active.length}</strong></div><div class="card"><span>Top group</span><strong>${top.group_name || '-'}</strong></div><div class="card"><span>Top prevalence</span><strong>${top.prevalence_display || '-'}</strong></div></div><p>Displayed totals across shown groups: ${positives}/${total}. Interpret these summaries as dataset-specific, not global prevalence.</p>`;
  renderMap(active);
  renderBars(active);
  renderTable(active);
}
for (const control of [databaseSelect, modeSelect, featureSelect, geoSelect, metricSelect, minSelect, displaySelect, warningSelect]) control.addEventListener('change', render);
updateFeatures(); render();
</script></body></html>
"""
    geographic_html = (
        geographic_html
        .replace("__DATABASE_ROWS__", json.dumps(database_rows))
        .replace("__FEATURE_ROWS__", json.dumps(report_feature_rows))
        .replace("__COUNTRY_COORDS__", json.dumps(COUNTRY_COORDS))
    )
    analysis_html_path = figures / "geographic_distribution.html"
    analysis_html_path.write_text(geographic_html, encoding="utf-8")
    html_path = figures / "geographic_distribution_map.html"
    html_path.write_text(geographic_html, encoding="utf-8")

    table_zip = important_dir / "geographic_tables.zip"
    figure_zip = important_dir / "geographic_figures.zip"
    table_files = [summary_path, feature_path, burden_path, warning_path, data_path]
    figure_files += [analysis_html_path, html_path, svg_path, png_path, figures / "geographic_distribution.data.tsv"]
    tables_zip = _write_zip_bundle(table_zip, table_files, important_dir)
    figures_zip = _write_zip_bundle(figure_zip, figure_files, important_dir)
    return {
        "important_geographic_distribution": str(data_path),
        "important_geographic_summary": str(summary_path),
        "important_geographic_feature_distribution": str(feature_path),
        "important_geographic_database_burden": str(burden_path),
        "important_geographic_warning_summary": str(warning_path),
        "important_geographic_analysis_html": str(analysis_html_path),
        "important_geographic_map_html": str(html_path),
        "important_geographic_map_svg": str(svg_path),
        "important_geographic_map_png": str(png_path),
        "important_geographic_tables_zip": tables_zip,
        "important_geographic_figures_zip": figures_zip,
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


def _metadata_support_label(group_n: int, outside_n: int) -> str:
    minimum = min(group_n, outside_n)
    if minimum < 5:
        return "insufficient_support"
    if minimum < 10:
        return "descriptive_only"
    if minimum < 30:
        return "exploratory"
    if minimum >= 100:
        return "strong_support"
    return "standard_support"


def _mean(values: list[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def _median(values: list[float]) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2.0


def _quantile(values: list[float], fraction: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * fraction
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def _format_pvalue(value: float | None) -> str:
    if value is None:
        return ""
    return f"{min(max(value, 0.0), 1.0):.6g}"


def _odds_ratio(a: int, b: int, c: int, d: int) -> float:
    return ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))


def _chi_square_2x2_p_value(a: int, b: int, c: int, d: int) -> float:
    n_total = a + b + c + d
    denominator = (a + b) * (c + d) * (a + c) * (b + d)
    if n_total == 0 or denominator == 0:
        return 1.0
    chi_square = n_total * ((a * d - b * c) ** 2) / denominator
    return math.erfc(math.sqrt(max(chi_square, 0.0) / 2.0))


def _binary_test_for_counts(a: int, b: int, c: int, d: int) -> tuple[str, float]:
    n_total = a + b + c + d
    expected = []
    for row_total in [a + b, c + d]:
        for column_total in [a + c, b + d]:
            expected.append(row_total * column_total / n_total if n_total else 0.0)
    if n_total >= 40 and min(expected or [0.0]) >= 5:
        return "chi_square_2x2", _chi_square_2x2_p_value(a, b, c, d)
    return "fisher_exact", fisher_exact_two_sided(a, b, c, d)


def _rank_values(values: list[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda item: item[1])
    ranks = [0.0] * len(values)
    idx = 0
    while idx < len(indexed):
        end = idx + 1
        while end < len(indexed) and indexed[end][1] == indexed[idx][1]:
            end += 1
        average_rank = (idx + 1 + end) / 2.0
        for original_index, _ in indexed[idx:end]:
            ranks[original_index] = average_rank
        idx = end
    return ranks


def _mann_whitney_u_p_value(group_values: list[float], outside_values: list[float]) -> float | None:
    n1 = len(group_values)
    n2 = len(outside_values)
    if n1 < 1 or n2 < 1:
        return None
    combined = group_values + outside_values
    ranks = _rank_values(combined)
    rank_sum_group = sum(ranks[:n1])
    u1 = rank_sum_group - n1 * (n1 + 1) / 2.0
    mean_u = n1 * n2 / 2.0
    sd_u = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
    if sd_u == 0:
        return 1.0
    z_score = (u1 - mean_u) / sd_u
    return math.erfc(abs(z_score) / math.sqrt(2.0))


def _regularized_gamma_q(a: float, x: float) -> float:
    if a <= 0:
        return 1.0
    if x <= 0:
        return 1.0
    gln = math.lgamma(a)
    if x < a + 1.0:
        ap = a
        delta = 1.0 / a
        total = delta
        for _ in range(100):
            ap += 1.0
            delta *= x / ap
            total += delta
            if abs(delta) < abs(total) * 1e-12:
                break
        p_value = total * math.exp(-x + a * math.log(x) - gln)
        return min(max(1.0 - p_value, 0.0), 1.0)
    tiny = 1e-300
    b = x + 1.0 - a
    c = 1.0 / tiny
    d = 1.0 / max(b, tiny)
    h = d
    for i in range(1, 101):
        an = -float(i) * (float(i) - a)
        b += 2.0
        d = an * d + b
        if abs(d) < tiny:
            d = tiny
        c = b + an / c
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-12:
            break
    q_value = math.exp(-x + a * math.log(x) - gln) * h
    return min(max(q_value, 0.0), 1.0)


def _chi_square_survival(value: float, degrees_freedom: int) -> float:
    if degrees_freedom <= 0:
        return 1.0
    if value <= 0:
        return 1.0
    return _regularized_gamma_q(degrees_freedom / 2.0, value / 2.0)


def _kruskal_wallis(groups: list[list[float]]) -> tuple[float, float] | None:
    groups = [group for group in groups if group]
    if len(groups) < 2:
        return None
    combined = [value for group in groups for value in group]
    n_total = len(combined)
    if n_total <= len(groups):
        return None
    ranks = _rank_values(combined)
    offset = 0
    rank_sums = []
    for group in groups:
        rank_sums.append(sum(ranks[offset:offset + len(group)]))
        offset += len(group)
    statistic = (12.0 / (n_total * (n_total + 1.0))) * sum((rank_sum ** 2) / len(group) for rank_sum, group in zip(rank_sums, groups)) - 3.0 * (n_total + 1.0)
    tie_counts = Counter(combined)
    tie_sum = sum(count ** 3 - count for count in tie_counts.values() if count > 1)
    if tie_sum:
        correction = 1.0 - tie_sum / (n_total ** 3 - n_total)
        if correction > 0:
            statistic /= correction
    statistic = max(statistic, 0.0)
    return statistic, _chi_square_survival(statistic, len(groups) - 1)


def _metadata_support_label_for_counts(counts: list[int]) -> str:
    if not counts:
        return "insufficient_support"
    minimum = min(counts)
    if minimum < 5:
        return "insufficient_support"
    if minimum < 10:
        return "descriptive_only"
    if minimum < 30:
        return "exploratory"
    if minimum >= 100:
        return "strong_support"
    return "standard_support"


def _dominance_flag(samples: set[str], metadata_by_sample: dict[str, dict[str, str]], column: str, flag: str) -> str:
    values = [metadata_by_sample.get(sample, {}).get(column, "") for sample in samples]
    values = [value for value in values if value and not is_missing_value(value)]
    if len(values) < 3:
        return ""
    most_common = Counter(values).most_common(1)[0]
    return flag if most_common[1] / len(values) >= 0.7 else ""


def _metadata_interpretation_label(q_value: float | None, group_n: int, outside_n: int, effect_size: float, warning_flags: str) -> str:
    support = _metadata_support_label(group_n, outside_n)
    severe_flags = {"bioproject_dominance", "lineage_dominance", "low_positive_count", "missing_metadata"}
    flags = {flag for flag in warning_flags.split(";") if flag}
    if support == "insufficient_support":
        return "insufficient_support"
    if support == "descriptive_only":
        return "descriptive_only"
    if support == "exploratory":
        return "exploratory"
    if q_value is not None and q_value <= 0.05 and abs(effect_size) >= 0.20 and not flags.intersection(severe_flags):
        return "strong_supported"
    if q_value is not None and q_value <= 0.10 and abs(effect_size) >= 0.10 and "bioproject_dominance" not in flags:
        return "moderate_supported"
    if q_value is not None and q_value <= 0.10:
        return "exploratory"
    return "descriptive_only"


def _write_metadata_volcano_svg(path: Path, rows: list[dict[str, str]], title: str) -> None:
    width, height = 920, 520
    left, top, plot_w, plot_h = 90, 70, 760, 360
    points = []
    for row in rows:
        x = _float_or_none(row.get("prevalence_difference", "")) or 0.0
        q_value = _float_or_none(row.get("q_value", ""))
        y = -math.log10(max(q_value or 1.0, 1e-12)) if q_value is not None else 0.0
        points.append((x, y, row))
    max_abs_x = max([abs(point[0]) for point in points] + [0.25])
    max_y = max([point[1] for point in points] + [1.0])
    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='#f8fafc'/>",
        f"<text x='20' y='32' font-family='Arial' font-size='20' font-weight='700' fill='#102a43'>{html.escape(title)}</text>",
        f"<line x1='{left}' y1='{top + plot_h}' x2='{left + plot_w}' y2='{top + plot_h}' stroke='#94a3b8'/>",
        f"<line x1='{left + plot_w / 2}' y1='{top}' x2='{left + plot_w / 2}' y2='{top + plot_h}' stroke='#cbd5e1' stroke-dasharray='4 4'/>",
        f"<line x1='{left}' y1='{top}' x2='{left}' y2='{top + plot_h}' stroke='#94a3b8'/>",
    ]
    for x, y, row in points[:600]:
        cx = left + (x + max_abs_x) / (2 * max_abs_x) * plot_w
        cy = top + plot_h - (y / max_y) * plot_h if max_y else top + plot_h
        label = row.get("interpretation_label", "")
        color = "#94a3b8"
        if label in {"strong_supported", "moderate_supported"} and x > 0:
            color = "#dc2626"
        elif label in {"strong_supported", "moderate_supported"} and x < 0:
            color = "#2563eb"
        tip = (
            f"{row.get('feature_id', '')}; diff={x:.3f}; q={row.get('q_value', '')}; "
            f"{row.get('positive_in_group', '')}/{row.get('group_n', '')} vs "
            f"{row.get('positive_outside_group', '')}/{row.get('outside_group_n', '')}; {row.get('warning_flags', '')}"
        )
        parts.append(f"<circle cx='{cx:.1f}' cy='{cy:.1f}' r='4' fill='{color}' fill-opacity='0.82'><title>{html.escape(tip)}</title></circle>")
    parts.extend([
        f"<text x='{left + plot_w / 2 - 90}' y='{height - 35}' font-family='Arial' font-size='12' fill='#52606d'>Prevalence difference (selected group - outside)</text>",
        f"<text x='18' y='{top + 20}' font-family='Arial' font-size='12' fill='#52606d' transform='rotate(-90 18,{top + 20})'>-log10(q-value)</text>",
        "</svg>\n",
    ])
    path.write_text("".join(parts), encoding="utf-8")


def _write_metadata_volcano_png(path: Path, rows: list[dict[str, str]]) -> None:
    width, height = 920, 520
    left, top, plot_w, plot_h = 90, 70, 760, 360
    pixels = [bytearray([248, 250, 252] * width) for _ in range(height)]
    points = []
    for row in rows:
        x = _float_or_none(row.get("prevalence_difference", "")) or 0.0
        q_value = _float_or_none(row.get("q_value", ""))
        y = -math.log10(max(q_value or 1.0, 1e-12)) if q_value is not None else 0.0
        points.append((x, y, row))
    max_abs_x = max([abs(point[0]) for point in points] + [0.25])
    max_y = max([point[1] for point in points] + [1.0])

    def set_pixel(x: int, y: int, color: tuple[int, int, int]) -> None:
        if 0 <= x < width and 0 <= y < height:
            idx = x * 3
            pixels[y][idx:idx + 3] = bytes(color)

    def circle(cx: int, cy: int, radius: int, color: tuple[int, int, int]) -> None:
        for yy in range(cy - radius, cy + radius + 1):
            for xx in range(cx - radius, cx + radius + 1):
                if (xx - cx) ** 2 + (yy - cy) ** 2 <= radius ** 2:
                    set_pixel(xx, yy, color)

    for x, y, row in points[:600]:
        cx = int(left + (x + max_abs_x) / (2 * max_abs_x) * plot_w)
        cy = int(top + plot_h - (y / max_y) * plot_h) if max_y else top + plot_h
        label = row.get("interpretation_label", "")
        color = (148, 163, 184)
        if label in {"strong_supported", "moderate_supported"} and x > 0:
            color = (220, 38, 38)
        elif label in {"strong_supported", "moderate_supported"} and x < 0:
            color = (37, 99, 235)
        circle(cx, cy, 4, color)
    _write_png(path, width, height, pixels)


def _write_diverging_heatmap_svg(path: Path, rows: list[dict[str, str]], title: str, row_field: str, column_field: str, value_field: str) -> None:
    row_labels = list(dict.fromkeys(row.get(row_field, "") for row in rows if row.get(row_field, "")))
    column_labels = sorted({row.get(column_field, "") for row in rows if row.get(column_field, "")})
    values = {(row.get(row_field, ""), row.get(column_field, "")): _float_or_none(row.get(value_field, "")) or 0.0 for row in rows}
    max_abs = max([abs(value) for value in values.values()] + [0.1])
    cell_w, cell_h = 92, 24
    left, top = 260, 84
    width = max(840, left + cell_w * max(len(column_labels), 1) + 60)
    height = max(220, top + cell_h * max(len(row_labels), 1) + 50)
    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='#f8fafc'/>",
        f"<text x='20' y='32' font-family='Arial' font-size='20' font-weight='700' fill='#102a43'>{html.escape(title)}</text>",
    ]
    for cidx, column in enumerate(column_labels):
        x = left + cidx * cell_w
        short = column if len(column) <= 18 else column[:15] + "..."
        parts.append(f"<text x='{x + 4}' y='{top - 14}' font-family='Arial' font-size='11' fill='#52606d'>{html.escape(short)}</text>")
    for ridx, label in enumerate(row_labels):
        y = top + ridx * cell_h
        short_label = label if len(label) <= 40 else label[:37] + "..."
        parts.append(f"<text x='20' y='{y + 16}' font-family='Arial' font-size='11' fill='#1f2933'>{html.escape(short_label)}</text>")
        for cidx, column in enumerate(column_labels):
            x = left + cidx * cell_w
            value = values.get((label, column), 0.0)
            intensity = min(abs(value) / max_abs, 1.0)
            if value > 0:
                fill = f"rgb(254,{int(226 - 130 * intensity)},{int(226 - 130 * intensity)})"
            elif value < 0:
                fill = f"rgb({int(219 - 135 * intensity)},{int(234 - 120 * intensity)},254)"
            else:
                fill = "#f1f5f9"
            parts.append(f"<rect x='{x}' y='{y}' width='{cell_w - 2}' height='{cell_h - 2}' fill='{fill}' stroke='#d9e2ec'><title>{html.escape(label)} / {html.escape(column)}: {value:.3f}</title></rect>")
    parts.append("</svg>\n")
    path.write_text("".join(parts), encoding="utf-8")


def _write_diverging_heatmap_png(path: Path, rows: list[dict[str, str]], row_field: str, column_field: str, value_field: str) -> None:
    row_labels = list(dict.fromkeys(row.get(row_field, "") for row in rows if row.get(row_field, "")))
    column_labels = sorted({row.get(column_field, "") for row in rows if row.get(column_field, "")})
    values = {(row.get(row_field, ""), row.get(column_field, "")): _float_or_none(row.get(value_field, "")) or 0.0 for row in rows}
    max_abs = max([abs(value) for value in values.values()] + [0.1])
    cell_w, cell_h = 92, 24
    left, top = 260, 84
    width = max(840, left + cell_w * max(len(column_labels), 1) + 60)
    height = max(220, top + cell_h * max(len(row_labels), 1) + 50)
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
            value = values.get((label, column), 0.0)
            intensity = min(abs(value) / max_abs, 1.0)
            if value > 0:
                color = (254, int(226 - 130 * intensity), int(226 - 130 * intensity))
            elif value < 0:
                color = (int(219 - 135 * intensity), int(234 - 120 * intensity), 254)
            else:
                color = (241, 245, 249)
            rect(x, y, x + cell_w - 2, y + cell_h - 2, color)
    _write_png(path, width, height, pixels)


def _write_burden_boxplot_svg(path: Path, rows: list[dict[str, str]], title: str) -> None:
    width = 940
    row_h = 42
    left, top, plot_w = 240, 64, 620
    height = max(220, top + row_h * max(len(rows), 1) + 50)
    max_value = max([_float_or_none(row.get("max", "")) or 0.0 for row in rows] + [1.0])
    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='#f8fafc'/>",
        f"<text x='20' y='32' font-family='Arial' font-size='20' font-weight='700' fill='#102a43'>{html.escape(title)}</text>",
    ]
    for idx, row in enumerate(rows):
        y = top + idx * row_h
        label = row.get("metadata_group", "")
        short = label if len(label) <= 32 else label[:29] + "..."
        values = {key: (_float_or_none(row.get(key, "")) or 0.0) for key in ["min", "q1", "median", "q3", "max", "mean"]}
        def xpos(value: float) -> float:
            return left + value / max_value * plot_w if max_value else left
        parts.append(f"<text x='20' y='{y + 23}' font-family='Arial' font-size='12' fill='#1f2933'>{html.escape(short)} (n={html.escape(row.get('n', ''))})</text>")
        parts.append(f"<line x1='{xpos(values['min']):.1f}' y1='{y + 18}' x2='{xpos(values['max']):.1f}' y2='{y + 18}' stroke='#64748b'/>")
        parts.append(f"<rect x='{xpos(values['q1']):.1f}' y='{y + 8}' width='{max(xpos(values['q3']) - xpos(values['q1']), 2):.1f}' height='20' fill='#bae6fd' stroke='#0369a1'/>")
        parts.append(f"<line x1='{xpos(values['median']):.1f}' y1='{y + 6}' x2='{xpos(values['median']):.1f}' y2='{y + 30}' stroke='#0f172a' stroke-width='2'/>")
        parts.append(f"<circle cx='{xpos(values['mean']):.1f}' cy='{y + 18}' r='4' fill='#0f766e'><title>mean {values['mean']:.2f}</title></circle>")
    parts.append("</svg>\n")
    path.write_text("".join(parts), encoding="utf-8")


def _write_burden_boxplot_png(path: Path, rows: list[dict[str, str]]) -> None:
    # Dependency-free raster companion. The SVG carries labels and tooltips.
    _write_bar_png(path, [{"metadata_group": row.get("metadata_group", ""), "median": row.get("median", "")} for row in rows], "median")


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


def _prevalence_label(percent: float) -> str:
    if percent >= 95:
        return "core"
    if percent >= 50:
        return "common"
    if percent >= 5:
        return "accessory"
    if percent > 0:
        return "rare"
    return "absent"


def _write_prevalence_bar_svg(path: Path, rows: list[dict[str, str]], title: str, value_field: str) -> None:
    width = 1040
    row_height = 31
    top = 58
    left = 285
    plot_width = 610
    height = max(190, top + row_height * max(len(rows), 1) + 44)
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
        label = row.get("feature_id", "") or row.get("database", "") or row.get("label", "")
        if len(label) > 42:
            label = label[:39] + "..."
        display = row.get("prevalence_display", "") or row.get(value_field, "")
        tip = (
            f"{row.get('database', '')}:{row.get('feature_id', '')}; "
            f"{row.get('positive_genomes', '')}/{row.get('total_genomes', row.get('sample_count', ''))}; "
            f"rows={row.get('feature_rows', '')}; label={row.get('prevalence_label', '')}"
        )
        parts.append(f"<text x='20' y='{y + 18}' font-family='Arial' font-size='12' fill='#1f2933'>{html.escape(label)}</text>")
        parts.append(f"<rect x='{left}' y='{y}' width='{bar_width:.1f}' height='20' fill='#0f766e'><title>{html.escape(tip)}</title></rect>")
        parts.append(f"<text x='{left + bar_width + 8:.1f}' y='{y + 15}' font-family='Arial' font-size='12' fill='#1f2933'>{html.escape(display)}</text>")
    parts.append(f"<text x='{left}' y='{height - 12}' font-family='Arial' font-size='12' fill='#52606d'>{html.escape(value_field)}</text>")
    parts.append("</svg>\n")
    path.write_text("".join(parts), encoding="utf-8")


def _write_prevalence_stacked_svg(path: Path, rows: list[dict[str, str]], title: str) -> None:
    width = 980
    row_height = 34
    top, left, plot_width = 76, 210, 620
    height = max(220, top + row_height * max(len(rows), 1) + 50)
    colors = {
        "core_features": "#0f766e",
        "common_features": "#38bdf8",
        "accessory_features": "#f59e0b",
        "rare_features": "#94a3b8",
    }
    labels = [
        ("core_features", "core"),
        ("common_features", "common"),
        ("accessory_features", "accessory"),
        ("rare_features", "rare"),
    ]
    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='#f8fafc'/>",
        f"<text x='20' y='30' font-family='Arial' font-size='20' font-weight='700' fill='#102a43'>{html.escape(title)}</text>",
    ]
    legend_x = 20
    for key, label in labels:
        parts.append(f"<rect x='{legend_x}' y='48' width='12' height='12' fill='{colors[key]}'/>")
        parts.append(f"<text x='{legend_x + 18}' y='59' font-family='Arial' font-size='12' fill='#334e68'>{label}</text>")
        legend_x += 112
    for idx, row in enumerate(rows):
        y = top + idx * row_height
        total = _float_or_none(row.get("total_unique_features", "")) or 0.0
        if total <= 0:
            total = sum((_float_or_none(row.get(key, "")) or 0.0) for key, _ in labels) or 1.0
        x = left
        database = row.get("database", "")
        parts.append(f"<text x='20' y='{y + 20}' font-family='Arial' font-size='12' fill='#1f2933'>{html.escape(database)}</text>")
        for key, label in labels:
            value = _float_or_none(row.get(key, "")) or 0.0
            width_part = value / total * plot_width
            parts.append(f"<rect x='{x:.1f}' y='{y}' width='{width_part:.1f}' height='22' fill='{colors[key]}'><title>{html.escape(database)} {label}: {value:g}</title></rect>")
            x += width_part
        parts.append(f"<text x='{left + plot_width + 10}' y='{y + 16}' font-family='Arial' font-size='12' fill='#1f2933'>{int(total)}</text>")
    parts.append("</svg>\n")
    path.write_text("".join(parts), encoding="utf-8")


def _write_prevalence_stacked_png(path: Path, rows: list[dict[str, str]]) -> None:
    width, height = 980, max(220, 76 + 34 * max(len(rows), 1) + 50)
    top, left, plot_width = 76, 210, 620
    colors = {
        "core_features": (15, 118, 110),
        "common_features": (56, 189, 248),
        "accessory_features": (245, 158, 11),
        "rare_features": (148, 163, 184),
    }
    pixels = [bytearray([248, 250, 252] * width) for _ in range(height)]

    def rect(x0: int, y0: int, x1: int, y1: int, color: tuple[int, int, int]) -> None:
        for y in range(max(0, y0), min(height, y1)):
            row_pixels = pixels[y]
            for x in range(max(0, x0), min(width, x1)):
                idx = x * 3
                row_pixels[idx:idx + 3] = bytes(color)

    for idx, row in enumerate(rows):
        y = top + idx * 34
        total = _float_or_none(row.get("total_unique_features", "")) or 0.0
        if total <= 0:
            total = sum((_float_or_none(row.get(key, "")) or 0.0) for key in colors) or 1.0
        x = left
        for key, color in colors.items():
            value = _float_or_none(row.get(key, "")) or 0.0
            width_part = int(value / total * plot_width)
            rect(int(x), y, int(x) + width_part, y + 22, color)
            x += width_part
    _write_png(path, width, height, pixels)


def write_important_prevalence_outputs(sample_dir: Path, out_dir: Path, important_dir: Path, top_n: int = 20) -> dict[str, str]:
    key_tables = important_dir / "key_tables"
    tables = important_dir / "tables"
    figures = important_dir / "figures"
    key_tables.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    metadata_rows = read_table(sample_dir / "basic" / "enriched_genome_dataset.csv")
    metadata_samples = [row.get("assembly_accession", "") or row.get("sample_id", "") for row in metadata_rows]
    features = [row for row in read_table(out_dir / "features" / "all_features.tsv") if row.get("presence", "1") != "0"]
    feature_samples = [row.get("assembly_accession", "") or row.get("sample_id", "") for row in features]
    samples = sorted({sample for sample in metadata_samples + feature_samples if sample})
    sample_count = len(samples)
    metadata_by_sample = {
        row.get("assembly_accession", "") or row.get("sample_id", ""): row
        for row in metadata_rows
        if row.get("assembly_accession", "") or row.get("sample_id", "")
    }
    by_feature: dict[tuple[str, str], set[str]] = defaultdict(set)
    rows_by_feature: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    sample_database_rows: Counter[tuple[str, str]] = Counter()
    sample_database_features: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in features:
        key = (row.get("database", ""), row.get("feature_id", ""))
        sample = row.get("assembly_accession", "") or row.get("sample_id", "")
        if key[0] and key[1] and sample:
            by_feature[key].add(sample)
            rows_by_feature[key].append(row)
            sample_database_rows[(sample, key[0])] += 1
            sample_database_features[(sample, key[0])].add(key[1])

    feature_rows = []
    for (database, feature_id), samples in sorted(by_feature.items()):
        rows = rows_by_feature[(database, feature_id)]
        category = first_value(rows[0], ["feature_category"], "")
        subcategory = first_value(rows[0], ["feature_subcategory"], "")
        identities = [_float_or_none(row.get("identity", "")) for row in rows]
        coverages = [_float_or_none(row.get("coverage", "")) for row in rows]
        identities = [value for value in identities if value is not None]
        coverages = [value for value in coverages if value is not None]
        identity_stats = _summary_stats(identities)
        coverage_stats = _summary_stats(coverages)
        prevalence = len(samples) / sample_count if sample_count else 0.0
        prevalence_percent = prevalence * 100
        mean_hits = len(rows) / len(samples) if samples else 0.0
        has_coordinates = any(row.get("contig") and row.get("start") and row.get("end") for row in rows)
        tools_detected = sorted({row.get("tool", "") for row in rows if row.get("tool", "")})
        database_versions = sorted({row.get("database_version", "") for row in rows if row.get("database_version", "")})
        label = _prevalence_label(prevalence_percent)
        warnings = []
        if prevalence_percent < 5:
            warnings.append("low_prevalence")
        if len(samples) == 1:
            warnings.append("single_genome_feature")
        if not category:
            warnings.append("missing_category_annotation")
        if len(rows) > len(samples):
            warnings.append("duplicate_feature_rows_detected")
        row_out = {
            "database": database,
            "feature_id": feature_id,
            "feature_name": first_value(rows[0], ["feature_name", "feature_id"], feature_id),
            "feature_category": category,
            "feature_subcategory": subcategory,
            "positive_genomes": str(len(samples)),
            "total_genomes": str(sample_count),
            "sample_count": str(sample_count),
            "prevalence": f"{prevalence:.4f}",
            "prevalence_percent": f"{prevalence_percent:.1f}",
            "prevalence_display": f"{prevalence_percent:.1f}% ({len(samples)}/{sample_count})",
            "feature_rows": str(len(rows)),
            "mean_hits_per_positive_genome": f"{mean_hits:.2f}" if samples else "",
            "median_identity": identity_stats["median"],
            "median_coverage": coverage_stats["median"],
            "has_coordinates": str(has_coordinates).lower(),
            "tools_detected": ";".join(tools_detected),
            "database_version": ";".join(database_versions),
            "prevalence_label": label,
            "warning_flags": ";".join(warnings),
        }
        feature_rows.append(row_out)
    feature_rows.sort(key=lambda row: (row["database"], -(_float_or_none(row["prevalence_percent"]) or 0.0), row["feature_id"]))
    feature_fields = [
        "database", "feature_id", "feature_name", "feature_category", "feature_subcategory",
        "positive_genomes", "total_genomes", "sample_count", "prevalence", "prevalence_percent", "prevalence_display",
        "feature_rows", "mean_hits_per_positive_genome", "median_identity", "median_coverage",
        "has_coordinates", "tools_detected", "database_version", "prevalence_label", "warning_flags",
    ]

    burden_rows = []
    databases = sorted({row["database"] for row in feature_rows})
    for database in databases:
        for sample in samples:
            meta = metadata_by_sample.get(sample, {})
            unique_features = len(sample_database_features.get((sample, database), set()))
            row_count = sample_database_rows.get((sample, database), 0)
            burden_rows.append({
                "database": database,
                "assembly_accession": sample,
                "sample_id": meta.get("sample_id", ""),
                "feature_rows": str(row_count),
                "unique_features": str(unique_features),
                "has_feature": str(unique_features > 0).lower(),
            })

    summary_by_database = []
    core_rows = []
    written_rows = []
    for database in databases:
        db_features = [row for row in feature_rows if row["database"] == database]
        db_burden = [row for row in burden_rows if row["database"] == database]
        unique_counts = [int(_float_or_none(row.get("unique_features", "")) or 0) for row in db_burden]
        unique_stats = _summary_stats_full([float(value) for value in unique_counts])
        positive_genomes = sum(1 for row in db_burden if row.get("has_feature") == "true")
        total_rows = sum(int(_float_or_none(row.get("feature_rows", "")) or 0) for row in db_features)
        top_feature = db_features[0] if db_features else {}
        positive_percent = positive_genomes / sample_count * 100 if sample_count else 0.0
        summary_by_database.append({
            "database": database,
            "total_feature_rows": str(total_rows),
            "unique_features": str(len(db_features)),
            "positive_genomes": str(positive_genomes),
            "total_genomes": str(sample_count),
            "genomes_positive_percent": f"{positive_percent:.1f}",
            "median_features_per_genome": "" if unique_stats["median"] is None else f"{unique_stats['median']:.2f}",
            "mean_features_per_genome": "" if unique_stats["mean"] is None else f"{unique_stats['mean']:.2f}",
            "max_features_per_genome": "" if unique_stats["max"] is None else f"{unique_stats['max']:.0f}",
            "top_feature_id": top_feature.get("feature_id", ""),
            "top_feature_prevalence_percent": top_feature.get("prevalence_percent", ""),
            "top_feature_positive_genomes": top_feature.get("positive_genomes", ""),
        })
        label_counts = Counter(row.get("prevalence_label", "") for row in db_features)
        core_rows.append({
            "database": database,
            "core_features": str(label_counts.get("core", 0)),
            "common_features": str(label_counts.get("common", 0)),
            "accessory_features": str(label_counts.get("accessory", 0)),
            "rare_features": str(label_counts.get("rare", 0)),
            "total_unique_features": str(len(db_features)),
        })
        top_text = ", ".join(
            f"{row.get('feature_id', '')} ({row.get('prevalence_percent', '')}%, {row.get('positive_genomes', '')}/{row.get('total_genomes', '')})"
            for row in db_features[:5]
        ) or "none"
        written_rows.append({
            "scope": "database",
            "database": database,
            "summary": (
                f"{database} features were detected in {positive_genomes} of {sample_count} genomes. "
                f"The run identified {total_rows} feature rows and {len(db_features)} unique features from this database. "
                f"The most prevalent features were {top_text}. "
                f"The median burden was {summary_by_database[-1]['median_features_per_genome']} features per genome."
            ),
        })
    total_feature_rows = sum(int(row["feature_rows"]) for row in feature_rows)
    total_unique_features = len(feature_rows)
    top_databases = sorted(summary_by_database, key=lambda row: -(int(_float_or_none(row.get("total_feature_rows", "")) or 0)))[:5]
    written_rows.insert(0, {
        "scope": "overall",
        "database": "all",
        "summary": (
            f"Across {sample_count} analyzed genomes, PanResistome detected {total_feature_rows} standardized feature rows "
            f"representing {total_unique_features} unique features from {len(databases)} databases. "
            f"The largest contributors by feature rows were {', '.join(row.get('database', '') for row in top_databases) or 'none'}. "
            "Prevalence is descriptive and reflects this dataset only."
        ),
    })

    top_rows = []
    for database in databases:
        top_rows.extend([row for row in feature_rows if row["database"] == database][:top_n])

    summary_path = key_tables / "feature_prevalence_summary.tsv"
    feature_prevalence_path = tables / "feature_prevalence.tsv"
    top_path = tables / "feature_prevalence_top.tsv"
    summary_db_path = tables / "prevalence_summary_by_database.tsv"
    core_path = tables / "prevalence_core_accessory_rare_summary.tsv"
    burden_path = tables / "prevalence_database_burden_by_sample.tsv"
    written_path = tables / "prevalence_written_summaries.tsv"
    write_rows(summary_path, feature_rows, feature_fields)
    write_rows(feature_prevalence_path, feature_rows, feature_fields)
    write_rows(top_path, top_rows, feature_fields)
    summary_db_fields = [
        "database", "total_feature_rows", "unique_features", "positive_genomes", "total_genomes",
        "genomes_positive_percent", "median_features_per_genome", "mean_features_per_genome",
        "max_features_per_genome", "top_feature_id", "top_feature_prevalence_percent", "top_feature_positive_genomes",
    ]
    write_rows(summary_db_path, summary_by_database, summary_db_fields)
    write_rows(core_path, core_rows, ["database", "core_features", "common_features", "accessory_features", "rare_features", "total_unique_features"])
    write_rows(burden_path, burden_rows, ["database", "assembly_accession", "sample_id", "feature_rows", "unique_features", "has_feature"])
    write_rows(written_path, written_rows, ["scope", "database", "summary"])

    figure_files: list[Path] = []
    counts_data = figures / "prevalence_feature_counts_by_database.data.tsv"
    counts_svg = figures / "prevalence_feature_counts_by_database.svg"
    counts_png = figures / "prevalence_feature_counts_by_database.png"
    counts_pdf = figures / "prevalence_feature_counts_by_database.pdf"
    write_rows(counts_data, summary_by_database, summary_db_fields)
    _write_bar_svg(counts_svg, summary_by_database, "Feature Counts By Database", "database", "unique_features", "Unique features")
    _write_bar_png(counts_png, summary_by_database, "unique_features")
    _write_simple_pdf(counts_pdf, "Feature Counts By Database", [f"{row['database']}: {row['unique_features']} unique, {row['total_feature_rows']} rows" for row in summary_by_database])
    figure_files.extend([counts_data, counts_svg, counts_png, counts_pdf])

    genomes_data = figures / "prevalence_genomes_positive_by_database.data.tsv"
    genomes_svg = figures / "prevalence_genomes_positive_by_database.svg"
    genomes_png = figures / "prevalence_genomes_positive_by_database.png"
    genomes_pdf = figures / "prevalence_genomes_positive_by_database.pdf"
    genome_rows = [
        {**row, "prevalence_display": f"{row.get('genomes_positive_percent', '')}% ({row.get('positive_genomes', '')}/{row.get('total_genomes', '')})"}
        for row in summary_by_database
    ]
    write_rows(genomes_data, genome_rows, summary_db_fields + ["prevalence_display"])
    _write_prevalence_bar_svg(genomes_svg, genome_rows, "Genomes Positive By Database", "genomes_positive_percent")
    _write_bar_png(genomes_png, genome_rows, "genomes_positive_percent")
    _write_simple_pdf(genomes_pdf, "Genomes Positive By Database", [f"{row['database']}: {row['prevalence_display']}" for row in genome_rows])
    figure_files.extend([genomes_data, genomes_svg, genomes_png, genomes_pdf])

    core_data = figures / "prevalence_core_accessory_rare_by_database.data.tsv"
    core_svg = figures / "prevalence_core_accessory_rare_by_database.svg"
    core_png = figures / "prevalence_core_accessory_rare_by_database.png"
    core_pdf = figures / "prevalence_core_accessory_rare_by_database.pdf"
    write_rows(core_data, core_rows, ["database", "core_features", "common_features", "accessory_features", "rare_features", "total_unique_features"])
    _write_prevalence_stacked_svg(core_svg, core_rows, "Core/Common/Accessory/Rare Features")
    _write_prevalence_stacked_png(core_png, core_rows)
    _write_simple_pdf(core_pdf, "Core/Common/Accessory/Rare Features", [f"{row['database']}: core={row['core_features']}, common={row['common_features']}, accessory={row['accessory_features']}, rare={row['rare_features']}" for row in core_rows])
    figure_files.extend([core_data, core_svg, core_png, core_pdf])

    burden_plot_rows = []
    for database in databases:
        counts = [float(int(_float_or_none(row.get("unique_features", "")) or 0)) for row in burden_rows if row["database"] == database]
        stats = _summary_stats_full(counts)
        burden_plot_rows.append({
            "metadata_group": database,
            "database": database,
            "n": str(len(counts)),
            "min": "" if stats["min"] is None else f"{stats['min']:.2f}",
            "q1": "" if stats["q1"] is None else f"{stats['q1']:.2f}",
            "median": "" if stats["median"] is None else f"{stats['median']:.2f}",
            "q3": "" if stats["q3"] is None else f"{stats['q3']:.2f}",
            "max": "" if stats["max"] is None else f"{stats['max']:.2f}",
            "mean": "" if stats["mean"] is None else f"{stats['mean']:.2f}",
        })
    burden_data = figures / "prevalence_database_burden_by_sample.data.tsv"
    burden_svg = figures / "prevalence_database_burden_by_sample.svg"
    burden_png = figures / "prevalence_database_burden_by_sample.png"
    burden_pdf = figures / "prevalence_database_burden_by_sample.pdf"
    write_rows(burden_data, burden_plot_rows, ["metadata_group", "database", "n", "min", "q1", "median", "q3", "max", "mean"])
    _write_burden_boxplot_svg(burden_svg, burden_plot_rows, "Database Burden By Sample")
    _write_burden_boxplot_png(burden_png, burden_plot_rows)
    _write_simple_pdf(burden_pdf, "Database Burden By Sample", [f"{row['database']}: median={row['median']}, n={row['n']}" for row in burden_plot_rows])
    figure_files.extend([burden_data, burden_svg, burden_png, burden_pdf])

    outputs = {
        "important_feature_prevalence_summary": str(summary_path),
        "important_feature_prevalence": str(feature_prevalence_path),
        "important_feature_prevalence_top": str(top_path),
        "important_prevalence_summary_by_database": str(summary_db_path),
        "important_prevalence_core_accessory_rare_summary": str(core_path),
        "important_prevalence_database_burden_by_sample": str(burden_path),
        "important_prevalence_written_summaries": str(written_path),
        "important_prevalence_feature_counts_svg": str(counts_svg),
        "important_prevalence_genomes_positive_svg": str(genomes_svg),
        "important_prevalence_core_accessory_rare_svg": str(core_svg),
        "important_prevalence_database_burden_svg": str(burden_svg),
    }
    for database in databases:
        safe_database = _safe_filename(database)
        db_rows = [row for row in feature_rows if row["database"] == database][:top_n]
        if not db_rows:
            continue
        data_path = figures / f"prevalence_top_features_{safe_database}.data.tsv"
        svg_path = figures / f"prevalence_top_features_{safe_database}.svg"
        png_path = figures / f"prevalence_top_features_{safe_database}.png"
        pdf_path = figures / f"prevalence_top_features_{safe_database}.pdf"
        write_rows(data_path, db_rows, feature_fields)
        _write_prevalence_bar_svg(svg_path, db_rows, f"{database} Top Feature Prevalence", "prevalence_percent")
        _write_bar_png(png_path, db_rows, "prevalence_percent")
        _write_simple_pdf(pdf_path, f"{database} Top Feature Prevalence", [f"{row['feature_id']}: {row['prevalence_display']}" for row in db_rows])
        figure_files.extend([data_path, svg_path, png_path, pdf_path])
        outputs[f"important_prevalence_{safe_database}_data"] = str(data_path)
        outputs[f"important_prevalence_{safe_database}_svg"] = str(svg_path)
        outputs[f"important_prevalence_{safe_database}_png"] = str(png_path)
        outputs[f"important_prevalence_{safe_database}_pdf"] = str(pdf_path)

    html_path = figures / "prevalence_analysis.html"
    prevalence_html = """<!doctype html>
<html><head><meta charset="utf-8"><title>Feature Prevalence</title>
<style>
body { font-family: Arial, sans-serif; margin: 1.5rem; color: #1f2933; background: #f8fafc; }
label { font-weight: 700; margin-right: 0.35rem; }
select, input { margin: 0 0.9rem 1rem 0; padding: 0.35rem; }
.warning { background: #fff7ed; border-left: 4px solid #c2410c; padding: 0.75rem; margin: 1rem 0; }
.cards { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 0.75rem; }
.card { border: 1px solid #d9e2ec; border-radius: 6px; padding: 0.75rem; background: white; }
.card span { display: block; color: #52606d; font-size: 0.8rem; }
.card strong { display: block; font-size: 1.25rem; margin-top: 0.25rem; }
.figure-box { max-width: 100%; max-height: 680px; overflow: auto; border: 1px solid #d9e2ec; background: white; padding: 0.75rem; }
svg { max-width: 100%; height: auto; }
table { border-collapse: collapse; width: 100%; margin-top: 1rem; background: white; font-size: 0.9rem; }
th, td { border: 1px solid #d9e2ec; padding: 0.4rem; text-align: left; }
th { background: #f0f4f8; position: sticky; top: 0; }
.downloads a { display: inline-block; margin: 0.25rem 0.45rem 0.25rem 0; color: #0f766e; font-weight: 700; }
</style></head><body>
<h1>Prevalence</h1>
<p>This section summarizes how frequently each detected feature appears across the analyzed genomes. Prevalence uses positive genome count; feature row count may be higher when a feature appears multiple times in one genome.</p>
<div class="warning">Prevalence reflects the analyzed dataset, not global prevalence. It is descriptive and does not test association or causality.</div>
<label for="database">Database</label><select id="database"></select>
<label for="view">View</label><select id="view"><option value="20">Top 20</option><option value="10">Top 10</option><option value="50">Top 50</option><option value="all">Complete</option></select>
<label for="metric">Metric</label><select id="metric"><option value="prevalence_percent">Genome prevalence %</option><option value="positive_genomes">Positive genome count</option><option value="feature_rows">Feature row count</option><option value="mean_hits_per_positive_genome">Mean hits per positive genome</option></select>
<label for="sort">Sort</label><select id="sort"><option value="prevalence">Prevalence descending</option><option value="feature">Feature name</option><option value="category">Feature category</option><option value="database">Database</option><option value="positive">Positive genome count</option><option value="rows">Row count</option></select>
<label for="minPrev">Minimum prevalence %</label><input id="minPrev" type="number" value="0" min="0" max="100" step="1">
<label for="minPositive">Minimum positive genomes</label><input id="minPositive" type="number" value="0" min="0" step="1">
<label><input id="coordsOnly" type="checkbox"> Show only features with coordinates</label>
<label><input id="categoryOnly" type="checkbox"> Show only features with category annotation</label>
<div id="summary"></div>
<p id="written"></p>
<div class="figure-box" id="figure"></div>
<div id="table"></div>
<div class="downloads">
  <a href="../tables/feature_prevalence.tsv">Download full feature prevalence</a>
  <a href="../tables/feature_prevalence_top.tsv">Download report-facing top features</a>
  <a href="../tables/prevalence_summary_by_database.tsv">Download database summary</a>
  <a href="../prevalence_tables.zip">Download prevalence tables ZIP</a>
  <a href="../prevalence_figures.zip">Download prevalence figures ZIP</a>
</div>
<script>
const rows = __FEATURE_ROWS__;
const databaseRows = __DATABASE_ROWS__;
const writtenRows = __WRITTEN_ROWS__;
function num(value) { const parsed = Number(value || 0); return Number.isFinite(parsed) ? parsed : 0; }
const dbSelect = document.getElementById('database');
for (const db of ['All', ...[...new Set(rows.map(r => r.database).filter(Boolean))].sort()]) {
  const option = document.createElement('option'); option.value = db; option.textContent = db; dbSelect.appendChild(option);
}
function activeRows() {
  const db = dbSelect.value, metric = document.getElementById('metric').value;
  const minPrev = num(document.getElementById('minPrev').value), minPositive = num(document.getElementById('minPositive').value);
  const coordsOnly = document.getElementById('coordsOnly').checked, categoryOnly = document.getElementById('categoryOnly').checked;
  let active = rows.filter(row => (db === 'All' || row.database === db) && num(row.prevalence_percent) >= minPrev && num(row.positive_genomes) >= minPositive);
  if (coordsOnly) active = active.filter(row => row.has_coordinates === 'true');
  if (categoryOnly) active = active.filter(row => row.feature_category);
  const sort = document.getElementById('sort').value;
  if (sort === 'feature') active = active.slice().sort((a, b) => (a.feature_id || '').localeCompare(b.feature_id || ''));
  else if (sort === 'category') active = active.slice().sort((a, b) => (a.feature_category || '').localeCompare(b.feature_category || ''));
  else if (sort === 'database') active = active.slice().sort((a, b) => (a.database || '').localeCompare(b.database || ''));
  else if (sort === 'positive') active = active.slice().sort((a, b) => num(b.positive_genomes) - num(a.positive_genomes));
  else if (sort === 'rows') active = active.slice().sort((a, b) => num(b.feature_rows) - num(a.feature_rows));
  else active = active.slice().sort((a, b) => num(b[metric]) - num(a[metric]));
  const view = document.getElementById('view').value;
  if (view === 'all') return active;
  return active.slice(0, Number(view || 20));
}
function renderCards(active) {
  const db = dbSelect.value;
  const dbRows = db === 'All' ? databaseRows : databaseRows.filter(row => row.database === db);
  const featureCount = active.length;
  const rowCount = active.reduce((a, r) => a + num(r.feature_rows), 0);
  const positive = Math.max(...active.map(r => num(r.positive_genomes)), 0);
  document.getElementById('summary').innerHTML = `<div class="cards">
    <div class="card"><span>Displayed features</span><strong>${featureCount}</strong></div>
    <div class="card"><span>Displayed feature rows</span><strong>${rowCount}</strong></div>
    <div class="card"><span>Max positive genomes</span><strong>${positive}</strong></div>
    <div class="card"><span>Databases represented</span><strong>${dbRows.length || new Set(active.map(r => r.database)).size}</strong></div>
  </div>`;
  const written = writtenRows.find(row => row.database === db) || writtenRows.find(row => row.scope === 'overall') || {};
  document.getElementById('written').textContent = written.summary || '';
}
function renderFigure(active) {
  const metric = document.getElementById('metric').value;
  const figureRows = active.slice(0, 50);
  const width = Math.max(900, figureRows.length * 42 + 180), height = 450;
  const left = 110, top = 42, plotHeight = 270, plotWidth = width - left - 42;
  const maxValue = Math.max(1, ...figureRows.map(row => num(row[metric])));
  let svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`;
  svg += `<rect width="100%" height="100%" fill="#f8fafc"/><text x="20" y="28" font-size="20" font-family="Arial" font-weight="700" fill="#102a43">Feature prevalence</text>`;
  svg += `<rect x="${left}" y="${top}" width="${plotWidth}" height="${plotHeight}" fill="#fff" stroke="#bcccdc"/>`;
  figureRows.forEach((row, idx) => {
    const value = num(row[metric]), x = left + 14 + idx * 42;
    const barHeight = value / maxValue * (plotHeight - 18), y = top + plotHeight - barHeight;
    svg += `<rect x="${x}" y="${y}" width="24" height="${barHeight}" fill="#0f766e"><title>${row.database}:${row.feature_id} ${row.prevalence_display}</title></rect>`;
    svg += `<text x="${x + 8}" y="${top + plotHeight + 14}" transform="rotate(70 ${x + 8} ${top + plotHeight + 14})" font-size="10" font-family="Arial" fill="#334e68">${(row.feature_id || '').slice(0, 24)}</text>`;
  });
  if (figureRows.length < active.length) svg += `<text x="20" y="${height - 18}" font-size="12" font-family="Arial" fill="#52606d">Figure capped at 50 features. Complete table and TSV preserve all rows.</text>`;
  svg += '</svg>';
  document.getElementById('figure').innerHTML = svg;
}
function renderTable(active) {
  const fields = ['database','feature_id','feature_category','positive_genomes','total_genomes','prevalence_display','feature_rows','mean_hits_per_positive_genome','prevalence_label','warning_flags'];
  const shown = active.slice(0, document.getElementById('view').value === 'all' ? active.length : 80);
  let table = '<table><thead><tr>' + fields.map(field => `<th>${field}</th>`).join('') + '</tr></thead><tbody>';
  for (const row of shown) table += '<tr>' + fields.map(field => `<td>${row[field] || ''}</td>`).join('') + '</tr>';
  table += '</tbody></table>';
  if (shown.length < active.length) table += `<p>Showing ${shown.length} of ${active.length} rows. Use Complete or download the full TSV.</p>`;
  document.getElementById('table').innerHTML = table;
}
function render() { const active = activeRows(); renderCards(active); renderFigure(active); renderTable(active); }
for (const id of ['database','view','metric','sort','minPrev','minPositive','coordsOnly','categoryOnly']) document.getElementById(id).addEventListener('change', render);
render();
</script></body></html>
"""
    html_path.write_text(
        prevalence_html
        .replace("__FEATURE_ROWS__", json.dumps(feature_rows))
        .replace("__DATABASE_ROWS__", json.dumps(summary_by_database))
        .replace("__WRITTEN_ROWS__", json.dumps(written_rows)),
        encoding="utf-8",
    )
    outputs["important_prevalence_analysis_html"] = str(html_path)
    table_zip = important_dir / "prevalence_tables.zip"
    figure_zip = important_dir / "prevalence_figures.zip"
    outputs["important_prevalence_tables_zip"] = _write_zip_bundle(
        table_zip,
        [feature_prevalence_path, top_path, summary_db_path, core_path, burden_path, written_path],
        important_dir,
    )
    outputs["important_prevalence_figures_zip"] = _write_zip_bundle(figure_zip, figure_files + [html_path], important_dir)
    return outputs


def _summary_stats_full(values: list[float]) -> dict[str, float | None]:
    if not values:
        return {"median": None, "min": None, "max": None, "q1": None, "q3": None, "iqr": None, "mean": None}
    values = sorted(values)
    mid = len(values) // 2
    median = values[mid] if len(values) % 2 else (values[mid - 1] + values[mid]) / 2
    q1 = values[len(values) // 4]
    q3 = values[(len(values) * 3) // 4]
    return {
        "median": median,
        "min": values[0],
        "max": values[-1],
        "q1": q1,
        "q3": q3,
        "iqr": q3 - q1,
        "mean": sum(values) / len(values),
    }


def _summary_stats(values: list[float]) -> dict[str, str]:
    stats = _summary_stats_full(values)
    return {
        "median": "" if stats["median"] is None else f"{stats['median']:.2f}",
        "min": "" if stats["min"] is None else f"{stats['min']:.2f}",
        "max": "" if stats["max"] is None else f"{stats['max']:.2f}",
        "iqr": "" if stats["iqr"] is None else f"{stats['iqr']:.2f}",
    }


def _alignment_length(row: dict[str, str]) -> str:
    direct = first_value(row, ["alignment_length", "alignment_len", "align_len", "length"], "")
    if direct:
        return direct
    start = _float_or_none(row.get("start", ""))
    end = _float_or_none(row.get("end", ""))
    if start is None or end is None:
        return ""
    length = abs(end - start) + 1
    return str(int(length)) if length > 0 else ""


def _variation_boxplot_rows(
    hit_rows: list[dict[str, str]],
    selected_features: set[tuple[str, str]],
    metric_field: str,
    metric_label: str,
) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str], list[float]] = defaultdict(list)
    for row in hit_rows:
        key = (row.get("database", ""), row.get("feature_id", ""))
        if key not in selected_features:
            continue
        value = _float_or_none(row.get(metric_field, ""))
        if value is not None:
            grouped[key].append(value)
    rows = []
    for (database, feature_id), values in sorted(grouped.items()):
        stats = _summary_stats_full(values)
        rows.append({
            "database": database,
            "feature_id": feature_id,
            "metadata_group": feature_id,
            "metric": metric_label,
            "n": str(len(values)),
            "min": "" if stats["min"] is None else f"{stats['min']:.2f}",
            "q1": "" if stats["q1"] is None else f"{stats['q1']:.2f}",
            "median": "" if stats["median"] is None else f"{stats['median']:.2f}",
            "q3": "" if stats["q3"] is None else f"{stats['q3']:.2f}",
            "max": "" if stats["max"] is None else f"{stats['max']:.2f}",
            "mean": "" if stats["mean"] is None else f"{stats['mean']:.2f}",
        })
    return rows


def _write_variation_scatter_svg(path: Path, rows: list[dict[str, str]], title: str) -> None:
    width, height = 900, 560
    left, top, plot_w, plot_h = 78, 62, 720, 410
    points = [
        (
            _float_or_none(row.get("identity", "")),
            _float_or_none(row.get("coverage", "")),
            row,
        )
        for row in rows
    ]
    points = [(x, y, row) for x, y, row in points if x is not None and y is not None]
    parts = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect width='100%' height='100%' fill='#f8fafc'/>",
        f"<text x='20' y='32' font-family='Arial' font-size='20' font-weight='700' fill='#102a43'>{html.escape(title)}</text>",
        f"<rect x='{left}' y='{top}' width='{plot_w}' height='{plot_h}' fill='#ffffff' stroke='#bcccdc'/>",
    ]
    for tick in range(0, 101, 20):
        x = left + tick / 100 * plot_w
        y = top + plot_h - tick / 100 * plot_h
        parts.append(f"<line x1='{x:.1f}' y1='{top}' x2='{x:.1f}' y2='{top + plot_h}' stroke='#e2e8f0'/>")
        parts.append(f"<line x1='{left}' y1='{y:.1f}' x2='{left + plot_w}' y2='{y:.1f}' stroke='#e2e8f0'/>")
        parts.append(f"<text x='{x - 8:.1f}' y='{top + plot_h + 22}' font-size='11' font-family='Arial' fill='#52606d'>{tick}</text>")
        parts.append(f"<text x='{left - 38}' y='{y + 4:.1f}' font-size='11' font-family='Arial' fill='#52606d'>{tick}</text>")
    colors = ["#0f766e", "#dc2626", "#7c3aed", "#d97706", "#2563eb", "#be123c", "#15803d"]
    feature_colors: dict[str, str] = {}
    for x_value, y_value, row in points[:600]:
        feature = row.get("feature_id", "")
        if feature not in feature_colors:
            feature_colors[feature] = colors[len(feature_colors) % len(colors)]
        cx = left + max(0, min(100, x_value)) / 100 * plot_w
        cy = top + plot_h - max(0, min(100, y_value)) / 100 * plot_h
        tip = f"{row.get('database', '')} {feature}; identity={x_value:.2f}; coverage={y_value:.2f}; sample={row.get('assembly_accession', '') or row.get('sample_id', '')}"
        parts.append(f"<circle cx='{cx:.1f}' cy='{cy:.1f}' r='4' fill='{feature_colors[feature]}' fill-opacity='0.72'><title>{html.escape(tip)}</title></circle>")
    parts.append(f"<text x='{left + plot_w / 2 - 55}' y='{height - 24}' font-size='13' font-family='Arial' fill='#334e68'>Identity (%)</text>")
    parts.append(f"<text x='20' y='{top + plot_h / 2}' font-size='13' font-family='Arial' fill='#334e68' transform='rotate(-90 20 {top + plot_h / 2})'>Coverage (%)</text>")
    parts.append("</svg>\n")
    path.write_text("".join(parts), encoding="utf-8")


def _write_variation_scatter_png(path: Path, rows: list[dict[str, str]]) -> None:
    width, height = 900, 560
    left, top, plot_w, plot_h = 78, 62, 720, 410
    pixels = [bytearray([248, 250, 252] * width) for _ in range(height)]

    def rect(x0: int, y0: int, x1: int, y1: int, color: tuple[int, int, int]) -> None:
        for y in range(max(0, y0), min(height, y1)):
            row_pixels = pixels[y]
            for x in range(max(0, x0), min(width, x1)):
                idx = x * 3
                row_pixels[idx:idx + 3] = bytes(color)

    def dot(cx: int, cy: int, color: tuple[int, int, int]) -> None:
        for y in range(cy - 3, cy + 4):
            for x in range(cx - 3, cx + 4):
                if (x - cx) * (x - cx) + (y - cy) * (y - cy) <= 9:
                    if 0 <= x < width and 0 <= y < height:
                        idx = x * 3
                        pixels[y][idx:idx + 3] = bytes(color)

    rect(left, top, left + plot_w, top + plot_h, (255, 255, 255))
    colors = [(15, 118, 110), (220, 38, 38), (124, 58, 237), (217, 119, 6), (37, 99, 235), (190, 18, 60), (21, 128, 61)]
    feature_colors: dict[str, tuple[int, int, int]] = {}
    for row in rows[:600]:
        identity = _float_or_none(row.get("identity", ""))
        coverage = _float_or_none(row.get("coverage", ""))
        if identity is None or coverage is None:
            continue
        feature = row.get("feature_id", "")
        if feature not in feature_colors:
            feature_colors[feature] = colors[len(feature_colors) % len(colors)]
        x = int(left + max(0, min(100, identity)) / 100 * plot_w)
        y = int(top + plot_h - max(0, min(100, coverage)) / 100 * plot_h)
        dot(x, y, feature_colors[feature])
    _write_png(path, width, height, pixels)


def write_important_variation_outputs(sample_dir: Path, out_dir: Path, important_dir: Path, top_n: int = 20) -> dict[str, str]:
    key_tables = important_dir / "key_tables"
    figures = important_dir / "figures"
    key_tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)
    metadata_rows = read_table(sample_dir / "basic" / "enriched_genome_dataset.csv")
    sample_count = len(metadata_rows)
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
            "alignment_length": _alignment_length(row),
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
        alignment_lengths = [_float_or_none(_alignment_length(row)) for row in rows]
        identities = [value for value in identities if value is not None]
        coverages = [value for value in coverages if value is not None]
        alignment_lengths = [value for value in alignment_lengths if value is not None]
        identity_stats = _summary_stats(identities)
        coverage_stats = _summary_stats(coverages)
        alignment_stats = _summary_stats(alignment_lengths)
        samples = {row.get("assembly_accession", "") or row.get("sample_id", "") for row in rows if row.get("assembly_accession", "") or row.get("sample_id", "")}
        low_identity = sum(1 for value in identities if value < 90)
        low_coverage = sum(1 for value in coverages if value < 80)
        warnings = []
        if low_identity:
            warnings.append("low_identity")
        if low_coverage:
            warnings.append("low_coverage")
            warnings.append("partial_feature")
        if len(rows) < 3:
            warnings.append("few_hits")
        iqr_identity = _float_or_none(identity_stats["iqr"]) or 0.0
        iqr_coverage = _float_or_none(coverage_stats["iqr"]) or 0.0
        variation_score = max(iqr_identity, iqr_coverage)
        label = "high_variation" if variation_score >= 10 else ("moderate_variation" if variation_score >= 3 else "low_variation")
        if label == "high_variation":
            warnings.append("high_variation")
        prevalence = len(samples) / sample_count if sample_count else 0.0
        summary_rows.append({
            "database": database,
            "feature_id": feature_id,
            "feature_name": first_value(rows[0], ["feature_name", "feature_id"], feature_id),
            "total_hits": str(len(rows)),
            "positive_genomes": str(len(samples)),
            "sample_count": str(sample_count),
            "prevalence_percent": f"{prevalence * 100:.1f}",
            "mean_hits_per_positive_genome": f"{(len(rows) / len(samples)):.2f}" if samples else "",
            "median_identity": identity_stats["median"],
            "min_identity": identity_stats["min"],
            "max_identity": identity_stats["max"],
            "iqr_identity": identity_stats["iqr"],
            "median_coverage": coverage_stats["median"],
            "min_coverage": coverage_stats["min"],
            "max_coverage": coverage_stats["max"],
            "iqr_coverage": coverage_stats["iqr"],
            "median_alignment_length": alignment_stats["median"],
            "min_alignment_length": alignment_stats["min"],
            "max_alignment_length": alignment_stats["max"],
            "iqr_alignment_length": alignment_stats["iqr"],
            "low_identity_hits": str(low_identity),
            "low_coverage_hits": str(low_coverage),
            "variation_score": f"{variation_score:.2f}",
            "variation_label": label,
            "warning_flags": ";".join(dict.fromkeys(warnings)),
        })
    summary_rows.sort(key=lambda row: (row["database"], -(_float_or_none(row["iqr_identity"]) or 0.0), row["feature_id"]))
    summary_fields = [
        "database", "feature_id", "feature_name", "total_hits", "positive_genomes", "sample_count", "prevalence_percent", "mean_hits_per_positive_genome",
        "median_identity", "min_identity", "max_identity", "iqr_identity",
        "median_coverage", "min_coverage", "max_coverage", "iqr_coverage",
        "median_alignment_length", "min_alignment_length", "max_alignment_length", "iqr_alignment_length",
        "low_identity_hits", "low_coverage_hits", "variation_score", "variation_label", "warning_flags",
    ]
    summary_path = key_tables / "feature_variation_summary.tsv"
    hits_path = key_tables / "feature_variation_hits.tsv"
    write_rows(summary_path, summary_rows, summary_fields)
    hit_fields = ["database", "feature_id", "assembly_accession", "sample_id", "identity", "coverage", "alignment_length", "contig", "start", "end", "tool", "source_file"]
    write_rows(hits_path, hit_rows, hit_fields)

    outputs = {
        "important_feature_variation_summary": str(summary_path),
        "important_feature_variation_hits": str(hits_path),
    }
    figure_files: list[Path] = []

    def write_pdf_for_rows(path: Path, title: str, rows: list[dict[str, str]], metric: str) -> None:
        _write_simple_pdf(
            path,
            title,
            [f"{row.get('feature_id', row.get('metadata_group', ''))}: {metric}={row.get(metric, row.get('median', ''))}, n={row.get('n', row.get('total_hits', ''))}" for row in rows],
        )

    for database in sorted({row["database"] for row in summary_rows}):
        safe_database = _safe_filename(database)
        db_summary = [row for row in summary_rows if row["database"] == database]
        db_rows = sorted(db_summary, key=lambda row: -(_float_or_none(row["iqr_identity"]) or 0.0))[:top_n]
        if not db_rows:
            continue
        for metric_name, metric_field, iqr_field, title_metric in [
            ("identity", "identity", "iqr_identity", "Identity"),
            ("coverage", "coverage", "iqr_coverage", "Coverage"),
        ]:
            metric_rows = sorted(db_summary, key=lambda row: -(_float_or_none(row.get(iqr_field, "")) or 0.0))[:top_n]
            selected = {(row["database"], row["feature_id"]) for row in metric_rows}
            box_rows = _variation_boxplot_rows(hit_rows, selected, metric_field, title_metric)
            data_path = figures / f"variation_{metric_name}_{safe_database}_top20.data.tsv"
            svg_path = figures / f"variation_{metric_name}_{safe_database}_top20.svg"
            png_path = figures / f"variation_{metric_name}_{safe_database}_top20.png"
            pdf_path = figures / f"variation_{metric_name}_{safe_database}_top20.pdf"
            write_rows(data_path, box_rows, ["database", "feature_id", "metadata_group", "metric", "n", "min", "q1", "median", "q3", "max", "mean"])
            _write_burden_boxplot_svg(svg_path, box_rows, f"{database} {title_metric} Variation Top {top_n}")
            _write_burden_boxplot_png(png_path, box_rows)
            _write_simple_pdf(pdf_path, f"{database} {title_metric} Variation Top {top_n}", [f"{row.get('feature_id', '')}: median={row.get('median', '')}, IQR={row.get('q1', '')}-{row.get('q3', '')}, n={row.get('n', '')}" for row in box_rows])
            figure_files.extend([data_path, svg_path, png_path, pdf_path])
            outputs[f"important_variation_{safe_database}_{metric_name}_data"] = str(data_path)
            outputs[f"important_variation_{safe_database}_{metric_name}_svg"] = str(svg_path)
            outputs[f"important_variation_{safe_database}_{metric_name}_png"] = str(png_path)
            outputs[f"important_variation_{safe_database}_{metric_name}_pdf"] = str(pdf_path)

        scatter_features = {(row["database"], row["feature_id"]) for row in sorted(db_summary, key=lambda row: -(_float_or_none(row.get("variation_score", "")) or 0.0))[:top_n]}
        scatter_rows = [row for row in hit_rows if (row.get("database", ""), row.get("feature_id", "")) in scatter_features]
        scatter_base = f"variation_identity_coverage_{safe_database}_top20"
        scatter_data = figures / f"{scatter_base}.data.tsv"
        scatter_svg = figures / f"{scatter_base}.svg"
        scatter_png = figures / f"{scatter_base}.png"
        scatter_pdf = figures / f"{scatter_base}.pdf"
        write_rows(scatter_data, scatter_rows, hit_fields)
        _write_variation_scatter_svg(scatter_svg, scatter_rows, f"{database} Identity vs Coverage")
        _write_variation_scatter_png(scatter_png, scatter_rows)
        _write_simple_pdf(scatter_pdf, f"{database} Identity vs Coverage", [f"{row.get('feature_id', '')}: identity={row.get('identity', '')}, coverage={row.get('coverage', '')}" for row in scatter_rows[:34]])
        figure_files.extend([scatter_data, scatter_svg, scatter_png, scatter_pdf])

        variable_rows = sorted(db_summary, key=lambda row: -(_float_or_none(row.get("variation_score", "")) or 0.0))[:top_n]
        variable_base = f"variation_top_variable_{safe_database}_top20"
        variable_data = figures / f"{variable_base}.data.tsv"
        variable_svg = figures / f"{variable_base}.svg"
        variable_png = figures / f"{variable_base}.png"
        variable_pdf = figures / f"{variable_base}.pdf"
        write_rows(variable_data, variable_rows, summary_fields)
        _write_bar_svg(variable_svg, variable_rows, f"{database} Top Variable Features", "feature_id", "variation_score", "Variation score")
        _write_bar_png(variable_png, variable_rows, "variation_score")
        write_pdf_for_rows(variable_pdf, f"{database} Top Variable Features", variable_rows, "variation_score")
        figure_files.extend([variable_data, variable_svg, variable_png, variable_pdf])
        outputs[f"important_variation_{safe_database}_scatter_svg"] = str(scatter_svg)
        outputs[f"important_variation_{safe_database}_scatter_png"] = str(scatter_png)
        outputs[f"important_variation_{safe_database}_scatter_pdf"] = str(scatter_pdf)
        outputs[f"important_variation_{safe_database}_top_variable_svg"] = str(variable_svg)
        outputs[f"important_variation_{safe_database}_top_variable_png"] = str(variable_png)
        outputs[f"important_variation_{safe_database}_top_variable_pdf"] = str(variable_pdf)

    summary_by_database = []
    for database in sorted({row["database"] for row in summary_rows}):
        db_rows = [row for row in summary_rows if row["database"] == database]
        total_hits = sum(int(_float_or_none(row.get("total_hits", "")) or 0) for row in db_rows)
        median_identity_values = [_float_or_none(row.get("median_identity", "")) for row in db_rows]
        median_coverage_values = [_float_or_none(row.get("median_coverage", "")) for row in db_rows]
        median_identity_values = [value for value in median_identity_values if value is not None]
        median_coverage_values = [value for value in median_coverage_values if value is not None]
        identity_summary = _summary_stats(median_identity_values)
        coverage_summary = _summary_stats(median_coverage_values)
        summary_by_database.append({
            "database": database,
            "unique_features": str(len(db_rows)),
            "total_hits": str(total_hits),
            "median_identity": identity_summary["median"],
            "median_coverage": coverage_summary["median"],
            "high_variation_features": str(sum(1 for row in db_rows if row.get("variation_label") == "high_variation")),
            "low_identity_hits": str(sum(int(_float_or_none(row.get("low_identity_hits", "")) or 0) for row in db_rows)),
            "low_coverage_hits": str(sum(int(_float_or_none(row.get("low_coverage_hits", "")) or 0) for row in db_rows)),
        })
    variation_summary_path = key_tables / "feature_variation_database_summary.tsv"
    write_rows(
        variation_summary_path,
        summary_by_database,
        ["database", "unique_features", "total_hits", "median_identity", "median_coverage", "high_variation_features", "low_identity_hits", "low_coverage_hits"],
    )
    outputs["important_feature_variation_database_summary"] = str(variation_summary_path)

    html_path = figures / "variation_analysis.html"
    variation_html = """<!doctype html>
<html><head><meta charset="utf-8"><title>Feature Variations</title>
<style>
body { font-family: Arial, sans-serif; margin: 1.5rem; color: #1f2933; background: #f8fafc; }
label { font-weight: 700; margin-right: 0.35rem; }
select { margin: 0 0.9rem 1rem 0; padding: 0.35rem; }
.warning { background: #fff7ed; border-left: 4px solid #c2410c; padding: 0.75rem; margin: 1rem 0; }
.cards { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 0.75rem; }
.card { border: 1px solid #d9e2ec; border-radius: 6px; padding: 0.75rem; background: white; }
.card span { display: block; color: #52606d; font-size: 0.8rem; }
.card strong { display: block; font-size: 1.25rem; margin-top: 0.25rem; }
.figure-box { max-width: 100%; max-height: 680px; overflow: auto; border: 1px solid #d9e2ec; background: white; padding: 0.75rem; }
svg { max-width: 100%; height: auto; }
table { border-collapse: collapse; width: 100%; margin-top: 1rem; background: white; font-size: 0.9rem; }
th, td { border: 1px solid #d9e2ec; padding: 0.4rem; text-align: left; }
th { background: #f0f4f8; position: sticky; top: 0; }
.downloads a { display: inline-block; margin: 0.25rem 0.45rem 0.25rem 0; color: #0f766e; font-weight: 700; }
</style></head><body>
<h1>Variations</h1>
<p>Use this view to compare identity, coverage, alignment length, and hit-count variation across detected features. Complete views are table-first so large feature sets do not force the full page to scroll sideways.</p>
<div class="warning">Low identity, low coverage, high variation, and few-hit labels are review flags. They may reflect divergent homologs, partial hits, fragmented assemblies, or database boundary effects.</div>
<label for="database">Database</label><select id="database"></select>
<label for="metric">Metric</label><select id="metric">
  <option value="identity">Identity</option>
  <option value="coverage">Coverage</option>
  <option value="alignment_length">Alignment length</option>
  <option value="hit_count">Feature count per genome</option>
</select>
<label for="view">View</label><select id="view">
  <option value="variable">Top variable features</option>
  <option value="prevalent">Top prevalent features</option>
  <option value="complete">Complete</option>
</select>
<label for="display">Display</label><select id="display">
  <option value="20">Top 20</option>
  <option value="10">Top 10</option>
  <option value="50">Top 50</option>
  <option value="all">Complete</option>
</select>
<label for="sort">Sort by</label><select id="sort">
  <option value="variation">variation</option>
  <option value="prevalence">prevalence</option>
  <option value="median_identity">median identity</option>
  <option value="low_identity">low-identity warnings</option>
  <option value="alphabetical">alphabetical</option>
</select>
<div id="summary"></div>
<div class="figure-box" id="figure"></div>
<div id="table"></div>
<div class="downloads">
  <a href="../key_tables/feature_variation_summary.tsv">Download full variation summary</a>
  <a href="../key_tables/feature_variation_hits.tsv">Download hit-level variation table</a>
  <a href="../key_tables/feature_variation_database_summary.tsv">Download database summary</a>
</div>
<script>
const summaryRows = __SUMMARY_ROWS__;
const databaseRows = __DATABASE_ROWS__;
function num(value) { const parsed = Number(value || 0); return Number.isFinite(parsed) ? parsed : 0; }
function metricField(metric) {
  if (metric === 'identity') return 'iqr_identity';
  if (metric === 'coverage') return 'iqr_coverage';
  if (metric === 'alignment_length') return 'iqr_alignment_length';
  return 'mean_hits_per_positive_genome';
}
function metricLabel(metric) {
  if (metric === 'identity') return 'Identity IQR';
  if (metric === 'coverage') return 'Coverage IQR';
  if (metric === 'alignment_length') return 'Alignment length IQR';
  return 'Mean hits per positive genome';
}
const dbSelect = document.getElementById('database');
for (const db of ['All', ...[...new Set(summaryRows.map(r => r.database).filter(Boolean))].sort()]) {
  const option = document.createElement('option'); option.value = db; option.textContent = db; dbSelect.appendChild(option);
}
function activeRows() {
  const database = dbSelect.value;
  const metric = document.getElementById('metric').value;
  const view = document.getElementById('view').value;
  const sort = document.getElementById('sort').value;
  const display = document.getElementById('display').value;
  let rows = summaryRows.filter(row => database === 'All' || row.database === database);
  if (view === 'prevalent') rows = rows.slice().sort((a, b) => num(b.prevalence_percent) - num(a.prevalence_percent));
  else if (sort === 'prevalence') rows = rows.slice().sort((a, b) => num(b.prevalence_percent) - num(a.prevalence_percent));
  else if (sort === 'median_identity') rows = rows.slice().sort((a, b) => num(a.median_identity) - num(b.median_identity));
  else if (sort === 'low_identity') rows = rows.slice().sort((a, b) => num(b.low_identity_hits) - num(a.low_identity_hits));
  else if (sort === 'alphabetical') rows = rows.slice().sort((a, b) => (a.feature_id || '').localeCompare(b.feature_id || ''));
  else rows = rows.slice().sort((a, b) => num(b[metricField(metric)]) - num(a[metricField(metric)]));
  if (view === 'complete' || display === 'all') return rows;
  return rows.slice(0, Number(display || 20));
}
function renderCards(rows) {
  const unique = rows.length;
  const hits = rows.reduce((a, r) => a + num(r.total_hits), 0);
  const high = rows.filter(r => r.variation_label === 'high_variation').length;
  const lowIdentity = rows.reduce((a, r) => a + num(r.low_identity_hits), 0);
  const lowCoverage = rows.reduce((a, r) => a + num(r.low_coverage_hits), 0);
  const medianIdentityValues = rows.map(r => num(r.median_identity)).filter(v => v > 0).sort((a, b) => a - b);
  const medianCoverageValues = rows.map(r => num(r.median_coverage)).filter(v => v > 0).sort((a, b) => a - b);
  const mid = values => values.length ? values[Math.floor(values.length / 2)].toFixed(1) : '';
  document.getElementById('summary').innerHTML = `<div class="cards">
    <div class="card"><span>Unique features analyzed</span><strong>${unique}</strong></div>
    <div class="card"><span>Total hits</span><strong>${hits}</strong></div>
    <div class="card"><span>Median identity</span><strong>${mid(medianIdentityValues)}</strong></div>
    <div class="card"><span>Median coverage</span><strong>${mid(medianCoverageValues)}</strong></div>
    <div class="card"><span>High-variation features</span><strong>${high}</strong></div>
    <div class="card"><span>Low identity / coverage hits</span><strong>${lowIdentity} / ${lowCoverage}</strong></div>
  </div>`;
}
function renderFigure(rows) {
  const metric = document.getElementById('metric').value;
  const field = metricField(metric);
  const figureRows = rows.slice(0, 50);
  const width = Math.max(860, figureRows.length * 46 + 170), height = 455;
  const plotLeft = 120, plotTop = 42, plotHeight = 270, plotWidth = width - plotLeft - 40;
  const maxValue = Math.max(1, ...figureRows.map(row => num(row[field])));
  let svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`;
  svg += `<rect width="100%" height="100%" fill="#f8fafc"/><text x="20" y="28" font-size="20" font-family="Arial" font-weight="700" fill="#102a43">${metricLabel(metric)}</text>`;
  svg += `<rect x="${plotLeft}" y="${plotTop}" width="${plotWidth}" height="${plotHeight}" fill="#fff" stroke="#bcccdc"/>`;
  figureRows.forEach((row, idx) => {
    const value = num(row[field]);
    const x = plotLeft + 18 + idx * 46;
    const barHeight = value / maxValue * (plotHeight - 20);
    const y = plotTop + plotHeight - barHeight;
    const color = row.warning_flags ? '#dc2626' : '#0f766e';
    svg += `<rect x="${x}" y="${y}" width="26" height="${barHeight}" fill="${color}"><title>${row.database} ${row.feature_id}: ${value}</title></rect>`;
    svg += `<text x="${x + 8}" y="${plotTop + plotHeight + 14}" transform="rotate(70 ${x + 8} ${plotTop + plotHeight + 14})" font-size="10" font-family="Arial" fill="#334e68">${(row.feature_id || '').slice(0, 24)}</text>`;
  });
  if (figureRows.length < rows.length) svg += `<text x="20" y="${height - 18}" font-size="12" font-family="Arial" fill="#52606d">Figure capped at 50 features. Use the table/download for the complete view.</text>`;
  svg += `</svg>`;
  document.getElementById('figure').innerHTML = svg;
}
function renderTable(rows) {
  const fields = ['database','feature_id','total_hits','positive_genomes','prevalence_percent','mean_hits_per_positive_genome','median_identity','iqr_identity','median_coverage','iqr_coverage','median_alignment_length','iqr_alignment_length','variation_label','warning_flags'];
  const shown = rows.slice(0, document.getElementById('display').value === 'all' ? rows.length : 80);
  let table = '<table><thead><tr>' + fields.map(field => `<th>${field}</th>`).join('') + '</tr></thead><tbody>';
  for (const row of shown) table += '<tr>' + fields.map(field => `<td>${row[field] || ''}</td>`).join('') + '</tr>';
  table += '</tbody></table>';
  if (shown.length < rows.length) table += `<p>Showing ${shown.length} of ${rows.length} rows. Use Complete to inspect all rows, or download the full TSV.</p>`;
  document.getElementById('table').innerHTML = table;
}
function render() {
  const rows = activeRows();
  renderCards(rows);
  renderFigure(rows);
  renderTable(rows);
}
for (const id of ['database','metric','view','display','sort']) document.getElementById(id).addEventListener('change', render);
render();
</script></body></html>
"""
    html_path.write_text(
        variation_html
        .replace("__SUMMARY_ROWS__", json.dumps(summary_rows))
        .replace("__DATABASE_ROWS__", json.dumps(summary_by_database)),
        encoding="utf-8",
    )
    outputs["important_variation_analysis_html"] = str(html_path)
    outputs["important_variation_figures_zip"] = _write_zip_bundle(important_dir / "variation_figures.zip", figure_files, important_dir)
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


def write_important_metadata_association_outputs(
    sample_dir: Path,
    out_dir: Path,
    important_dir: Path,
    top_n: int = 20,
    max_features_per_database: int = 200,
) -> dict[str, str]:
    tables = important_dir / "tables"
    figures = important_dir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    features = read_table(out_dir / "features" / "all_features.tsv")
    metadata_rows = read_table(sample_dir / "basic" / "enriched_genome_dataset.csv")
    if not metadata_rows:
        metadata_rows = normalize_metadata_rows(load_metadata_rows(sample_dir))
    metadata_by_sample = {row.get("assembly_accession", ""): row for row in metadata_rows if row.get("assembly_accession")}
    samples = sorted(set(metadata_by_sample) | {row.get("assembly_accession", "") for row in features if row.get("assembly_accession")})
    sample_count = len(samples)
    presence = feature_presence(features)
    category_by_feature: dict[tuple[str, str], str] = {}
    for row in features:
        key = (row.get("database", ""), row.get("feature_id", ""))
        if key[0] and key[1]:
            category_by_feature[key] = row.get("feature_category", "")

    by_database: dict[str, list[tuple[tuple[str, str], set[str]]]] = defaultdict(list)
    for key, present_samples in presence.items():
        by_database[key[0]].append((key, present_samples))
    limited_presence: dict[tuple[str, str], set[str]] = {}
    capped_databases = []
    for database, items in by_database.items():
        ranked = sorted(items, key=lambda item: (-len(item[1]), item[0][1]))[:max_features_per_database]
        if len(items) > max_features_per_database:
            capped_databases.append(database)
        for key, present_samples in ranked:
            limited_presence[key] = present_samples

    preferred_columns = [
        "isolation_source",
        "country",
        "continent",
        "subcontinent",
        "region",
        "host",
        "host_group",
        "sample_type",
        "environment_medium",
        "collection_year",
        "bioproject",
        "biosample",
        "mlst_ST",
        "ani_cluster",
        "mash_cluster",
        "dominant_lineage_label",
        "assembly_level",
    ]
    excluded_metadata_columns = {
        "assembly_accession",
        "sample_id",
        "organism_name",
        "taxid",
        "ftp_path_refseq",
        "ftp_path_genbank",
        "qc_status",
        "qc_pass",
        "qc_fail_reasons",
        "quast_status",
        "ani_status",
        "ani_species_match",
        "mash_status",
        "genome_size",
        "contig_count",
        "n50",
        "gc_percent",
        "checkm2_completeness",
        "checkm2_contamination",
        "amr_genes",
        "amrfinderplus_genes",
        "drug_classes",
        "resistance_mechanisms",
        "vfdb_genes",
        "vfdb_categories",
        "plasmid_replicons",
        "integron_features",
        "mobsuite_plasmid_types",
        "genomad_regions",
        "defense_systems",
        "mobile_elements",
        "features_detected_databases",
        "modules_run",
        "modules_failed",
        "modules_warning",
        "panresistome_version",
        "feature_contract_version",
    }

    def metadata_candidate_columns() -> list[str]:
        all_columns = sorted({column for row in metadata_rows for column in row})
        candidates = []
        for column in [*preferred_columns, *all_columns]:
            if column in candidates or column in excluded_metadata_columns:
                continue
            if column.endswith("_count") or column.endswith("_feature_count") or column.endswith("_gene_count"):
                continue
            candidates.append(column)
        return candidates

    usability_rows = []
    for column in metadata_candidate_columns():
        values = [row.get(column, "") for row in metadata_rows]
        non_missing_values = [value for value in values if value and not is_missing_value(value)]
        counts = Counter(non_missing_values)
        largest_group, largest_group_count = ("", 0)
        if counts:
            largest_group, largest_group_count = max(counts.items(), key=lambda item: (item[1], item[0]))
        non_missing_count = len(non_missing_values)
        missing_count = max(len(values) - non_missing_count, 0)
        missing_fraction = (missing_count / len(values)) if values else 1.0
        unique_values = len(counts)
        largest_group_fraction = (largest_group_count / non_missing_count) if non_missing_count else 0.0
        eligible = unique_values >= 2 and non_missing_count >= 2
        reasons = []
        if not non_missing_count:
            reasons.append("all_missing")
        if unique_values < 2:
            reasons.append("fewer_than_two_groups")
        if non_missing_count < 10:
            reasons.append("small_non_missing_count")
        if missing_fraction >= 0.50:
            reasons.append("high_missingness")
        if largest_group_fraction >= 0.80 and unique_values >= 2:
            reasons.append("dominant_group")
        if not reasons:
            reasons.append("usable")
        if not eligible:
            recommended_use = "exclude"
        elif non_missing_count < 10 or missing_fraction >= 0.50 or largest_group_fraction >= 0.80:
            recommended_use = "descriptive_or_bias_check"
        else:
            recommended_use = "association_testing"
        usability_rows.append({
            "metadata_column": column,
            "non_missing_count": str(non_missing_count),
            "missing_count": str(missing_count),
            "missing_fraction": f"{missing_fraction:.4f}",
            "unique_values": str(unique_values),
            "largest_group": largest_group,
            "largest_group_count": str(largest_group_count),
            "largest_group_fraction": f"{largest_group_fraction:.4f}",
            "eligible_for_testing": "true" if eligible else "false",
            "reason": ";".join(reasons),
            "recommended_use": recommended_use,
        })
    usability_path = tables / "metadata_usability_summary.tsv"
    usability_fields = [
        "metadata_column",
        "non_missing_count",
        "missing_count",
        "missing_fraction",
        "unique_values",
        "largest_group",
        "largest_group_count",
        "largest_group_fraction",
        "eligible_for_testing",
        "reason",
        "recommended_use",
    ]
    write_rows(usability_path, usability_rows, usability_fields)

    metadata_columns = []
    missing_fraction_by_column = {}
    for column in metadata_candidate_columns():
        values = [row.get(column, "") for row in metadata_rows]
        non_missing = [value for value in values if value and not is_missing_value(value)]
        if len(set(non_missing)) >= 2:
            metadata_columns.append(column)
            missing_fraction_by_column[column] = 1.0 - (len(non_missing) / len(values) if values else 0.0)

    def groups_for_column(column: str) -> dict[str, set[str]]:
        groups: dict[str, set[str]] = defaultdict(set)
        for sample in samples:
            value = metadata_by_sample.get(sample, {}).get(column, "")
            if value and not is_missing_value(value):
                groups[value].add(sample)
        return dict(groups)

    feature_rows = []
    for metadata_column in metadata_columns:
        groups = groups_for_column(metadata_column)
        if len(groups) < 2:
            continue
        all_column_samples = set().union(*groups.values()) if groups else set()
        for (database, feature_id), present_samples in sorted(limited_presence.items()):
            for group, group_samples in sorted(groups.items(), key=lambda item: (-len(item[1]), item[0])):
                outside_samples = all_column_samples - group_samples
                group_n = len(group_samples)
                outside_n = len(outside_samples)
                if group_n == 0 or outside_n == 0:
                    continue
                positive_group_samples = present_samples & group_samples
                positive_outside_samples = present_samples & outside_samples
                a = len(positive_group_samples)
                b = group_n - a
                c = len(positive_outside_samples)
                d = outside_n - c
                prevalence_group = a / group_n if group_n else 0.0
                prevalence_outside = c / outside_n if outside_n else 0.0
                diff = prevalence_group - prevalence_outside
                support_label = _metadata_support_label(group_n, outside_n)
                test_name = "descriptive_only"
                p_value: float | None = None
                if support_label not in {"insufficient_support", "descriptive_only"} and (a + c) >= 3:
                    test_name, p_value = _binary_test_for_counts(a, b, c, d)
                flags = ["multiple_testing", "exploratory_only"]
                if support_label == "insufficient_support":
                    flags.append("small_group_warning")
                elif support_label == "descriptive_only":
                    flags.append("descriptive_only")
                elif support_label == "exploratory":
                    flags.append("exploratory_support")
                if (a + c) < 3:
                    flags.append("low_positive_count")
                if missing_fraction_by_column.get(metadata_column, 0.0) >= 0.30:
                    flags.append("missing_metadata")
                bioproject_flag = _dominance_flag(positive_group_samples, metadata_by_sample, "bioproject", "bioproject_dominance")
                lineage_flag = _dominance_flag(positive_group_samples, metadata_by_sample, "mlst_ST", "lineage_dominance")
                if bioproject_flag:
                    flags.append(bioproject_flag)
                if lineage_flag:
                    flags.append(lineage_flag)
                feature_rows.append({
                    "database": database,
                    "feature_id": feature_id,
                    "feature_category": category_by_feature.get((database, feature_id), ""),
                    "metadata_column": metadata_column,
                    "metadata_group": group,
                    "group_n": str(group_n),
                    "outside_group_n": str(outside_n),
                    "positive_in_group": str(a),
                    "positive_outside_group": str(c),
                    "prevalence_in_group": f"{prevalence_group:.4f}",
                    "prevalence_outside_group": f"{prevalence_outside:.4f}",
                    "prevalence_difference": f"{diff:.4f}",
                    "prevalence_difference_percent": f"{diff * 100:.2f}",
                    "odds_ratio": f"{_odds_ratio(a, b, c, d):.6g}",
                    "test_name": test_name,
                    "p_value": _format_pvalue(p_value),
                    "q_value": "",
                    "effect_size_label": "enriched" if diff > 0 else ("depleted" if diff < 0 else "neutral"),
                    "support_label": support_label,
                    "warning_flags": _flag_string(flags),
                    "interpretation_label": "",
                })
    add_bh_qvalues(feature_rows)
    for row in feature_rows:
        q_value = _float_or_none(row.get("q_value", ""))
        diff = _float_or_none(row.get("prevalence_difference", "")) or 0.0
        row["interpretation_label"] = _metadata_interpretation_label(
            q_value,
            int(row.get("group_n", "0") or 0),
            int(row.get("outside_group_n", "0") or 0),
            diff,
            row.get("warning_flags", ""),
        )
    feature_rows = sorted(
        feature_rows,
        key=lambda row: (
            row.get("interpretation_label", "") not in {"strong_supported", "moderate_supported"},
            -abs(_float_or_none(row.get("prevalence_difference", "")) or 0.0),
            row.get("database", ""),
            row.get("metadata_column", ""),
            row.get("feature_id", ""),
        ),
    )
    feature_path = tables / "metadata_feature_enrichment.tsv"
    feature_fields = [
        "database",
        "feature_id",
        "feature_category",
        "metadata_column",
        "metadata_group",
        "group_n",
        "outside_group_n",
        "positive_in_group",
        "positive_outside_group",
        "prevalence_in_group",
        "prevalence_outside_group",
        "prevalence_difference",
        "prevalence_difference_percent",
        "odds_ratio",
        "test_name",
        "p_value",
        "q_value",
        "effect_size_label",
        "support_label",
        "warning_flags",
        "interpretation_label",
    ]
    write_rows(feature_path, feature_rows, feature_fields)

    features_by_sample_database: dict[tuple[str, str], set[str]] = defaultdict(set)
    features_by_sample_database_category: dict[tuple[str, str, str], set[str]] = defaultdict(set)
    for row in features:
        sample = row.get("assembly_accession", "") or row.get("sample_id", "")
        database = row.get("database", "")
        feature_id = row.get("feature_id", "")
        category = row.get("feature_category", "")
        if sample and database and feature_id:
            features_by_sample_database[(sample, database)].add(feature_id)
            if category:
                features_by_sample_database_category[(sample, database, category)].add(feature_id)
    databases = sorted({database for _, database in features_by_sample_database})
    categories = sorted({(database, category) for _, database, category in features_by_sample_database_category})

    def burden_values(database: str, sample_set: set[str]) -> list[float]:
        return [float(len(features_by_sample_database.get((sample, database), set()))) for sample in sorted(sample_set)]

    def category_values(database: str, category: str, sample_set: set[str]) -> list[float]:
        return [float(len(features_by_sample_database_category.get((sample, database, category), set()))) for sample in sorted(sample_set)]

    burden_rows = []
    category_rows = []
    burden_omnibus_rows = []
    category_omnibus_rows = []
    for metadata_column in metadata_columns:
        groups = groups_for_column(metadata_column)
        all_column_samples = set().union(*groups.values()) if groups else set()
        ordered_groups = sorted(groups.items(), key=lambda item: (-len(item[1]), item[0]))
        for database in databases:
            burden_group_values = []
            for group, group_samples in ordered_groups:
                outside_samples = all_column_samples - group_samples
                group_n = len(group_samples)
                outside_n = len(outside_samples)
                if not group_n or not outside_n:
                    continue
                group_values = burden_values(database, group_samples)
                burden_group_values.append((group, group_values))
                outside_values = burden_values(database, outside_samples)
                support_label = _metadata_support_label(group_n, outside_n)
                p_value = _mann_whitney_u_p_value(group_values, outside_values) if support_label not in {"insufficient_support", "descriptive_only"} else None
                diff = _median(group_values) - _median(outside_values)
                flags = ["multiple_testing", "exploratory_only"]
                if support_label == "insufficient_support":
                    flags.append("small_group_warning")
                elif support_label == "descriptive_only":
                    flags.append("descriptive_only")
                elif support_label == "exploratory":
                    flags.append("exploratory_support")
                if missing_fraction_by_column.get(metadata_column, 0.0) >= 0.30:
                    flags.append("missing_metadata")
                burden_rows.append({
                    "database": database,
                    "metadata_column": metadata_column,
                    "metadata_group": group,
                    "group_n": str(group_n),
                    "outside_group_n": str(outside_n),
                    "median_burden_group": f"{_median(group_values):.4f}",
                    "median_burden_outside": f"{_median(outside_values):.4f}",
                    "mean_burden_group": f"{_mean(group_values):.4f}",
                    "mean_burden_outside": f"{_mean(outside_values):.4f}",
                    "burden_difference": f"{diff:.4f}",
                    "test_name": "mann_whitney_u" if p_value is not None else "descriptive_only",
                    "p_value": _format_pvalue(p_value),
                    "q_value": "",
                    "support_label": support_label,
                    "warning_flags": _flag_string(flags),
                    "interpretation_label": "",
                })
            if len(burden_group_values) >= 2:
                group_counts = [len(values) for _, values in burden_group_values]
                support_label = _metadata_support_label_for_counts(group_counts)
                test_result = _kruskal_wallis([values for _, values in burden_group_values]) if support_label not in {"insufficient_support", "descriptive_only"} else None
                medians = {group: _median(values) for group, values in burden_group_values}
                means = {group: _mean(values) for group, values in burden_group_values}
                median_values = list(medians.values())
                burden_range = (max(median_values) - min(median_values)) if median_values else 0.0
                flags = ["multiple_testing", "exploratory_only"]
                if support_label == "insufficient_support":
                    flags.append("insufficient_group_size")
                elif support_label == "descriptive_only":
                    flags.append("descriptive_only")
                elif support_label == "exploratory":
                    flags.append("exploratory_support")
                if missing_fraction_by_column.get(metadata_column, 0.0) >= 0.30:
                    flags.append("missing_metadata")
                statistic, p_value = test_result if test_result else (None, None)
                burden_omnibus_rows.append({
                    "database": database,
                    "metadata_column": metadata_column,
                    "groups_tested": str(len(burden_group_values)),
                    "samples_tested": str(sum(group_counts)),
                    "group_sizes": ";".join(f"{group}:{len(values)}" for group, values in burden_group_values),
                    "median_burden_by_group": ";".join(f"{group}:{medians[group]:.4f}" for group, _ in burden_group_values),
                    "mean_burden_by_group": ";".join(f"{group}:{means[group]:.4f}" for group, _ in burden_group_values),
                    "burden_range_median": f"{burden_range:.4f}",
                    "test_name": "kruskal_wallis" if test_result else "descriptive_only",
                    "test_statistic": f"{statistic:.6g}" if statistic is not None else "",
                    "p_value": _format_pvalue(p_value),
                    "q_value": "",
                    "support_label": support_label,
                    "warning_flags": _flag_string(flags),
                    "interpretation_label": "",
                })
        for database, category in categories:
            category_group_values = []
            for group, group_samples in ordered_groups:
                outside_samples = all_column_samples - group_samples
                group_n = len(group_samples)
                outside_n = len(outside_samples)
                if not group_n or not outside_n:
                    continue
                group_values = category_values(database, category, group_samples)
                category_group_values.append((group, group_values))
                outside_values = category_values(database, category, outside_samples)
                support_label = _metadata_support_label(group_n, outside_n)
                p_value = _mann_whitney_u_p_value(group_values, outside_values) if support_label not in {"insufficient_support", "descriptive_only"} else None
                diff = _median(group_values) - _median(outside_values)
                flags = ["multiple_testing", "exploratory_only"]
                if support_label == "insufficient_support":
                    flags.append("small_group_warning")
                elif support_label == "descriptive_only":
                    flags.append("descriptive_only")
                elif support_label == "exploratory":
                    flags.append("exploratory_support")
                category_rows.append({
                    "database": database,
                    "feature_category": category,
                    "metadata_column": metadata_column,
                    "metadata_group": group,
                    "group_n": str(group_n),
                    "outside_group_n": str(outside_n),
                    "median_burden_group": f"{_median(group_values):.4f}",
                    "median_burden_outside": f"{_median(outside_values):.4f}",
                    "mean_burden_group": f"{_mean(group_values):.4f}",
                    "mean_burden_outside": f"{_mean(outside_values):.4f}",
                    "burden_difference": f"{diff:.4f}",
                    "test_name": "mann_whitney_u" if p_value is not None else "descriptive_only",
                    "p_value": _format_pvalue(p_value),
                    "q_value": "",
                    "support_label": support_label,
                    "warning_flags": _flag_string(flags),
                    "interpretation_label": "",
                })
            if len(category_group_values) >= 2:
                group_counts = [len(values) for _, values in category_group_values]
                support_label = _metadata_support_label_for_counts(group_counts)
                test_result = _kruskal_wallis([values for _, values in category_group_values]) if support_label not in {"insufficient_support", "descriptive_only"} else None
                medians = {group: _median(values) for group, values in category_group_values}
                means = {group: _mean(values) for group, values in category_group_values}
                median_values = list(medians.values())
                burden_range = (max(median_values) - min(median_values)) if median_values else 0.0
                flags = ["multiple_testing", "exploratory_only"]
                if support_label == "insufficient_support":
                    flags.append("insufficient_group_size")
                elif support_label == "descriptive_only":
                    flags.append("descriptive_only")
                elif support_label == "exploratory":
                    flags.append("exploratory_support")
                if missing_fraction_by_column.get(metadata_column, 0.0) >= 0.30:
                    flags.append("missing_metadata")
                statistic, p_value = test_result if test_result else (None, None)
                category_omnibus_rows.append({
                    "database": database,
                    "feature_category": category,
                    "metadata_column": metadata_column,
                    "groups_tested": str(len(category_group_values)),
                    "samples_tested": str(sum(group_counts)),
                    "group_sizes": ";".join(f"{group}:{len(values)}" for group, values in category_group_values),
                    "median_burden_by_group": ";".join(f"{group}:{medians[group]:.4f}" for group, _ in category_group_values),
                    "mean_burden_by_group": ";".join(f"{group}:{means[group]:.4f}" for group, _ in category_group_values),
                    "burden_range_median": f"{burden_range:.4f}",
                    "test_name": "kruskal_wallis" if test_result else "descriptive_only",
                    "test_statistic": f"{statistic:.6g}" if statistic is not None else "",
                    "p_value": _format_pvalue(p_value),
                    "q_value": "",
                    "support_label": support_label,
                    "warning_flags": _flag_string(flags),
                    "interpretation_label": "",
                })
    add_bh_qvalues(burden_rows)
    add_bh_qvalues(category_rows)
    add_bh_qvalues(burden_omnibus_rows)
    add_bh_qvalues(category_omnibus_rows)
    for row in burden_rows + category_rows:
        q_value = _float_or_none(row.get("q_value", ""))
        diff = _float_or_none(row.get("burden_difference", "")) or 0.0
        row["interpretation_label"] = _metadata_interpretation_label(
            q_value,
            int(row.get("group_n", "0") or 0),
            int(row.get("outside_group_n", "0") or 0),
            diff,
            row.get("warning_flags", ""),
        )
    for row in burden_omnibus_rows + category_omnibus_rows:
        q_value = _float_or_none(row.get("q_value", ""))
        effect = _float_or_none(row.get("burden_range_median", "")) or 0.0
        flags = {flag for flag in row.get("warning_flags", "").split(";") if flag}
        if row.get("support_label") == "insufficient_support":
            row["interpretation_label"] = "insufficient_support"
        elif row.get("support_label") == "descriptive_only":
            row["interpretation_label"] = "descriptive_only"
        elif row.get("support_label") == "exploratory":
            row["interpretation_label"] = "exploratory"
        elif q_value is not None and q_value <= 0.05 and effect >= 1.0 and not flags.intersection({"missing_metadata", "bioproject_dominance", "lineage_dominance"}):
            row["interpretation_label"] = "strong_supported"
        elif q_value is not None and q_value <= 0.10 and effect >= 0.5:
            row["interpretation_label"] = "moderate_supported"
        elif q_value is not None and q_value <= 0.10:
            row["interpretation_label"] = "exploratory"
        else:
            row["interpretation_label"] = "descriptive_only"
    burden_rows = sorted(burden_rows, key=lambda row: (-abs(_float_or_none(row.get("burden_difference", "")) or 0.0), row.get("database", ""), row.get("metadata_column", "")))
    category_rows = sorted(category_rows, key=lambda row: (-abs(_float_or_none(row.get("burden_difference", "")) or 0.0), row.get("database", ""), row.get("feature_category", "")))
    burden_omnibus_rows = sorted(burden_omnibus_rows, key=lambda row: (-abs(_float_or_none(row.get("burden_range_median", "")) or 0.0), row.get("database", ""), row.get("metadata_column", "")))
    category_omnibus_rows = sorted(category_omnibus_rows, key=lambda row: (-abs(_float_or_none(row.get("burden_range_median", "")) or 0.0), row.get("database", ""), row.get("feature_category", "")))

    burden_path = tables / "metadata_burden_associations.tsv"
    burden_fields = [
        "database",
        "metadata_column",
        "metadata_group",
        "group_n",
        "outside_group_n",
        "median_burden_group",
        "median_burden_outside",
        "mean_burden_group",
        "mean_burden_outside",
        "burden_difference",
        "test_name",
        "p_value",
        "q_value",
        "support_label",
        "warning_flags",
        "interpretation_label",
    ]
    write_rows(burden_path, burden_rows, burden_fields)

    category_path = tables / "metadata_category_enrichment.tsv"
    category_fields = ["database", "feature_category", *burden_fields[1:]]
    write_rows(category_path, category_rows, category_fields)

    omnibus_fields = [
        "database",
        "metadata_column",
        "groups_tested",
        "samples_tested",
        "group_sizes",
        "median_burden_by_group",
        "mean_burden_by_group",
        "burden_range_median",
        "test_name",
        "test_statistic",
        "p_value",
        "q_value",
        "support_label",
        "warning_flags",
        "interpretation_label",
    ]
    burden_omnibus_path = tables / "metadata_burden_omnibus.tsv"
    category_omnibus_path = tables / "metadata_category_omnibus.tsv"
    write_rows(burden_omnibus_path, burden_omnibus_rows, omnibus_fields)
    write_rows(category_omnibus_path, category_omnibus_rows, ["database", "feature_category", *omnibus_fields[1:]])

    default_database = "amr" if any(row.get("database") == "amr" for row in feature_rows) else (feature_rows[0]["database"] if feature_rows else (databases[0] if databases else "all"))
    default_column = "isolation_source" if "isolation_source" in metadata_columns else ("country" if "country" in metadata_columns else (metadata_columns[0] if metadata_columns else "metadata"))
    default_group = ""
    for row in feature_rows:
        if row.get("database") == default_database and row.get("metadata_column") == default_column:
            default_group = row.get("metadata_group", "")
            break
    if not default_group and default_column in metadata_columns:
        group_counts = groups_for_column(default_column)
        default_group = max(group_counts.items(), key=lambda item: len(item[1]))[0] if group_counts else ""

    volcano_rows = [
        row for row in feature_rows
        if row.get("database") == default_database and row.get("metadata_column") == default_column and row.get("metadata_group") == default_group
    ][:max_features_per_database]
    volcano_base = f"metadata_volcano_{_safe_filename(default_database)}_{_safe_filename(default_column)}_{_safe_filename(default_group)}"
    volcano_data = figures / f"{volcano_base}.data.tsv"
    volcano_svg = figures / f"{volcano_base}.svg"
    volcano_png = figures / f"{volcano_base}.png"
    volcano_pdf = figures / f"{volcano_base}.pdf"
    write_rows(volcano_data, volcano_rows, feature_fields)
    _write_metadata_volcano_svg(volcano_svg, volcano_rows, f"{default_database} Enrichment In {default_column}={default_group}")
    _write_metadata_volcano_png(volcano_png, volcano_rows)
    _write_simple_pdf(
        volcano_pdf,
        f"{default_database} Enrichment In {default_column}={default_group}",
        [
            f"{row.get('feature_id', '')}: diff={row.get('prevalence_difference', '')}, q={row.get('q_value', '')}, {row.get('interpretation_label', '')}"
            for row in volcano_rows[:30]
        ] or ["No feature-level metadata enrichment rows were available."],
    )

    heatmap_features = {
        row.get("feature_id", "") for row in sorted(
            [row for row in feature_rows if row.get("database") == default_database and row.get("metadata_column") == default_column],
            key=lambda row: -abs(_float_or_none(row.get("prevalence_difference", "")) or 0.0),
        )[:top_n]
    }
    heatmap_rows = [
        {
            **row,
            "feature_label": f"{row.get('database', '')}:{row.get('feature_id', '')}",
        }
        for row in feature_rows
        if row.get("database") == default_database and row.get("metadata_column") == default_column and row.get("feature_id") in heatmap_features
    ]
    heatmap_base = f"metadata_enrichment_heatmap_{_safe_filename(default_database)}_{_safe_filename(default_column)}"
    heatmap_data = figures / f"{heatmap_base}.data.tsv"
    heatmap_svg = figures / f"{heatmap_base}.svg"
    heatmap_png = figures / f"{heatmap_base}.png"
    heatmap_pdf = figures / f"{heatmap_base}.pdf"
    write_rows(heatmap_data, heatmap_rows, [*feature_fields, "feature_label"])
    _write_diverging_heatmap_svg(heatmap_svg, heatmap_rows, f"{default_database} Metadata Enrichment Heatmap", "feature_label", "metadata_group", "prevalence_difference")
    _write_diverging_heatmap_png(heatmap_png, heatmap_rows, "feature_label", "metadata_group", "prevalence_difference")
    _write_simple_pdf(
        heatmap_pdf,
        f"{default_database} Metadata Enrichment Heatmap",
        [f"{row.get('feature_label', '')} / {row.get('metadata_group', '')}: {row.get('prevalence_difference', '')}" for row in heatmap_rows[:34]]
        or ["No metadata enrichment heatmap rows were available."],
    )

    boxplot_groups = []
    if default_column in metadata_columns and default_database in databases:
        for group, group_samples in sorted(groups_for_column(default_column).items(), key=lambda item: (-len(item[1]), item[0]))[:20]:
            values = burden_values(default_database, group_samples)
            boxplot_groups.append({
                "metadata_group": group,
                "n": str(len(values)),
                "min": f"{(min(values) if values else 0.0):.4f}",
                "q1": f"{_quantile(values, 0.25):.4f}",
                "median": f"{_median(values):.4f}",
                "q3": f"{_quantile(values, 0.75):.4f}",
                "max": f"{(max(values) if values else 0.0):.4f}",
                "mean": f"{_mean(values):.4f}",
            })
    boxplot_base = f"metadata_burden_boxplot_{_safe_filename(default_database)}_{_safe_filename(default_column)}"
    boxplot_data = figures / f"{boxplot_base}.data.tsv"
    boxplot_svg = figures / f"{boxplot_base}.svg"
    boxplot_png = figures / f"{boxplot_base}.png"
    boxplot_pdf = figures / f"{boxplot_base}.pdf"
    write_rows(boxplot_data, boxplot_groups, ["metadata_group", "n", "min", "q1", "median", "q3", "max", "mean"])
    _write_burden_boxplot_svg(boxplot_svg, boxplot_groups, f"{default_database} Burden By {default_column}")
    _write_burden_boxplot_png(boxplot_png, boxplot_groups)
    _write_simple_pdf(
        boxplot_pdf,
        f"{default_database} Burden By {default_column}",
        [f"{row.get('metadata_group', '')}: median={row.get('median', '')}, n={row.get('n', '')}" for row in boxplot_groups],
    )

    category_plot_rows = [
        {
            **row,
            "category_label": f"{row.get('feature_category', '')} / {row.get('metadata_group', '')}",
            "abs_burden_difference": f"{abs(_float_or_none(row.get('burden_difference', '')) or 0.0):.4f}",
        }
        for row in category_rows
        if row.get("database") == default_database and row.get("metadata_column") == default_column
    ][:top_n]
    category_base = f"metadata_category_enrichment_{_safe_filename(default_database)}_{_safe_filename(default_column)}"
    category_data = figures / f"{category_base}.data.tsv"
    category_svg = figures / f"{category_base}.svg"
    category_png = figures / f"{category_base}.png"
    category_pdf = figures / f"{category_base}.pdf"
    write_rows(category_data, category_plot_rows, [*category_fields, "category_label", "abs_burden_difference"])
    _write_bar_svg(category_svg, category_plot_rows, f"{default_database} Category Enrichment", "category_label", "abs_burden_difference", "Absolute median burden difference")
    _write_bar_png(category_png, category_plot_rows, "abs_burden_difference")
    _write_simple_pdf(
        category_pdf,
        f"{default_database} Category Enrichment",
        [f"{row.get('category_label', '')}: {row.get('burden_difference', '')}" for row in category_plot_rows],
    )

    strong_features = sum(1 for row in feature_rows if row.get("interpretation_label") == "strong_supported")
    moderate_features = sum(1 for row in feature_rows if row.get("interpretation_label") == "moderate_supported")
    exploratory_features = sum(1 for row in feature_rows if row.get("interpretation_label") == "exploratory")
    strong_burdens = sum(1 for row in burden_rows if row.get("interpretation_label") == "strong_supported")
    strong_omnibus = sum(1 for row in burden_omnibus_rows + category_omnibus_rows if row.get("interpretation_label") == "strong_supported")
    warning_count = sum(1 for row in feature_rows + burden_rows + category_rows + burden_omnibus_rows + category_omnibus_rows if row.get("warning_flags"))
    usable_metadata_columns = sum(1 for row in usability_rows if row.get("recommended_use") == "association_testing")
    sparse_metadata_columns = sum(1 for row in usability_rows if row.get("recommended_use") == "descriptive_or_bias_check")
    excluded_metadata_columns_count = sum(1 for row in usability_rows if row.get("recommended_use") == "exclude")
    summary_rows = [
        {"metric": "samples", "value": str(sample_count), "message": "Samples represented in the important metadata association screen."},
        {"metric": "metadata_columns_screened", "value": str(len(metadata_columns)), "message": "Metadata columns with at least two non-missing groups."},
        {"metric": "metadata_columns_usable", "value": str(usable_metadata_columns), "message": "Metadata columns recommended for association testing in the usability summary."},
        {"metric": "metadata_columns_sparse_or_biased", "value": str(sparse_metadata_columns), "message": "Metadata columns recommended only for descriptive review or bias checks."},
        {"metric": "metadata_columns_excluded", "value": str(excluded_metadata_columns_count), "message": "Metadata columns excluded from testing because they lack usable group structure."},
        {"metric": "feature_enrichment_rows", "value": str(len(feature_rows)), "message": "Feature-by-group comparisons in the report-facing metadata enrichment table."},
        {"metric": "strong_feature_associations", "value": str(strong_features), "message": "Feature rows with strong_supported interpretation labels."},
        {"metric": "moderate_feature_associations", "value": str(moderate_features), "message": "Feature rows with moderate_supported interpretation labels."},
        {"metric": "exploratory_feature_associations", "value": str(exploratory_features), "message": "Feature rows retained as exploratory signals."},
        {"metric": "database_burden_associations", "value": str(len(burden_rows)), "message": "Database-burden group comparisons."},
        {"metric": "strong_burden_associations", "value": str(strong_burdens), "message": "Database-burden rows with strong_supported interpretation labels."},
        {"metric": "burden_omnibus_tests", "value": str(len(burden_omnibus_rows)), "message": "Kruskal-Wallis multi-group database-burden tests."},
        {"metric": "category_omnibus_tests", "value": str(len(category_omnibus_rows)), "message": "Kruskal-Wallis multi-group category-burden tests."},
        {"metric": "strong_omnibus_tests", "value": str(strong_omnibus), "message": "Omnibus rows with strong_supported interpretation labels."},
        {"metric": "warning_rows", "value": str(warning_count), "message": "Rows carrying at least one caution flag."},
        {"metric": "default_view", "value": f"{default_database}|{default_column}|{default_group}", "message": "Default database, metadata variable, and group shown in figures."},
        {"metric": "feature_cap", "value": str(max_features_per_database), "message": f"Report-facing feature screening cap per database. Capped databases: {','.join(sorted(capped_databases)) or 'none'}."},
    ]
    summary_path = tables / "metadata_association_summary.tsv"
    write_rows(summary_path, summary_rows, ["metric", "value", "message"])

    metadata_tables_zip = important_dir / "metadata_association_tables.zip"
    metadata_figures_zip = important_dir / "metadata_association_figures.zip"
    _write_zip_bundle(metadata_tables_zip, [usability_path, feature_path, burden_path, category_path, burden_omnibus_path, category_omnibus_path, summary_path], important_dir)
    _write_zip_bundle(
        metadata_figures_zip,
        [
            volcano_svg, volcano_png, volcano_pdf, volcano_data,
            heatmap_svg, heatmap_png, heatmap_pdf, heatmap_data,
            boxplot_svg, boxplot_png, boxplot_pdf, boxplot_data,
            category_svg, category_png, category_pdf, category_data,
        ],
        important_dir,
    )

    interactive_html = figures / "metadata_associations.html"
    interactive_html.write_text(
        f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Metadata Associations</title>
<style>
body {{ font-family: Arial, sans-serif; color: #1f2933; margin: 1.5rem; }}
label {{ font-weight: 700; margin-right: 0.35rem; }}
select {{ margin: 0 1rem 0.75rem 0; padding: 0.35rem; }}
.warning {{ background: #fff7ed; border-left: 4px solid #c2410c; padding: 0.75rem; margin: 0.75rem 0; }}
.panel {{ border: 1px solid #d9e2ec; background: #f8fafc; padding: 0.75rem; margin: 0.75rem 0; }}
.figure-row {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 1rem; }}
.figure-row img {{ max-width: 100%; border: 1px solid #d9e2ec; background: white; }}
table {{ border-collapse: collapse; width: 100%; font-size: 0.88rem; }}
th, td {{ border: 1px solid #d9e2ec; padding: 0.35rem; text-align: left; }}
th {{ background: #f0f4f8; }}
</style></head><body>
<h1>Metadata Associations</h1>
<div class="warning">Metadata associations are exploratory and may reflect sampling, BioProject structure, lineage composition, geography, collection year, or missing metadata. Enriched/depleted means more/less prevalent in this dataset, not causal.</div>
<label for="database">Database</label><select id="database"></select>
<label for="associationType">Association type</label><select id="associationType"><option value="feature" selected>Feature prevalence</option><option value="burden">Database burden</option><option value="category">Category burden</option></select>
<label for="metadataColumn">Metadata variable</label><select id="metadataColumn"></select>
<label for="metadataGroup">Group</label><select id="metadataGroup"></select>
<label for="minGroup">Minimum group size</label><select id="minGroup"><option value="0">All</option><option value="5">n >= 5</option><option value="10" selected>n >= 10</option><option value="30">n >= 30</option></select>
<label for="significance">Significance</label><select id="significance"><option value="q0.05">q < 0.05</option><option value="q0.10" selected>q < 0.10</option><option value="p0.05">p < 0.05</option><option value="all">show all</option></select>
<label for="effectSize">Effect size</label><select id="effectSize"><option value="any" selected>any</option><option value="diff10">difference >= 10%</option><option value="diff20">difference >= 20%</option><option value="or2">odds ratio >= 2 or <= 0.5</option><option value="or5">odds ratio >= 5 or <= 0.2</option><option value="burden1">burden difference >= 1</option><option value="burden5">burden difference >= 5</option></select>
<label for="warningFilter">Warning filter</label><select id="warningFilter"><option value="all" selected>show all</option><option value="hide_weak">hide weak support</option><option value="strong">strong only</option></select>
<label for="displayCount">Display</label><select id="displayCount"><option value="10">Top 10</option><option value="20" selected>Top 20</option><option value="50">Top 50</option><option value="99999">Complete</option></select>
<div id="summary"></div>
<div class="panel"><h2>Default Visuals</h2><div class="figure-row">
<div><h3>Volcano plot</h3><img src="{html.escape(volcano_base)}.svg" alt="Metadata volcano plot"><p><a href="{html.escape(volcano_base)}.png">PNG</a> | <a href="{html.escape(volcano_base)}.svg">SVG</a> | <a href="{html.escape(volcano_base)}.pdf">PDF</a> | <a href="{html.escape(volcano_base)}.data.tsv">Data TSV</a></p></div>
<div><h3>Enrichment heatmap</h3><img src="{html.escape(heatmap_base)}.svg" alt="Metadata enrichment heatmap"><p><a href="{html.escape(heatmap_base)}.png">PNG</a> | <a href="{html.escape(heatmap_base)}.svg">SVG</a> | <a href="{html.escape(heatmap_base)}.pdf">PDF</a> | <a href="{html.escape(heatmap_base)}.data.tsv">Data TSV</a></p></div>
<div><h3>Burden boxplot</h3><img src="{html.escape(boxplot_base)}.svg" alt="Metadata burden boxplot"><p><a href="{html.escape(boxplot_base)}.png">PNG</a> | <a href="{html.escape(boxplot_base)}.svg">SVG</a> | <a href="{html.escape(boxplot_base)}.pdf">PDF</a> | <a href="{html.escape(boxplot_base)}.data.tsv">Data TSV</a></p></div>
<div><h3>Category enrichment</h3><img src="{html.escape(category_base)}.svg" alt="Metadata category enrichment"><p><a href="{html.escape(category_base)}.png">PNG</a> | <a href="{html.escape(category_base)}.svg">SVG</a> | <a href="{html.escape(category_base)}.pdf">PDF</a> | <a href="{html.escape(category_base)}.data.tsv">Data TSV</a></p></div>
</div></div>
<h2>Filtered Results</h2><div id="table"></div>
<p><a href="../tables/metadata_usability_summary.tsv">Download metadata usability summary</a> | <a href="../tables/metadata_feature_enrichment.tsv">Download feature enrichment table</a> | <a href="../tables/metadata_burden_associations.tsv">Download burden associations</a> | <a href="../tables/metadata_category_enrichment.tsv">Download category enrichment</a> | <a href="../tables/metadata_burden_omnibus.tsv">Download burden omnibus tests</a> | <a href="../tables/metadata_category_omnibus.tsv">Download category omnibus tests</a> | <a href="../metadata_association_tables.zip">Download all metadata association tables ZIP</a> | <a href="../metadata_association_figures.zip">Download all metadata association figures ZIP</a></p>
<script>
const featureRows = {json.dumps(feature_rows[:4000])};
const burdenRows = {json.dumps(burden_rows[:3000])};
const categoryRows = {json.dumps(category_rows[:3000])};
const burdenOmnibusRows = {json.dumps(burden_omnibus_rows[:1000])};
const categoryOmnibusRows = {json.dumps(category_omnibus_rows[:1000])};
const defaultDatabase = {json.dumps(default_database)};
const defaultColumn = {json.dumps(default_column)};
const defaultGroup = {json.dumps(default_group)};
function num(value) {{ const n = Number(value); return Number.isFinite(n) ? n : 0; }}
const database = document.getElementById('database'), associationType = document.getElementById('associationType'), metadataColumn = document.getElementById('metadataColumn'), metadataGroup = document.getElementById('metadataGroup'), minGroup = document.getElementById('minGroup'), significance = document.getElementById('significance'), effectSize = document.getElementById('effectSize'), warningFilter = document.getElementById('warningFilter'), displayCount = document.getElementById('displayCount');
function currentRows() {{ return associationType.value === 'burden' ? burdenRows : (associationType.value === 'category' ? categoryRows : featureRows); }}
function fillSelect(select, values, preferred) {{ select.innerHTML = values.map(v => `<option value="${{v}}">${{v}}</option>`).join(''); if (values.includes(preferred)) select.value = preferred; }}
function refreshControls() {{
  const rows = currentRows();
  fillSelect(database, Array.from(new Set(rows.map(r => r.database))).filter(Boolean).sort(), defaultDatabase);
  fillSelect(metadataColumn, Array.from(new Set(rows.map(r => r.metadata_column))).filter(Boolean).sort(), defaultColumn);
  refreshGroups();
}}
function refreshGroups() {{
  const rows = currentRows().filter(r => r.database === database.value && r.metadata_column === metadataColumn.value);
  fillSelect(metadataGroup, Array.from(new Set(rows.map(r => r.metadata_group))).filter(Boolean).sort(), defaultGroup);
}}
function passes(row) {{
  if (row.database !== database.value || row.metadata_column !== metadataColumn.value || row.metadata_group !== metadataGroup.value) return false;
  if (num(row.group_n) < num(minGroup.value)) return false;
  if (significance.value === 'q0.05' && !(num(row.q_value) > 0 && num(row.q_value) < 0.05)) return false;
  if (significance.value === 'q0.10' && !(num(row.q_value) > 0 && num(row.q_value) < 0.10)) return false;
  if (significance.value === 'p0.05' && !(num(row.p_value) > 0 && num(row.p_value) < 0.05)) return false;
  const diff = Math.abs(num(row.prevalence_difference || row.burden_difference));
  const odds = num(row.odds_ratio);
  if (effectSize.value === 'diff10' && diff < 0.10) return false;
  if (effectSize.value === 'diff20' && diff < 0.20) return false;
  if (effectSize.value === 'or2' && row.odds_ratio && !(odds >= 2 || (odds > 0 && odds <= 0.5))) return false;
  if (effectSize.value === 'or5' && row.odds_ratio && !(odds >= 5 || (odds > 0 && odds <= 0.2))) return false;
  if (effectSize.value === 'burden1' && Math.abs(num(row.burden_difference)) < 1) return false;
  if (effectSize.value === 'burden5' && Math.abs(num(row.burden_difference)) < 5) return false;
  if (warningFilter.value === 'hide_weak' && ['insufficient_support','descriptive_only'].includes(row.interpretation_label)) return false;
  if (warningFilter.value === 'strong' && row.interpretation_label !== 'strong_supported') return false;
  return true;
}}
function render() {{
  const rows = currentRows().filter(passes).slice(0, Number(displayCount.value));
  const cols = associationType.value === 'feature'
    ? ['database','feature_id','metadata_column','metadata_group','group_n','positive_in_group','prevalence_in_group','prevalence_difference','odds_ratio','q_value','effect_size_label','support_label','interpretation_label','warning_flags']
    : (associationType.value === 'category'
      ? ['database','feature_category','metadata_column','metadata_group','group_n','median_burden_group','median_burden_outside','burden_difference','q_value','support_label','interpretation_label','warning_flags']
      : ['database','metadata_column','metadata_group','group_n','median_burden_group','median_burden_outside','burden_difference','q_value','support_label','interpretation_label','warning_flags']);
  const omnibus = associationType.value === 'category' ? categoryOmnibusRows : burdenOmnibusRows;
  const omnibusMatches = omnibus.filter(r => r.database === database.value && r.metadata_column === metadataColumn.value).slice(0, 10);
  const omnibusText = associationType.value === 'feature' ? '' : `<p>Omnibus Kruskal-Wallis tests for this view: ${{omnibusMatches.map(r => `${{r.test_name}} q=${{r.q_value || 'NA'}} support=${{r.support_label}}`).join('; ') || 'none available'}}.</p>`;
  document.getElementById('summary').innerHTML = `<p>${{rows.length}} rows shown after filters. Statistical labels use tiered sample-size support; rows with q-values are FDR-corrected within this report-facing screen.</p>${{omnibusText}}`;
  document.getElementById('table').innerHTML = '<table><thead><tr>' + cols.map(c => `<th>${{c}}</th>`).join('') + '</tr></thead><tbody>' + rows.map(r => '<tr>' + cols.map(c => `<td>${{r[c] || ''}}</td>`).join('') + '</tr>').join('') + '</tbody></table>';
}}
associationType.addEventListener('change', () => {{ refreshControls(); render(); }});
database.addEventListener('change', () => {{ refreshGroups(); render(); }});
metadataColumn.addEventListener('change', () => {{ refreshGroups(); render(); }});
[metadataGroup, minGroup, significance, effectSize, warningFilter, displayCount].forEach(el => el.addEventListener('change', render));
refreshControls(); render();
</script></body></html>
""",
        encoding="utf-8",
    )

    return {
        "important_metadata_usability_summary": str(usability_path),
        "important_metadata_feature_enrichment": str(feature_path),
        "important_metadata_burden_associations": str(burden_path),
        "important_metadata_category_enrichment": str(category_path),
        "important_metadata_burden_omnibus": str(burden_omnibus_path),
        "important_metadata_category_omnibus": str(category_omnibus_path),
        "important_metadata_association_summary": str(summary_path),
        "important_metadata_volcano_svg": str(volcano_svg),
        "important_metadata_volcano_png": str(volcano_png),
        "important_metadata_volcano_pdf": str(volcano_pdf),
        "important_metadata_volcano_data": str(volcano_data),
        "important_metadata_enrichment_heatmap_svg": str(heatmap_svg),
        "important_metadata_enrichment_heatmap_png": str(heatmap_png),
        "important_metadata_enrichment_heatmap_pdf": str(heatmap_pdf),
        "important_metadata_enrichment_heatmap_data": str(heatmap_data),
        "important_metadata_burden_boxplot_svg": str(boxplot_svg),
        "important_metadata_burden_boxplot_png": str(boxplot_png),
        "important_metadata_burden_boxplot_pdf": str(boxplot_pdf),
        "important_metadata_burden_boxplot_data": str(boxplot_data),
        "important_metadata_category_enrichment_svg": str(category_svg),
        "important_metadata_category_enrichment_png": str(category_png),
        "important_metadata_category_enrichment_pdf": str(category_pdf),
        "important_metadata_category_enrichment_data": str(category_data),
        "important_metadata_associations_html": str(interactive_html),
        "important_metadata_association_tables_zip": str(metadata_tables_zip),
        "important_metadata_association_figures_zip": str(metadata_figures_zip),
    }


def _lineage_support_label(lineage_n: int, outside_n: int | None = None) -> str:
    minimum = lineage_n if outside_n is None else min(lineage_n, outside_n)
    if minimum < 5:
        return "insufficient_support"
    if minimum < 10:
        return "descriptive_only"
    if minimum < 30:
        return "exploratory"
    if minimum >= 100:
        return "strong_support"
    return "standard_support"


def _lineage_dominance_label(fraction: float | None) -> str:
    if fraction is None:
        return "insufficient_lineage_data"
    if fraction >= 0.90:
        return "severe_lineage_confounding"
    if fraction >= 0.70:
        return "strong_lineage_confounding"
    if fraction >= 0.50:
        return "possible_lineage_confounding"
    return "not_lineage_dominated"


def _lineage_interpretation_label(
    q_value: float | None,
    lineage_n: int,
    outside_n: int,
    prevalence_in: float,
    prevalence_outside: float,
    warning_flags: str,
) -> str:
    support = _lineage_support_label(lineage_n, outside_n)
    flags = {flag for flag in str(warning_flags or "").split(";") if flag}
    if support == "insufficient_support":
        return "insufficient_support"
    if support == "descriptive_only":
        return "descriptive_only"
    if prevalence_in >= 0.80 and prevalence_outside <= 0.20 and lineage_n >= 10:
        return "feature_lineage_specific"
    if q_value is not None and q_value <= 0.05 and abs(prevalence_in - prevalence_outside) >= 0.20 and not {"BioProject_lineage_overlap", "small_lineage_group"}.intersection(flags):
        return "feature_lineage_enriched"
    if support == "exploratory" or (q_value is not None and q_value <= 0.10):
        return "exploratory"
    return "descriptive_only"


def write_important_lineage_outputs(
    sample_dir: Path,
    out_dir: Path,
    important_dir: Path,
    top_n: int = 20,
    max_features_per_database: int = 200,
) -> dict[str, str]:
    tables = important_dir / "tables"
    figures = important_dir / "figures"
    tables.mkdir(parents=True, exist_ok=True)
    figures.mkdir(parents=True, exist_ok=True)

    features = read_table(out_dir / "features" / "all_features.tsv")
    metadata_rows = read_table(sample_dir / "basic" / "enriched_genome_dataset.csv")
    if not metadata_rows:
        metadata_rows = normalize_metadata_rows(load_metadata_rows(sample_dir))
    metadata_by_sample = {row.get("assembly_accession", ""): row for row in metadata_rows if row.get("assembly_accession")}
    samples = sorted(set(metadata_by_sample) | {row.get("assembly_accession", "") or row.get("sample_id", "") for row in features if row.get("assembly_accession") or row.get("sample_id")})
    for sample in samples:
        metadata_by_sample.setdefault(sample, {"assembly_accession": sample, "sample_id": sample})
    sample_count = len(samples)

    def clean_metadata_value(value: str) -> str:
        text = str(value or "").strip()
        return "" if is_missing_value(text) else text

    def extract_mlst_st(feature_id: str) -> str:
        text = str(feature_id or "").strip()
        if is_placeholder_mlst_feature(text):
            return ""
        if text.startswith("ST_"):
            return text
        match = re.search(r"ST[-_: ]*([A-Za-z0-9]+)", text, flags=re.IGNORECASE)
        if match:
            st_value = match.group(1).strip("_:- ")
            if st_value and not is_missing_value(st_value):
                return f"ST_{st_value}"
        return ""

    mlst_by_sample: dict[str, str] = {}
    for sample, row in metadata_by_sample.items():
        value = clean_metadata_value(row.get("mlst_ST", ""))
        if value:
            mlst_by_sample[sample] = value
    for row in features:
        if row.get("database") != "mlst" or row.get("presence", "1") != "1":
            continue
        sample = row.get("assembly_accession", "") or row.get("sample_id", "")
        feature_id = row.get("feature_id", "")
        category = row.get("feature_category", "")
        st_value = extract_mlst_st(feature_id)
        if sample and st_value and (category == "sequence_type" or "ST" in feature_id.upper()):
            mlst_by_sample.setdefault(sample, st_value)

    ani_by_sample = {sample: clean_metadata_value(row.get("ani_cluster", "")) for sample, row in metadata_by_sample.items() if clean_metadata_value(row.get("ani_cluster", ""))}
    ani_by_sample.update({sample: cluster for sample, cluster in _ani_cluster_by_sample(out_dir).items() if cluster})
    bioproject_by_sample = {sample: clean_metadata_value(row.get("bioproject", "")) for sample, row in metadata_by_sample.items() if clean_metadata_value(row.get("bioproject", ""))}
    combined_by_sample = {
        sample: mlst_by_sample.get(sample) or ani_by_sample.get(sample) or bioproject_by_sample.get(sample, "")
        for sample in samples
    }
    lineage_maps = {
        "mlst_ST": mlst_by_sample,
        "ani_cluster": ani_by_sample,
        "bioproject": bioproject_by_sample,
        "combined_lineage_label": {sample: value for sample, value in combined_by_sample.items() if value},
    }
    available_lineage_types = [
        lineage_type
        for lineage_type in ["mlst_ST", "ani_cluster", "bioproject", "combined_lineage_label"]
        if any(lineage_maps[lineage_type].get(sample, "") for sample in samples)
    ]
    default_lineage_type = "mlst_ST" if "mlst_ST" in available_lineage_types else ("ani_cluster" if "ani_cluster" in available_lineage_types else (available_lineage_types[0] if available_lineage_types else "bioproject"))

    summary_rows = []
    for sample in samples:
        mlst = mlst_by_sample.get(sample, "")
        ani = ani_by_sample.get(sample, "")
        bioproject = bioproject_by_sample.get(sample, "")
        combined = combined_by_sample.get(sample, "")
        status_flags = []
        if not mlst:
            status_flags.append("missing_MLST")
        if not ani:
            status_flags.append("missing_ANI")
        if not combined:
            status_flags.append("lineage_data_unavailable")
        summary_rows.append({
            "assembly_accession": sample,
            "sample_id": metadata_by_sample.get(sample, {}).get("sample_id", sample),
            "mlst_ST": mlst,
            "ani_cluster": ani,
            "bioproject": bioproject,
            "combined_lineage_label": combined,
            "lineage_data_status": "lineage_context_available" if combined else "lineage_data_unavailable",
            "warning_flags": _flag_string(status_flags),
        })
    summary_path = tables / "lineage_summary.tsv"
    summary_fields = ["assembly_accession", "sample_id", "mlst_ST", "ani_cluster", "bioproject", "combined_lineage_label", "lineage_data_status", "warning_flags"]
    write_rows(summary_path, summary_rows, summary_fields)

    def group_samples_by_lineage(lineage_type: str) -> dict[str, set[str]]:
        groups: dict[str, set[str]] = defaultdict(set)
        values = lineage_maps.get(lineage_type, {})
        for sample in samples:
            value = values.get(sample, "")
            if value and not is_missing_value(value):
                groups[value].add(sample)
        return dict(groups)

    def top_value(sample_set: set[str], column: str) -> tuple[str, int, float]:
        values = [metadata_by_sample.get(sample, {}).get(column, "") for sample in sample_set]
        values = [value for value in values if value and not is_missing_value(value)]
        if not values:
            return "", 0, 0.0
        value, count = Counter(values).most_common(1)[0]
        return value, count, count / len(values) if values else 0.0

    distribution_rows = []
    for lineage_type in available_lineage_types:
        groups = group_samples_by_lineage(lineage_type)
        for lineage_id, lineage_samples in sorted(groups.items(), key=lambda item: (-len(item[1]), item[0])):
            lineage_n = len(lineage_samples)
            top_country, _, top_country_fraction = top_value(lineage_samples, "country")
            top_source, _, top_source_fraction = top_value(lineage_samples, "isolation_source")
            top_project, _, top_project_fraction = top_value(lineage_samples, "bioproject")
            dataset_fraction = lineage_n / sample_count if sample_count else 0.0
            flags = []
            if lineage_n < 5:
                flags.append("small_lineage_group")
            if dataset_fraction >= 0.50:
                flags.append("single_lineage_dominance")
                if lineage_type == "mlst_ST":
                    flags.append("dominant_ST")
                if lineage_type == "ani_cluster":
                    flags.append("dominant_ANI_cluster")
            if top_project_fraction >= 0.70 and lineage_n >= 3:
                flags.append("BioProject_lineage_overlap")
            distribution_rows.append({
                "lineage_type": lineage_type,
                "lineage_id": lineage_id,
                "total_genomes": str(lineage_n),
                "fraction_of_dataset": f"{dataset_fraction:.4f}",
                "fraction_display": f"{dataset_fraction * 100:.1f}% ({lineage_n}/{sample_count})" if sample_count else "",
                "typed_or_clustered_status": "typed_or_clustered",
                "top_country": top_country,
                "top_country_fraction": f"{top_country_fraction:.4f}",
                "top_source": top_source,
                "top_source_fraction": f"{top_source_fraction:.4f}",
                "top_bioproject": top_project,
                "top_bioproject_fraction": f"{top_project_fraction:.4f}",
                "warning_flags": _flag_string(flags),
            })
    distribution_path = tables / "lineage_distribution.tsv"
    distribution_fields = [
        "lineage_type", "lineage_id", "total_genomes", "fraction_of_dataset", "fraction_display", "typed_or_clustered_status",
        "top_country", "top_country_fraction", "top_source", "top_source_fraction", "top_bioproject", "top_bioproject_fraction", "warning_flags",
    ]
    write_rows(distribution_path, distribution_rows, distribution_fields)

    metadata_columns = [
        column
        for column in ["country", "continent", "host", "isolation_source", "sample_type", "collection_year", "bioproject"]
        if any(clean_metadata_value(metadata_by_sample.get(sample, {}).get(column, "")) for sample in samples)
    ]
    default_metadata_column = "isolation_source" if "isolation_source" in metadata_columns else ("country" if "country" in metadata_columns else (metadata_columns[0] if metadata_columns else "bioproject"))

    overlap_rows = []
    for lineage_type in available_lineage_types:
        values = lineage_maps.get(lineage_type, {})
        for metadata_column in metadata_columns:
            groups: dict[str, set[str]] = defaultdict(set)
            for sample in samples:
                group_value = clean_metadata_value(metadata_by_sample.get(sample, {}).get(metadata_column, ""))
                if group_value:
                    groups[group_value].add(sample)
            for metadata_group, group_samples in sorted(groups.items(), key=lambda item: (-len(item[1]), item[0])):
                lineage_counts = Counter(values.get(sample, "") for sample in group_samples if values.get(sample, ""))
                if not lineage_counts:
                    continue
                dominant_lineage, dominant_count = lineage_counts.most_common(1)[0]
                dominant_fraction = dominant_count / len(group_samples) if group_samples else 0.0
                dominance_label = _lineage_dominance_label(dominant_fraction)
                for lineage_id, lineage_count in sorted(lineage_counts.items(), key=lambda item: (-item[1], item[0])):
                    fraction = lineage_count / len(group_samples) if group_samples else 0.0
                    flags = []
                    if dominance_label != "not_lineage_dominated":
                        flags.append("metadata_lineage_confounding")
                    if dominance_label in {"strong_lineage_confounding", "severe_lineage_confounding"}:
                        flags.append("single_lineage_dominance")
                    if len(group_samples) < 5:
                        flags.append("small_lineage_group")
                    overlap_rows.append({
                        "lineage_type": lineage_type,
                        "metadata_column": metadata_column,
                        "metadata_group": metadata_group,
                        "lineage_id": lineage_id,
                        "group_total_genomes": str(len(group_samples)),
                        "lineage_genomes_in_group": str(lineage_count),
                        "lineage_fraction_in_group": f"{fraction:.4f}",
                        "lineage_fraction_percent": f"{fraction * 100:.2f}",
                        "dominant_lineage": dominant_lineage,
                        "dominant_lineage_fraction": f"{dominant_fraction:.4f}",
                        "dominance_label": dominance_label,
                        "warning_flags": _flag_string(flags),
                    })
    overlap_path = tables / "lineage_metadata_overlap.tsv"
    overlap_fields = [
        "lineage_type", "metadata_column", "metadata_group", "lineage_id", "group_total_genomes", "lineage_genomes_in_group",
        "lineage_fraction_in_group", "lineage_fraction_percent", "dominant_lineage", "dominant_lineage_fraction", "dominance_label", "warning_flags",
    ]
    write_rows(overlap_path, overlap_rows, overlap_fields)

    presence = feature_presence(features)
    feature_name_by_key = {}
    category_by_key = {}
    feature_rows_by_sample_database: dict[tuple[str, str], int] = defaultdict(int)
    feature_rows_by_sample_database_feature: dict[tuple[str, str, str], int] = defaultdict(int)
    feature_samples_by_database: dict[tuple[str, str], set[str]] = defaultdict(set)
    for row in features:
        if row.get("presence", "1") != "1":
            continue
        sample = row.get("assembly_accession", "") or row.get("sample_id", "")
        database = row.get("database", "")
        feature_id = row.get("feature_id", "")
        if not sample or not database or not feature_id:
            continue
        feature_name_by_key[(database, feature_id)] = row.get("feature_name", "")
        category_by_key[(database, feature_id)] = row.get("feature_category", "")
        feature_rows_by_sample_database[(sample, database)] += 1
        feature_rows_by_sample_database_feature[(sample, database, feature_id)] += 1
        feature_samples_by_database[(sample, database)].add(feature_id)
    databases = sorted({database for database, _ in presence})
    default_database = "amr" if "amr" in databases else (databases[0] if databases else "")

    by_database: dict[str, list[tuple[tuple[str, str], set[str]]]] = defaultdict(list)
    for key, present_samples in presence.items():
        by_database[key[0]].append((key, present_samples))
    limited_presence: dict[tuple[str, str], set[str]] = {}
    capped_databases = []
    for database, items in by_database.items():
        ranked = sorted(items, key=lambda item: (-len(item[1]), item[0][1]))
        if len(ranked) > max_features_per_database:
            capped_databases.append(database)
        for key, present_samples in ranked[:max_features_per_database]:
            limited_presence[key] = present_samples

    burden_rows = []
    for lineage_type in available_lineage_types:
        groups = group_samples_by_lineage(lineage_type)
        for lineage_id, lineage_samples in groups.items():
            lineage_n = len(lineage_samples)
            for database in databases:
                counts = [feature_rows_by_sample_database.get((sample, database), 0) for sample in sorted(lineage_samples)]
                positive_samples = {sample for sample in lineage_samples if feature_samples_by_database.get((sample, database), set())}
                unique_features = set()
                for sample in lineage_samples:
                    unique_features.update(feature_samples_by_database.get((sample, database), set()))
                flags = []
                if lineage_n < 5:
                    flags.append("small_lineage_group")
                burden_rows.append({
                    "lineage_type": lineage_type,
                    "lineage_id": lineage_id,
                    "database": database,
                    "lineage_n": str(lineage_n),
                    "total_feature_rows": str(sum(counts)),
                    "unique_features": str(len(unique_features)),
                    "positive_genomes": str(len(positive_samples)),
                    "prevalence_percent": f"{(len(positive_samples) / lineage_n * 100) if lineage_n else 0.0:.2f}",
                    "mean_features_per_genome": f"{_mean([float(value) for value in counts]):.4f}",
                    "median_features_per_genome": f"{_median([float(value) for value in counts]):.4f}",
                    "max_features_per_genome": str(max(counts) if counts else 0),
                    "support_label": _lineage_support_label(lineage_n),
                    "warning_flags": _flag_string(flags),
                })
    burden_path = tables / "lineage_feature_burden.tsv"
    burden_fields = [
        "lineage_type", "lineage_id", "database", "lineage_n", "total_feature_rows", "unique_features", "positive_genomes",
        "prevalence_percent", "mean_features_per_genome", "median_features_per_genome", "max_features_per_genome", "support_label", "warning_flags",
    ]
    write_rows(burden_path, burden_rows, burden_fields)

    presence_rows = []
    for lineage_type in available_lineage_types:
        groups = group_samples_by_lineage(lineage_type)
        for lineage_id, lineage_samples in groups.items():
            lineage_n = len(lineage_samples)
            for (database, feature_id), present_samples in sorted(presence.items(), key=lambda item: (item[0][0], item[0][1])):
                positive_samples = present_samples & lineage_samples
                if not positive_samples:
                    continue
                feature_rows_count = sum(feature_rows_by_sample_database_feature.get((sample, database, feature_id), 0) for sample in positive_samples)
                top_country, _, _ = top_value(positive_samples, "country")
                top_source, _, _ = top_value(positive_samples, "isolation_source")
                top_project, _, top_project_fraction = top_value(positive_samples, "bioproject")
                flags = []
                if lineage_n < 5:
                    flags.append("small_lineage_group")
                if top_project_fraction >= 0.70 and len(positive_samples) >= 3:
                    flags.append("BioProject_lineage_overlap")
                prevalence_percent = (len(positive_samples) / lineage_n * 100) if lineage_n else 0.0
                presence_rows.append({
                    "database": database,
                    "feature_id": feature_id,
                    "feature_name": feature_name_by_key.get((database, feature_id), ""),
                    "lineage_type": lineage_type,
                    "lineage_id": lineage_id,
                    "lineage_n": str(lineage_n),
                    "positive_genomes": str(len(positive_samples)),
                    "prevalence_percent": f"{prevalence_percent:.2f}",
                    "prevalence_display": f"{prevalence_percent:.1f}% ({len(positive_samples)}/{lineage_n})" if lineage_n else "",
                    "feature_rows": str(feature_rows_count),
                    "mean_hits_per_positive_genome": f"{(feature_rows_count / len(positive_samples)) if positive_samples else 0.0:.4f}",
                    "top_country": top_country,
                    "top_source": top_source,
                    "top_bioproject": top_project,
                    "warning_flags": _flag_string(flags),
                })
    presence_path = tables / "lineage_feature_presence.tsv"
    presence_fields = [
        "database", "feature_id", "feature_name", "lineage_type", "lineage_id", "lineage_n", "positive_genomes",
        "prevalence_percent", "prevalence_display", "feature_rows", "mean_hits_per_positive_genome", "top_country", "top_source", "top_bioproject", "warning_flags",
    ]
    write_rows(presence_path, presence_rows, presence_fields)

    enrichment_rows = []
    for lineage_type in available_lineage_types:
        groups = group_samples_by_lineage(lineage_type)
        lineage_samples_any = set().union(*groups.values()) if groups else set()
        for lineage_id, lineage_samples in groups.items():
            outside_samples = lineage_samples_any - lineage_samples
            lineage_n = len(lineage_samples)
            outside_n = len(outside_samples)
            if lineage_n == 0 or outside_n == 0:
                continue
            for (database, feature_id), present_samples in sorted(limited_presence.items(), key=lambda item: (item[0][0], item[0][1])):
                positive_lineage = present_samples & lineage_samples
                positive_outside = present_samples & outside_samples
                a = len(positive_lineage)
                b = lineage_n - a
                c = len(positive_outside)
                d = outside_n - c
                if a + c < 1:
                    continue
                prevalence_in = a / lineage_n if lineage_n else 0.0
                prevalence_out = c / outside_n if outside_n else 0.0
                support_label = _lineage_support_label(lineage_n, outside_n)
                p_value = None
                test_name = "descriptive_only"
                if support_label not in {"insufficient_support", "descriptive_only"} and (a + c) >= 3:
                    test_name, p_value = _binary_test_for_counts(a, b, c, d)
                flags = ["multiple_testing", "exploratory_only"]
                if lineage_n < 5:
                    flags.append("small_lineage_group")
                if a + c < 3:
                    flags.append("low_positive_count")
                top_project, _, top_project_fraction = top_value(positive_lineage, "bioproject")
                if top_project_fraction >= 0.70 and len(positive_lineage) >= 3:
                    flags.append("BioProject_lineage_overlap")
                enrichment_rows.append({
                    "lineage_type": lineage_type,
                    "lineage_id": lineage_id,
                    "database": database,
                    "feature_id": feature_id,
                    "feature_name": feature_name_by_key.get((database, feature_id), ""),
                    "feature_category": category_by_key.get((database, feature_id), ""),
                    "lineage_n": str(lineage_n),
                    "outside_lineage_n": str(outside_n),
                    "positive_in_lineage": str(a),
                    "positive_outside_lineage": str(c),
                    "prevalence_in_lineage": f"{prevalence_in:.4f}",
                    "prevalence_outside_lineage": f"{prevalence_out:.4f}",
                    "prevalence_difference": f"{prevalence_in - prevalence_out:.4f}",
                    "prevalence_difference_percent": f"{(prevalence_in - prevalence_out) * 100:.2f}",
                    "odds_ratio": f"{_odds_ratio(a, b, c, d):.6g}",
                    "test_name": test_name,
                    "p_value": _format_pvalue(p_value),
                    "q_value": "",
                    "support_label": support_label,
                    "warning_flags": _flag_string(flags),
                    "interpretation_label": "",
                })
    add_bh_qvalues(enrichment_rows)
    for row in enrichment_rows:
        q_value = _float_or_none(row.get("q_value", ""))
        row["interpretation_label"] = _lineage_interpretation_label(
            q_value,
            int(row.get("lineage_n", "0") or 0),
            int(row.get("outside_lineage_n", "0") or 0),
            _float_or_none(row.get("prevalence_in_lineage", "")) or 0.0,
            _float_or_none(row.get("prevalence_outside_lineage", "")) or 0.0,
            row.get("warning_flags", ""),
        )
    enrichment_rows = sorted(
        enrichment_rows,
        key=lambda row: (
            row.get("interpretation_label", "") not in {"feature_lineage_specific", "feature_lineage_enriched"},
            -abs(_float_or_none(row.get("prevalence_difference", "")) or 0.0),
            row.get("lineage_type", ""),
            row.get("database", ""),
            row.get("feature_id", ""),
        ),
    )
    enrichment_path = tables / "lineage_feature_enrichment.tsv"
    enrichment_fields = [
        "lineage_type", "lineage_id", "database", "feature_id", "feature_name", "feature_category", "lineage_n", "outside_lineage_n",
        "positive_in_lineage", "positive_outside_lineage", "prevalence_in_lineage", "prevalence_outside_lineage", "prevalence_difference",
        "prevalence_difference_percent", "odds_ratio", "test_name", "p_value", "q_value", "support_label", "warning_flags", "interpretation_label",
    ]
    write_rows(enrichment_path, enrichment_rows, enrichment_fields)

    adjusted_source = read_table(out_dir / "metadata_feature_analysis" / "lineage_adjusted_warnings.tsv")
    adjusted_rows = []
    for idx, row in enumerate(adjusted_source, start=1):
        st_fraction = _float_or_none(row.get("dominant_ST_fraction", ""))
        ani_fraction = _float_or_none(row.get("dominant_ani_cluster_fraction", "") or row.get("dominant_ANI_cluster_fraction", ""))
        best_fraction = max([value for value in [st_fraction, ani_fraction] if value is not None] or [0.0])
        lineage_label = _lineage_dominance_label(best_fraction)
        flags = row.get("lineage_warning_flags", "")
        adjusted_rows.append({
            "finding_id": str(idx),
            "database": row.get("database", ""),
            "feature_id": row.get("feature_id", ""),
            "metadata_column": row.get("metadata_column", ""),
            "metadata_group": row.get("metadata_value", "") or row.get("metadata_group", ""),
            "original_interpretation_label": row.get("finding_type", ""),
            "supporting_samples": row.get("supporting_samples", ""),
            "dominant_ST": row.get("dominant_ST", ""),
            "dominant_ST_fraction": row.get("dominant_ST_fraction", ""),
            "dominant_ANI_cluster": row.get("dominant_ani_cluster", "") or row.get("dominant_ANI_cluster", ""),
            "dominant_ANI_cluster_fraction": row.get("dominant_ani_cluster_fraction", "") or row.get("dominant_ANI_cluster_fraction", ""),
            "dominant_BioProject": "",
            "dominant_BioProject_fraction": "",
            "lineage_warning_flags": flags,
            "lineage_adjusted_interpretation": lineage_label,
            "recommended_interpretation": "treat_as_lineage_caution" if lineage_label != "not_lineage_dominated" else "lineage_not_dominant_in_supporting_samples",
        })
    adjusted_path = tables / "lineage_adjusted_top_findings.tsv"
    adjusted_fields = [
        "finding_id", "database", "feature_id", "metadata_column", "metadata_group", "original_interpretation_label", "supporting_samples",
        "dominant_ST", "dominant_ST_fraction", "dominant_ANI_cluster", "dominant_ANI_cluster_fraction", "dominant_BioProject", "dominant_BioProject_fraction",
        "lineage_warning_flags", "lineage_adjusted_interpretation", "recommended_interpretation",
    ]
    write_rows(adjusted_path, adjusted_rows, adjusted_fields)

    figure_paths: list[Path] = []

    def write_bar_figure(base_name: str, title: str, rows: list[dict[str, str]], label_field: str, value_field: str, fields: list[str]) -> tuple[Path, Path, Path, Path]:
        data_path = figures / f"{base_name}.data.tsv"
        svg_path = figures / f"{base_name}.svg"
        png_path = figures / f"{base_name}.png"
        pdf_path = figures / f"{base_name}.pdf"
        write_rows(data_path, rows, fields)
        _write_bar_svg(svg_path, rows, title, label_field, value_field)
        _write_bar_png(png_path, rows, value_field)
        _write_simple_pdf(pdf_path, title, [f"{row.get(label_field, '')}: {row.get(value_field, '')}" for row in rows[:60]])
        figure_paths.extend([data_path, svg_path, png_path, pdf_path])
        return data_path, svg_path, png_path, pdf_path

    default_distribution_svg = figures / f"lineage_distribution_{_safe_filename(default_lineage_type)}.svg"
    for lineage_type in available_lineage_types or [default_lineage_type]:
        rows = [
            row for row in distribution_rows
            if row.get("lineage_type") == lineage_type
        ][:top_n]
        write_bar_figure(
            f"lineage_distribution_{_safe_filename(lineage_type)}",
            f"Lineage Distribution: {lineage_type}",
            rows,
            "lineage_id",
            "total_genomes",
            distribution_fields,
        )

    default_overlap_svg = figures / f"lineage_metadata_overlap_{_safe_filename(default_lineage_type)}_{_safe_filename(default_metadata_column)}.svg"
    for lineage_type in available_lineage_types or [default_lineage_type]:
        rows = [
            row for row in overlap_rows
            if row.get("lineage_type") == lineage_type and row.get("metadata_column") == default_metadata_column
        ]
        top_groups = {row.get("metadata_group", "") for row in sorted(rows, key=lambda row: (-int(row.get("group_total_genomes", "0") or 0), row.get("metadata_group", "")))[:10]}
        top_lineages = {row.get("lineage_id", "") for row in sorted(rows, key=lambda row: (-int(row.get("lineage_genomes_in_group", "0") or 0), row.get("lineage_id", "")))[:10]}
        plot_rows = [row for row in rows if row.get("metadata_group") in top_groups and row.get("lineage_id") in top_lineages]
        base = f"lineage_metadata_overlap_{_safe_filename(lineage_type)}_{_safe_filename(default_metadata_column)}"
        data_path = figures / f"{base}.data.tsv"
        svg_path = figures / f"{base}.svg"
        png_path = figures / f"{base}.png"
        pdf_path = figures / f"{base}.pdf"
        write_rows(data_path, plot_rows, overlap_fields)
        _write_heatmap_svg(svg_path, plot_rows, f"{lineage_type} vs {default_metadata_column}", "metadata_group", "lineage_id", "lineage_fraction_percent")
        _write_heatmap_png(png_path, plot_rows, "metadata_group", "lineage_id", "lineage_fraction_percent")
        _write_simple_pdf(pdf_path, f"{lineage_type} vs {default_metadata_column}", [f"{row.get('metadata_group', '')}/{row.get('lineage_id', '')}: {row.get('lineage_fraction_percent', '')}%" for row in plot_rows[:60]])
        figure_paths.extend([data_path, svg_path, png_path, pdf_path])

    default_burden_svg = figures / f"lineage_database_burden_{_safe_filename(default_database)}_{_safe_filename(default_lineage_type)}.svg"
    for lineage_type in available_lineage_types or [default_lineage_type]:
        groups = group_samples_by_lineage(lineage_type)
        ordered_groups = sorted(groups.items(), key=lambda item: (-len(item[1]), item[0]))[:top_n]
        boxplot_rows = []
        for lineage_id, lineage_samples in ordered_groups:
            values = [float(feature_rows_by_sample_database.get((sample, default_database), 0)) for sample in sorted(lineage_samples)]
            stats = _summary_stats_full(values)
            boxplot_rows.append({
                "metadata_group": lineage_id,
                "lineage_type": lineage_type,
                "lineage_id": lineage_id,
                "database": default_database,
                "n": str(len(values)),
                "min": "" if stats["min"] is None else f"{stats['min']:.4f}",
                "q1": "" if stats["q1"] is None else f"{stats['q1']:.4f}",
                "median": "" if stats["median"] is None else f"{stats['median']:.4f}",
                "q3": "" if stats["q3"] is None else f"{stats['q3']:.4f}",
                "max": "" if stats["max"] is None else f"{stats['max']:.4f}",
                "mean": "" if stats["mean"] is None else f"{stats['mean']:.4f}",
            })
        base = f"lineage_database_burden_{_safe_filename(default_database)}_{_safe_filename(lineage_type)}"
        data_path = figures / f"{base}.data.tsv"
        svg_path = figures / f"{base}.svg"
        png_path = figures / f"{base}.png"
        pdf_path = figures / f"{base}.pdf"
        boxplot_fields = ["metadata_group", "lineage_type", "lineage_id", "database", "n", "min", "q1", "median", "q3", "max", "mean"]
        write_rows(data_path, boxplot_rows, boxplot_fields)
        _write_burden_boxplot_svg(svg_path, boxplot_rows, f"{default_database} Burden by {lineage_type}")
        _write_burden_boxplot_png(png_path, boxplot_rows)
        _write_simple_pdf(
            pdf_path,
            f"{default_database} Burden by {lineage_type}",
            [f"{row.get('lineage_id', '')}: median={row.get('median', '')}, n={row.get('n', '')}" for row in boxplot_rows],
        )
        figure_paths.extend([data_path, svg_path, png_path, pdf_path])

    heatmap_feature_rows = [
        row for row in presence_rows
        if row.get("database") == default_database and row.get("lineage_type") == default_lineage_type
    ]
    top_heatmap_features = {row.get("feature_id", "") for row in sorted(heatmap_feature_rows, key=lambda row: (-(_float_or_none(row.get("prevalence_percent", "")) or 0.0), row.get("feature_id", "")))[:top_n]}
    top_heatmap_lineages = {row.get("lineage_id", "") for row in sorted(heatmap_feature_rows, key=lambda row: (-int(row.get("lineage_n", "0") or 0), row.get("lineage_id", "")))[:top_n]}
    heatmap_rows = [row for row in heatmap_feature_rows if row.get("feature_id") in top_heatmap_features and row.get("lineage_id") in top_heatmap_lineages]
    heatmap_base = f"lineage_feature_heatmap_{_safe_filename(default_database)}_{_safe_filename(default_lineage_type)}"
    heatmap_data = figures / f"{heatmap_base}.data.tsv"
    heatmap_svg = figures / f"{heatmap_base}.svg"
    heatmap_png = figures / f"{heatmap_base}.png"
    heatmap_pdf = figures / f"{heatmap_base}.pdf"
    write_rows(heatmap_data, heatmap_rows, presence_fields)
    _write_heatmap_svg(heatmap_svg, heatmap_rows, f"{default_database} Feature Prevalence by {default_lineage_type}", "feature_id", "lineage_id", "prevalence_percent")
    _write_heatmap_png(heatmap_png, heatmap_rows, "feature_id", "lineage_id", "prevalence_percent")
    _write_simple_pdf(heatmap_pdf, f"{default_database} Feature Prevalence by {default_lineage_type}", [f"{row.get('feature_id', '')}/{row.get('lineage_id', '')}: {row.get('prevalence_display', '')}" for row in heatmap_rows[:60]])
    figure_paths.extend([heatmap_data, heatmap_svg, heatmap_png, heatmap_pdf])

    enrichment_plot_rows = [
        {
            **row,
            "feature_lineage": f"{row.get('feature_id', '')} / {row.get('lineage_id', '')}",
            "abs_prevalence_difference_percent": f"{abs(_float_or_none(row.get('prevalence_difference_percent', '')) or 0.0):.2f}",
        }
        for row in enrichment_rows
        if row.get("database") == default_database and row.get("lineage_type") == default_lineage_type
    ][:top_n]
    enrichment_base = f"lineage_feature_enrichment_{_safe_filename(default_database)}_{_safe_filename(default_lineage_type)}"
    enrichment_data = figures / f"{enrichment_base}.data.tsv"
    enrichment_svg = figures / f"{enrichment_base}.svg"
    enrichment_png = figures / f"{enrichment_base}.png"
    enrichment_pdf = figures / f"{enrichment_base}.pdf"
    enrichment_plot_fields = [*enrichment_fields, "feature_lineage", "abs_prevalence_difference_percent"]
    write_rows(enrichment_data, enrichment_plot_rows, enrichment_plot_fields)
    _write_metadata_volcano_svg(enrichment_svg, enrichment_plot_rows, f"{default_database} Lineage Enrichment")
    _write_metadata_volcano_png(enrichment_png, enrichment_plot_rows)
    _write_simple_pdf(enrichment_pdf, f"{default_database} Lineage Enrichment", [f"{row.get('feature_lineage', '')}: {row.get('prevalence_difference_percent', '')} pp" for row in enrichment_plot_rows])
    figure_paths.extend([enrichment_data, enrichment_svg, enrichment_png, enrichment_pdf])

    default_feature = ""
    for row in sorted(presence_rows, key=lambda row: (row.get("database", "") != default_database, -(_float_or_none(row.get("positive_genomes", "")) or 0.0), row.get("feature_id", ""))):
        if row.get("database") == default_database:
            default_feature = row.get("feature_id", "")
            break
    selected_presence_rows = [
        row for row in presence_rows
        if row.get("database") == default_database and row.get("feature_id") == default_feature and row.get("lineage_type") == default_lineage_type
    ]
    selected_presence_rows = sorted(selected_presence_rows, key=lambda row: (-(_float_or_none(row.get("prevalence_percent", "")) or 0.0), row.get("lineage_id", "")))[:top_n]
    selected_base = f"lineage_feature_presence_{_safe_filename(default_database)}_{_safe_filename(default_feature or 'feature')}_{_safe_filename(default_lineage_type)}"
    selected_data = figures / f"{selected_base}.data.tsv"
    selected_svg = figures / f"{selected_base}.svg"
    selected_png = figures / f"{selected_base}.png"
    selected_pdf = figures / f"{selected_base}.pdf"
    write_rows(selected_data, selected_presence_rows, presence_fields)
    _write_bar_svg(selected_svg, selected_presence_rows, f"{default_feature or default_database} by {default_lineage_type}", "lineage_id", "prevalence_percent", "Prevalence (%)")
    _write_bar_png(selected_png, selected_presence_rows, "prevalence_percent")
    _write_simple_pdf(selected_pdf, f"{default_feature or default_database} by {default_lineage_type}", [f"{row.get('lineage_id', '')}: {row.get('prevalence_display', '')}" for row in selected_presence_rows])
    figure_paths.extend([selected_data, selected_svg, selected_png, selected_pdf])

    confounding_plot_rows = []
    severity_score = {
        "not_lineage_dominated": 0,
        "possible_lineage_confounding": 1,
        "strong_lineage_confounding": 2,
        "severe_lineage_confounding": 3,
        "insufficient_lineage_data": 1,
    }
    for row in adjusted_rows[:top_n]:
        label = row.get("lineage_adjusted_interpretation", "")
        confounding_plot_rows.append({
            **row,
            "finding_label": f"{row.get('feature_id', '')} / {row.get('metadata_group', '')}"[:80],
            "lineage_warning_score": str(severity_score.get(label, 0)),
        })
    confounding_base = "lineage_confounding_top_findings"
    confounding_data = figures / f"{confounding_base}.data.tsv"
    confounding_svg = figures / f"{confounding_base}.svg"
    confounding_png = figures / f"{confounding_base}.png"
    confounding_pdf = figures / f"{confounding_base}.pdf"
    write_rows(confounding_data, confounding_plot_rows, [*adjusted_fields, "finding_label", "lineage_warning_score"])
    _write_bar_svg(confounding_svg, confounding_plot_rows, "Lineage Confounding In Top Findings", "finding_label", "lineage_warning_score", "Warning severity")
    _write_bar_png(confounding_png, confounding_plot_rows, "lineage_warning_score")
    _write_simple_pdf(confounding_pdf, "Lineage Confounding In Top Findings", [f"{row.get('finding_label', '')}: {row.get('lineage_adjusted_interpretation', '')}" for row in confounding_plot_rows])
    figure_paths.extend([confounding_data, confounding_svg, confounding_png, confounding_pdf])

    lineage_warning_count = sum(1 for row in distribution_rows + overlap_rows + burden_rows + presence_rows + enrichment_rows + adjusted_rows if row.get("warning_flags") or row.get("lineage_warning_flags"))
    dominant_st = next((row for row in distribution_rows if row.get("lineage_type") == "mlst_ST"), {})
    dominant_ani = next((row for row in distribution_rows if row.get("lineage_type") == "ani_cluster"), {})
    dominant_project = next((row for row in distribution_rows if row.get("lineage_type") == "bioproject"), {})
    top_overlap = next(iter(sorted(overlap_rows, key=lambda row: (-(_float_or_none(row.get("dominant_lineage_fraction", "")) or 0.0), row.get("metadata_column", ""), row.get("metadata_group", "")))), {})
    top_burden = next(iter(sorted(burden_rows, key=lambda row: (-(_float_or_none(row.get("median_features_per_genome", "")) or 0.0), row.get("database", ""), row.get("lineage_id", "")))), {})
    top_selected_presence = next(iter(sorted(selected_presence_rows, key=lambda row: (-(_float_or_none(row.get("prevalence_percent", "")) or 0.0), row.get("lineage_id", "")))), {})
    possible_confounding = sum(1 for row in adjusted_rows if row.get("lineage_adjusted_interpretation") == "possible_lineage_confounding")
    strong_confounding = sum(1 for row in adjusted_rows if row.get("lineage_adjusted_interpretation") == "strong_lineage_confounding")
    severe_confounding = sum(1 for row in adjusted_rows if row.get("lineage_adjusted_interpretation") == "severe_lineage_confounding")
    mlst_genome_count = sum(1 for sample in samples if mlst_by_sample.get(sample))
    ani_genome_count = sum(1 for sample in samples if ani_by_sample.get(sample))
    if not mlst_genome_count and not ani_genome_count:
        availability_note = "MLST and ANI cluster information were not available. Lineage-aware interpretation could not be performed; BioProject context is still reported where possible."
    elif not ani_genome_count:
        availability_note = "ANI cluster information was not available. Lineage interpretation is based on MLST and BioProject context where available."
    elif not mlst_genome_count:
        availability_note = "MLST sequence type information was not available. Lineage interpretation is based on ANI clusters and BioProject context where available."
    else:
        availability_note = "MLST and ANI cluster context were both available for lineage-aware interpretation."
    written_rows = [
        {
            "section": "overall_lineage_summary",
            "summary": (
                f"Lineage analysis identified {len(set(mlst_by_sample.values()))} MLST sequence type(s) among {mlst_genome_count} typed genome(s) "
                f"and {len(set(ani_by_sample.values()))} ANI cluster(s) among {ani_genome_count} genome(s). "
                f"The largest MLST group was {dominant_st.get('lineage_id', 'not available')} ({dominant_st.get('fraction_display', 'not available')}); "
                f"the largest ANI cluster was {dominant_ani.get('lineage_id', 'not available')} ({dominant_ani.get('fraction_display', 'not available')}); "
                f"the largest BioProject group was {dominant_project.get('lineage_id', 'not available')} ({dominant_project.get('fraction_display', 'not available')}). "
                f"{availability_note}"
            ),
        },
        {
            "section": "metadata_lineage_overlap_summary",
            "summary": (
                (
                    f"Metadata-lineage overlap screening found {len(overlap_rows)} group-by-lineage rows. "
                    f"The strongest default imbalance was {top_overlap.get('metadata_column', 'metadata')}={top_overlap.get('metadata_group', 'not available')} "
                    f"dominated by {top_overlap.get('dominant_lineage', 'not available')} ({top_overlap.get('dominant_lineage_fraction', '')}). "
                    "Dominated metadata groups should be treated as possible clonal-structure confounding."
                )
                if overlap_rows else "Metadata-lineage overlap could not be summarized because no lineage group and metadata group overlap was available."
            ),
        },
        {
            "section": "feature_burden_by_lineage_summary",
            "summary": (
                (
                    f"Feature burden by lineage was summarized for {len(burden_rows)} database-lineage combinations. "
                    f"The highest median report-facing burden row was {top_burden.get('database', 'database')} in {top_burden.get('lineage_id', 'not available')} "
                    f"with median {top_burden.get('median_features_per_genome', '0')} feature row(s) per genome."
                )
                if burden_rows else "Feature burden by lineage could not be summarized because no lineage-feature overlap was available."
            ),
        },
        {
            "section": "selected_feature_lineage_summary",
            "summary": (
                (
                    f"The default selected feature view is {default_database}:{default_feature} by {default_lineage_type}. "
                    f"The top lineage carrying it was {top_selected_presence.get('lineage_id', 'not available')} "
                    f"({top_selected_presence.get('prevalence_display', 'not available')})."
                )
                if default_feature else "No default selected-feature lineage prevalence view was available because no feature-lineage overlap was detected."
            ),
        },
        {
            "section": "lineage_adjusted_top_findings_summary",
            "summary": (
                f"Among lineage-adjusted top findings, {possible_confounding} were flagged for possible lineage confounding, "
                f"{strong_confounding} for strong lineage confounding, and {severe_confounding} for severe lineage confounding. "
                "These findings remain useful for exploration but should be checked with lineage-aware or phylogeny-aware analysis before strong interpretation."
            ),
        },
    ]
    written_path = tables / "lineage_written_summaries.tsv"
    write_rows(written_path, written_rows, ["section", "summary"])
    summary_stat_rows = [
        {"metric": "samples", "value": str(sample_count), "message": "Samples represented in the important lineage report."},
        {"metric": "lineage_types_available", "value": ",".join(available_lineage_types) or "none", "message": "Lineage/context views available for this run."},
        {"metric": "genomes_with_mlst_ST", "value": str(sum(1 for sample in samples if mlst_by_sample.get(sample))), "message": "Genomes with MLST sequence type context."},
        {"metric": "unique_STs", "value": str(len(set(mlst_by_sample.values()))), "message": "Unique MLST sequence types detected."},
        {"metric": "dominant_ST", "value": dominant_st.get("lineage_id", ""), "message": dominant_st.get("fraction_display", "")},
        {"metric": "genomes_with_ANI_cluster", "value": str(sum(1 for sample in samples if ani_by_sample.get(sample))), "message": "Genomes with ANI/skani cluster context."},
        {"metric": "unique_ANI_clusters", "value": str(len(set(ani_by_sample.values()))), "message": "Unique ANI clusters detected."},
        {"metric": "dominant_ANI_cluster", "value": dominant_ani.get("lineage_id", ""), "message": dominant_ani.get("fraction_display", "")},
        {"metric": "bioprojects_detected", "value": str(len(set(bioproject_by_sample.values()))), "message": "BioProject identifiers represented in metadata."},
        {"metric": "lineage_distribution_rows", "value": str(len(distribution_rows)), "message": "Lineage distribution rows."},
        {"metric": "metadata_lineage_overlap_rows", "value": str(len(overlap_rows)), "message": "Metadata-lineage overlap rows."},
        {"metric": "lineage_feature_burden_rows", "value": str(len(burden_rows)), "message": "Database-burden by lineage rows."},
        {"metric": "lineage_feature_enrichment_rows", "value": str(len(enrichment_rows)), "message": "Feature-lineage enrichment comparisons."},
        {"metric": "lineage_adjusted_top_findings_rows", "value": str(len(adjusted_rows)), "message": "Top findings with lineage-aware caution labels."},
        {"metric": "lineage_written_summaries", "value": str(len(written_rows)), "message": "Auto-written lineage interpretation summaries."},
        {"metric": "warning_rows", "value": str(lineage_warning_count), "message": "Rows carrying lineage/context warning flags."},
        {"metric": "feature_cap", "value": str(max_features_per_database), "message": f"Lineage enrichment feature cap per database. Capped databases: {','.join(sorted(capped_databases)) or 'none'}."},
        {"metric": "default_view", "value": f"{default_lineage_type}|{default_database}|{default_metadata_column}|{default_feature}", "message": "Default lineage report controls."},
    ]
    summary_stats_path = tables / "lineage_report_summary.tsv"
    write_rows(summary_stats_path, summary_stat_rows, ["metric", "value", "message"])

    table_paths = [summary_path, distribution_path, overlap_path, burden_path, enrichment_path, adjusted_path, presence_path, summary_stats_path, written_path]
    tables_zip = important_dir / "lineage_tables.zip"
    figures_zip = important_dir / "lineage_figures.zip"
    _write_zip_bundle(tables_zip, table_paths, important_dir)
    _write_zip_bundle(figures_zip, figure_paths, important_dir)

    written_html = "".join(
        f"<p><strong>{html.escape(row.get('section', '').replace('_', ' ').title())}:</strong> {html.escape(row.get('summary', ''))}</p>"
        for row in written_rows
    )
    interactive_html = figures / "lineage_clonal_structure.html"
    interactive_html.write_text(
        f"""<!doctype html>
<html><head><meta charset="utf-8"><title>Lineage / Clonal Structure</title>
<style>
body {{ font-family: Arial, sans-serif; color: #1f2933; margin: 1.5rem; }}
label {{ font-weight: 700; margin-right: 0.35rem; }}
select, input {{ margin: 0 1rem 0.75rem 0; padding: 0.35rem; }}
.warning {{ background: #fff7ed; border-left: 4px solid #c2410c; padding: 0.75rem; margin: 0.75rem 0; }}
.panel {{ border: 1px solid #d9e2ec; background: #f8fafc; padding: 0.75rem; margin: 0.75rem 0; }}
.figure-row {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 1rem; }}
.figure-row img {{ max-width: 100%; border: 1px solid #d9e2ec; background: white; }}
table {{ border-collapse: collapse; width: 100%; font-size: 0.88rem; }}
th, td {{ border: 1px solid #d9e2ec; padding: 0.35rem; text-align: left; }}
th {{ background: #f0f4f8; }}
.scrollbox {{ max-height: 620px; overflow: auto; border: 1px solid #d9e2ec; }}
</style></head><body>
<h1>Lineage / Clonal Structure</h1>
<div class="warning">Lineage summaries are exploratory and do not replace phylogenetic analysis. Apparent metadata associations may reflect clonal structure, BioProject sampling, geography, or temporal sampling.</div>
<label for="lineageType">Lineage type</label><select id="lineageType"></select>
<label for="database">Database</label><select id="database"></select>
<label for="featureMode">Feature mode</label><select id="featureMode"><option value="burden" selected>Database burden</option><option value="feature">Individual feature</option><option value="enrichment">Top lineage-enriched features</option></select>
<label for="featureSearch">Search feature</label><input id="featureSearch" type="search" placeholder="type feature name">
<label for="feature">Feature</label><select id="feature"></select>
<label for="metadataOverlay">Metadata overlay</label><select id="metadataOverlay"></select>
<label for="minLineage">Minimum lineage size</label><select id="minLineage"><option value="0">All</option><option value="3" selected>n >= 3</option><option value="5">n >= 5</option><option value="10">n >= 10</option><option value="custom">custom</option></select>
<input id="customMinLineage" type="number" min="0" value="3" style="width:5rem" aria-label="Custom minimum lineage size">
<label for="displayCount">Display</label><select id="displayCount"><option value="10">Top 10</option><option value="20" selected>Top 20</option><option value="50">Top 50</option><option value="99999">Complete</option></select>
<div id="summary"></div>
<div class="panel"><h2>Written Summaries</h2>{written_html}</div>
<div class="panel"><h2>Default Visuals</h2><div class="figure-row">
<div><h3>Lineage distribution</h3><img src="{html.escape(default_distribution_svg.name)}" alt="Lineage distribution"></div>
<div><h3>Metadata-lineage overlap</h3><img src="{html.escape(default_overlap_svg.name)}" alt="Metadata-lineage overlap"></div>
<div><h3>Database burden by lineage</h3><img src="{html.escape(default_burden_svg.name)}" alt="Database burden by lineage"></div>
<div><h3>Feature prevalence heatmap</h3><img src="{html.escape(heatmap_svg.name)}" alt="Feature prevalence by lineage"></div>
<div><h3>Lineage enrichment</h3><img src="{html.escape(enrichment_svg.name)}" alt="Lineage enrichment"></div>
<div><h3>Selected feature prevalence</h3><img src="{html.escape(selected_svg.name)}" alt="Selected feature lineage prevalence"></div>
<div><h3>Lineage-adjusted top findings</h3><img src="{html.escape(confounding_svg.name)}" alt="Lineage-adjusted top findings"></div>
</div></div>
<h2>Filtered Results</h2><div id="table" class="scrollbox"></div>
<p><a href="../tables/lineage_summary.tsv">Download lineage summary</a> | <a href="../tables/lineage_distribution.tsv">Download distribution</a> | <a href="../tables/lineage_metadata_overlap.tsv">Download metadata overlap</a> | <a href="../tables/lineage_feature_burden.tsv">Download burden table</a> | <a href="../tables/lineage_feature_enrichment.tsv">Download enrichment table</a> | <a href="../tables/lineage_feature_presence.tsv">Download selected-feature lineage table</a> | <a href="../tables/lineage_adjusted_top_findings.tsv">Download lineage-adjusted top findings</a> | <a href="../tables/lineage_written_summaries.tsv">Download written summaries</a> | <a href="../lineage_tables.zip">Download lineage tables ZIP</a> | <a href="../lineage_figures.zip">Download lineage figures ZIP</a></p>
<script>
const distributionRows = {json.dumps(distribution_rows[:3000])};
const overlapRows = {json.dumps(overlap_rows[:5000])};
const burdenRows = {json.dumps(burden_rows[:5000])};
const enrichmentRows = {json.dumps(enrichment_rows[:5000])};
const presenceRows = {json.dumps(presence_rows[:5000])};
const adjustedRows = {json.dumps(adjusted_rows[:2000])};
const defaultLineageType = {json.dumps(default_lineage_type)};
const defaultDatabase = {json.dumps(default_database)};
const defaultFeature = {json.dumps(default_feature)};
const defaultMetadataColumn = {json.dumps(default_metadata_column)};
function num(value) {{ const n = Number(value); return Number.isFinite(n) ? n : 0; }}
const lineageType = document.getElementById('lineageType'), database = document.getElementById('database'), featureMode = document.getElementById('featureMode'), featureSearch = document.getElementById('featureSearch'), feature = document.getElementById('feature'), metadataOverlay = document.getElementById('metadataOverlay'), minLineage = document.getElementById('minLineage'), customMinLineage = document.getElementById('customMinLineage'), displayCount = document.getElementById('displayCount');
function fillSelect(select, values, preferred) {{ const unique = Array.from(new Set(values.filter(Boolean))).sort(); select.innerHTML = unique.map(v => `<option value="${{v}}">${{v}}</option>`).join(''); if (unique.includes(preferred)) select.value = preferred; }}
function refreshControls() {{
  fillSelect(lineageType, distributionRows.map(r => r.lineage_type), defaultLineageType);
  fillSelect(database, burdenRows.map(r => r.database), defaultDatabase);
  fillSelect(metadataOverlay, overlapRows.map(r => r.metadata_column), defaultMetadataColumn);
  refreshFeatures();
}}
function refreshFeatures() {{
  const query = featureSearch.value.toLowerCase();
  fillSelect(feature, presenceRows.filter(r => r.database === database.value && (!query || String(r.feature_id || '').toLowerCase().includes(query))).map(r => r.feature_id), defaultFeature);
}}
function selectedMinLineage() {{ return minLineage.value === 'custom' ? num(customMinLineage.value) : num(minLineage.value); }}
function currentRows() {{
  if (featureMode.value === 'feature') return presenceRows.filter(r => r.feature_id === feature.value);
  if (featureMode.value === 'enrichment') return enrichmentRows;
  return burdenRows;
}}
function passes(row) {{
  if (row.lineage_type !== lineageType.value) return false;
  if (row.database && row.database !== database.value) return false;
  if (num(row.lineage_n || row.total_genomes) < selectedMinLineage()) return false;
  return true;
}}
function render() {{
  const rows = currentRows().filter(passes).slice(0, Number(displayCount.value));
  const cols = featureMode.value === 'feature'
    ? ['database','feature_id','lineage_type','lineage_id','lineage_n','positive_genomes','prevalence_display','feature_rows','top_country','top_source','top_bioproject','warning_flags']
    : (featureMode.value === 'enrichment'
      ? ['lineage_type','lineage_id','database','feature_id','lineage_n','positive_in_lineage','prevalence_in_lineage','prevalence_difference','odds_ratio','q_value','support_label','interpretation_label','warning_flags']
      : ['lineage_type','lineage_id','database','lineage_n','total_feature_rows','unique_features','positive_genomes','prevalence_percent','median_features_per_genome','support_label','warning_flags']);
  const overlap = overlapRows.filter(r => r.lineage_type === lineageType.value && r.metadata_column === metadataOverlay.value).slice(0, 10);
  document.getElementById('summary').innerHTML = `<p>${{rows.length}} rows shown. Metadata overlay preview: ${{overlap.map(r => `${{r.metadata_group}}/${{r.lineage_id}}=${{r.lineage_fraction_percent}}%`).join('; ') || 'none available'}}.</p>`;
  document.getElementById('table').innerHTML = '<table><thead><tr>' + cols.map(c => `<th>${{c}}</th>`).join('') + '</tr></thead><tbody>' + rows.map(r => '<tr>' + cols.map(c => `<td>${{r[c] || ''}}</td>`).join('') + '</tr>').join('') + '</tbody></table>';
}}
featureMode.addEventListener('change', render);
database.addEventListener('change', () => {{ refreshFeatures(); render(); }});
featureSearch.addEventListener('input', () => {{ refreshFeatures(); render(); }});
[lineageType, feature, metadataOverlay, minLineage, customMinLineage, displayCount].forEach(el => el.addEventListener('change', render));
refreshControls(); render();
</script></body></html>
""",
        encoding="utf-8",
    )

    return {
        "important_lineage_summary": str(summary_path),
        "important_lineage_distribution": str(distribution_path),
        "important_lineage_metadata_overlap": str(overlap_path),
        "important_lineage_feature_burden": str(burden_path),
        "important_lineage_feature_enrichment": str(enrichment_path),
        "important_lineage_adjusted_top_findings": str(adjusted_path),
        "important_lineage_feature_presence": str(presence_path),
        "important_lineage_report_summary": str(summary_stats_path),
        "important_lineage_written_summaries": str(written_path),
        "important_lineage_distribution_svg": str(default_distribution_svg),
        "important_lineage_metadata_overlap_svg": str(default_overlap_svg),
        "important_lineage_database_burden_svg": str(default_burden_svg),
        "important_lineage_feature_heatmap_svg": str(heatmap_svg),
        "important_lineage_feature_enrichment_svg": str(enrichment_svg),
        "important_lineage_feature_presence_svg": str(selected_svg),
        "important_lineage_confounding_top_findings_svg": str(confounding_svg),
        "important_lineage_html": str(interactive_html),
        "important_lineage_tables_zip": str(tables_zip),
        "important_lineage_figures_zip": str(figures_zip),
    }


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
    metadata_association_outputs: dict[str, str],
    lineage_outputs: dict[str, str],
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
    prevalence_rows = read_table(important_dir / "tables" / "feature_prevalence.tsv") or read_table(important_dir / "key_tables" / "feature_prevalence_summary.tsv")
    prevalence_database_rows = read_table(important_dir / "tables" / "prevalence_summary_by_database.tsv")
    prevalence_core_rows = read_table(important_dir / "tables" / "prevalence_core_accessory_rare_summary.tsv")
    prevalence_written_rows = read_table(important_dir / "tables" / "prevalence_written_summaries.tsv")
    geographic_summary_rows = read_table(important_dir / "tables" / "geographic_distribution_summary.tsv")
    geographic_feature_rows = read_table(important_dir / "tables" / "geographic_feature_distribution.tsv")
    geographic_burden_rows = read_table(important_dir / "tables" / "geographic_database_burden.tsv")
    geographic_warning_rows = read_table(important_dir / "tables" / "geographic_warning_summary.tsv")
    variation_rows = read_table(important_dir / "key_tables" / "feature_variation_summary.tsv")
    variation_database_rows = read_table(important_dir / "key_tables" / "feature_variation_database_summary.tsv")
    temporal_rows = read_table(important_dir / "key_tables" / "temporal_trend_summary.tsv")
    cooccurrence_rows = read_table(important_dir / "tables" / "cooccurrence_pair_summary.tsv")
    cooccurrence_summary_rows = read_table(important_dir / "tables" / "cooccurrence_context_summary.tsv")
    context_rows = read_table(important_dir / "tables" / "genomic_context_evidence.tsv")
    metadata_feature_rows = read_table(important_dir / "tables" / "metadata_feature_enrichment.tsv")
    metadata_burden_rows = read_table(important_dir / "tables" / "metadata_burden_associations.tsv")
    metadata_category_rows = read_table(important_dir / "tables" / "metadata_category_enrichment.tsv")
    metadata_burden_omnibus_rows = read_table(important_dir / "tables" / "metadata_burden_omnibus.tsv")
    metadata_category_omnibus_rows = read_table(important_dir / "tables" / "metadata_category_omnibus.tsv")
    metadata_summary_rows = read_table(important_dir / "tables" / "metadata_association_summary.tsv")
    metadata_usability_rows = read_table(important_dir / "tables" / "metadata_usability_summary.tsv")
    lineage_summary_rows = read_table(important_dir / "tables" / "lineage_report_summary.tsv")
    lineage_distribution_rows = read_table(important_dir / "tables" / "lineage_distribution.tsv")
    lineage_overlap_rows = read_table(important_dir / "tables" / "lineage_metadata_overlap.tsv")
    lineage_burden_rows = read_table(important_dir / "tables" / "lineage_feature_burden.tsv")
    lineage_enrichment_rows = read_table(important_dir / "tables" / "lineage_feature_enrichment.tsv")
    lineage_presence_rows = read_table(important_dir / "tables" / "lineage_feature_presence.tsv")
    lineage_adjusted_rows = read_table(important_dir / "tables" / "lineage_adjusted_top_findings.tsv")
    lineage_written_rows = read_table(important_dir / "tables" / "lineage_written_summaries.tsv")
    top_prevalence = sorted(prevalence_rows, key=lambda row: (row.get("database", ""), -(_float_or_none(row.get("prevalence_percent", "")) or 0.0), row.get("feature_id", "")))[:20]
    top_geographic_burden = sorted(
        [
            row for row in geographic_burden_rows
            if row.get("geo_level") == "country" and row.get("group_name") not in {"missing", "unknown", "missing (unknown)"}
        ],
        key=lambda row: (-(_float_or_none(row.get("prevalence_percent", "")) or 0.0), row.get("database", ""), row.get("group_name", "")),
    )[:20]
    top_geographic_features = sorted(
        [
            row for row in geographic_feature_rows
            if row.get("geo_level") == "country" and row.get("group_name") not in {"missing", "unknown", "missing (unknown)"}
        ],
        key=lambda row: (-(_float_or_none(row.get("prevalence_percent", "")) or 0.0), row.get("database", ""), row.get("feature_id", ""), row.get("group_name", "")),
    )[:20]
    top_variation = sorted(variation_rows, key=lambda row: (-(_float_or_none(row.get("iqr_identity", "")) or 0.0), row.get("database", ""), row.get("feature_id", "")))[:20]
    top_temporal = sorted(temporal_rows, key=lambda row: (-abs(_float_or_none(row.get("change_percent_points", "")) or 0.0), row.get("database", ""), row.get("feature_id", "")))[:20]
    top_cooccurrence = sorted(cooccurrence_rows, key=lambda row: (-abs(_float_or_none(row.get("phi_correlation", "")) or 0.0), row.get("feature_a_database", ""), row.get("feature_a_id", "")))[:20]
    top_context = context_rows[:20]
    top_metadata_features = sorted(
        metadata_feature_rows,
        key=lambda row: (
            row.get("interpretation_label", "") not in {"strong_supported", "moderate_supported"},
            -abs(_float_or_none(row.get("prevalence_difference", "")) or 0.0),
            row.get("database", ""),
            row.get("feature_id", ""),
        ),
    )[:20]
    top_metadata_burden = sorted(
        metadata_burden_rows,
        key=lambda row: (-abs(_float_or_none(row.get("burden_difference", "")) or 0.0), row.get("database", ""), row.get("metadata_column", "")),
    )[:20]
    top_metadata_omnibus = sorted(
        metadata_burden_omnibus_rows + metadata_category_omnibus_rows,
        key=lambda row: (-abs(_float_or_none(row.get("burden_range_median", "")) or 0.0), row.get("database", ""), row.get("metadata_column", "")),
    )[:20]
    top_lineage_distribution = sorted(
        lineage_distribution_rows,
        key=lambda row: (-int(_float_or_none(row.get("total_genomes", "")) or 0), row.get("lineage_type", ""), row.get("lineage_id", "")),
    )[:20]
    top_lineage_overlap = sorted(
        lineage_overlap_rows,
        key=lambda row: (-(_float_or_none(row.get("dominant_lineage_fraction", "")) or 0.0), row.get("metadata_column", ""), row.get("metadata_group", "")),
    )[:20]
    top_lineage_burden = sorted(
        lineage_burden_rows,
        key=lambda row: (-(_float_or_none(row.get("median_features_per_genome", "")) or 0.0), row.get("database", ""), row.get("lineage_id", "")),
    )[:20]
    top_lineage_enrichment = sorted(
        lineage_enrichment_rows,
        key=lambda row: (
            row.get("interpretation_label", "") not in {"feature_lineage_specific", "feature_lineage_enriched"},
            -abs(_float_or_none(row.get("prevalence_difference", "")) or 0.0),
            row.get("database", ""),
            row.get("feature_id", ""),
        ),
    )[:20]
    top_lineage_presence = sorted(
        lineage_presence_rows,
        key=lambda row: (-(_float_or_none(row.get("prevalence_percent", "")) or 0.0), row.get("database", ""), row.get("feature_id", "")),
    )[:20]
    qc_table_html = _html_table(qc_steps, ["step_order", "qc_step", "tool", "enabled", "pass", "warning", "fail", "skipped", "status", "notes"], max_rows=20)
    prevalence_table_html = _html_table(top_prevalence, ["database", "feature_id", "feature_category", "positive_genomes", "total_genomes", "prevalence_display", "feature_rows", "mean_hits_per_positive_genome", "prevalence_label", "warning_flags"], max_rows=20)
    variation_table_html = _html_table(top_variation, ["database", "feature_id", "total_hits", "positive_genomes", "mean_hits_per_positive_genome", "median_identity", "iqr_identity", "median_coverage", "iqr_coverage", "variation_label", "warning_flags"], max_rows=20)
    temporal_table_html = _html_table(top_temporal, ["database", "feature_id", "first_year", "last_year", "first_year_prevalence_percent", "last_year_prevalence_percent", "change_percent_points", "correlation", "trend_label", "support_label", "warning_flags"], max_rows=20)
    cooccurrence_table_html = _html_table(top_cooccurrence, ["feature_a_database", "feature_a_id", "feature_b_database", "feature_b_id", "n_total", "n_both_present", "phi_correlation", "q_value", "direction", "significance_label", "warning_flags"], max_rows=20)
    context_table_html = _html_table(top_context, ["selected_database", "selected_feature", "context_database", "context_feature", "assembly_accession", "contig", "distance_bp", "evidence_level", "warning_flags"], max_rows=20)
    metadata_feature_table_html = _html_table(top_metadata_features, ["database", "feature_id", "metadata_column", "metadata_group", "group_n", "positive_in_group", "prevalence_in_group", "prevalence_difference", "odds_ratio", "q_value", "effect_size_label", "support_label", "interpretation_label", "warning_flags"], max_rows=20)
    metadata_burden_table_html = _html_table(top_metadata_burden, ["database", "metadata_column", "metadata_group", "group_n", "median_burden_group", "median_burden_outside", "burden_difference", "q_value", "support_label", "interpretation_label", "warning_flags"], max_rows=20)
    metadata_omnibus_table_html = _html_table(top_metadata_omnibus, ["database", "feature_category", "metadata_column", "groups_tested", "samples_tested", "burden_range_median", "test_name", "test_statistic", "q_value", "support_label", "interpretation_label", "warning_flags"], max_rows=20)
    lineage_distribution_table_html = _html_table(top_lineage_distribution, ["lineage_type", "lineage_id", "total_genomes", "fraction_display", "top_country", "top_source", "top_bioproject", "warning_flags"], max_rows=20)
    lineage_overlap_table_html = _html_table(top_lineage_overlap, ["lineage_type", "metadata_column", "metadata_group", "dominant_lineage", "dominant_lineage_fraction", "dominance_label", "warning_flags"], max_rows=20)
    lineage_burden_table_html = _html_table(top_lineage_burden, ["lineage_type", "lineage_id", "database", "lineage_n", "positive_genomes", "prevalence_percent", "median_features_per_genome", "support_label", "warning_flags"], max_rows=20)
    lineage_enrichment_table_html = _html_table(top_lineage_enrichment, ["lineage_type", "lineage_id", "database", "feature_id", "lineage_n", "positive_in_lineage", "prevalence_in_lineage", "prevalence_difference", "odds_ratio", "q_value", "support_label", "interpretation_label", "warning_flags"], max_rows=20)
    lineage_presence_table_html = _html_table(top_lineage_presence, ["database", "feature_id", "lineage_type", "lineage_id", "lineage_n", "positive_genomes", "prevalence_display", "top_country", "top_source", "top_bioproject", "warning_flags"], max_rows=20)
    prevalence_total_rows = sum(int(_float_or_none(row.get("feature_rows", "")) or 0) for row in prevalence_rows)
    prevalence_unique_features = len(prevalence_rows)
    prevalence_databases = len({row.get("database", "") for row in prevalence_rows if row.get("database", "")})
    prevalence_core_features = sum(int(_float_or_none(row.get("core_features", "")) or 0) for row in prevalence_core_rows)
    prevalence_common_features = sum(int(_float_or_none(row.get("common_features", "")) or 0) for row in prevalence_core_rows)
    prevalence_accessory_features = sum(int(_float_or_none(row.get("accessory_features", "")) or 0) for row in prevalence_core_rows)
    prevalence_rare_features = sum(int(_float_or_none(row.get("rare_features", "")) or 0) for row in prevalence_core_rows)
    prevalence_cards_html = (
        "<div class='cards'>"
        f"<div class='card'><span>Unique features</span><strong>{prevalence_unique_features}</strong></div>"
        f"<div class='card'><span>Feature rows</span><strong>{prevalence_total_rows}</strong></div>"
        f"<div class='card'><span>Databases</span><strong>{prevalence_databases}</strong></div>"
        f"<div class='card'><span>Core / common</span><strong>{prevalence_core_features} / {prevalence_common_features}</strong></div>"
        f"<div class='card'><span>Accessory / rare</span><strong>{prevalence_accessory_features} / {prevalence_rare_features}</strong></div>"
        "</div>"
    )
    prevalence_written_html = "".join(f"<p>{html.escape(row.get('summary', ''))}</p>" for row in prevalence_written_rows[:2]) or "<p>No prevalence written summary was generated.</p>"
    prevalence_database_table_html = _html_table(
        prevalence_database_rows,
        ["database", "total_feature_rows", "unique_features", "positive_genomes", "total_genomes", "genomes_positive_percent", "median_features_per_genome", "top_feature_id", "top_feature_prevalence_percent"],
        max_rows=20,
    )
    prevalence_figure_items = []
    for figure_name, title in [
        ("prevalence_feature_counts_by_database", "Feature counts by database"),
        ("prevalence_genomes_positive_by_database", "Genomes positive by database"),
        ("prevalence_core_accessory_rare_by_database", "Core/common/accessory/rare features"),
        ("prevalence_database_burden_by_sample", "Database burden by sample"),
    ]:
        svg_path = important_dir / "figures" / f"{figure_name}.svg"
        if not svg_path.exists():
            continue
        prevalence_figure_items.append(
            f"<div><h3>{html.escape(title)}</h3><img src='figures/{html.escape(svg_path.name)}' alt='{html.escape(title)}'>"
            f"<p><a href='figures/{html.escape(figure_name)}.png'>PNG</a> | <a href='figures/{html.escape(figure_name)}.svg'>SVG</a> | <a href='figures/{html.escape(figure_name)}.pdf'>PDF</a> | <a href='figures/{html.escape(figure_name)}.data.tsv'>Data TSV</a></p></div>"
        )
    prevalence_figures = []
    for path in sorted((important_dir / "figures").glob("prevalence_top_features_*.svg")):
        database = path.name.replace("prevalence_top_features_", "").replace(".svg", "")
        stem = path.stem
        prevalence_figures.append(
            f"<div><h3>{html.escape(database)} prevalence</h3><img src='figures/{html.escape(path.name)}' alt='{html.escape(database)} prevalence'>"
            f"<p><a href='figures/{html.escape(stem)}.png'>PNG</a> | <a href='figures/{html.escape(path.name)}'>SVG</a> | <a href='figures/{html.escape(stem)}.pdf'>PDF</a> | <a href='figures/{html.escape(stem)}.data.tsv'>Data TSV</a></p></div>"
        )
    prevalence_figures_html = "<div class='figure-row'>" + "".join((prevalence_figure_items + prevalence_figures)[:8]) + "</div>" if (prevalence_figure_items or prevalence_figures) else "<p>No prevalence figures were generated because no feature rows were available.</p>"
    geographic_databases = len({row.get("database", "") for row in geographic_burden_rows if row.get("database", "")})
    geographic_country_groups = len({row.get("group_name", "") for row in geographic_burden_rows if row.get("geo_level") == "country" and row.get("group_name") not in {"", "missing", "unknown", "missing (unknown)"}})
    geographic_missing_country = max([int(_float_or_none(row.get("missing_country_count", "")) or 0) for row in geographic_summary_rows] or [0])
    geographic_warning_count = sum(1 for row in geographic_warning_rows if row.get("warning_flags", ""))
    geographic_cards_html = (
        "<div class='cards'>"
        f"<div class='card'><span>Databases</span><strong>{geographic_databases}</strong></div>"
        f"<div class='card'><span>Country groups</span><strong>{geographic_country_groups}</strong></div>"
        f"<div class='card'><span>Missing country metadata</span><strong>{geographic_missing_country}</strong></div>"
        f"<div class='card'><span>Warning groups</span><strong>{geographic_warning_count}</strong></div>"
        "</div>"
    )
    best_geo = next((row for row in geographic_summary_rows if row.get("mode") == "database_burden" and row.get("geo_level") == "country" and row.get("top_group")), {})
    geographic_summary_html = (
        "<p>"
        f"Geographic summaries cover {geographic_databases} detected database(s) across {geographic_country_groups} country group(s). "
        f"The default country-level view highlights {html.escape(best_geo.get('database', 'detected features'))}; "
        f"the top group is {html.escape(best_geo.get('top_group', 'not available'))} "
        f"({html.escape(best_geo.get('top_group_prevalence_percent', ''))}% "
        f"{html.escape(best_geo.get('top_group_positive_genomes', ''))}/{html.escape(best_geo.get('top_group_total_genomes', ''))}). "
        f"Country metadata was missing for {geographic_missing_country} genome(s)."
        "</p>"
        if geographic_summary_rows else "<p>No geographic summaries were generated because country metadata or feature rows were unavailable.</p>"
    )
    geographic_burden_table_html = _html_table(
        top_geographic_burden,
        ["database", "geo_level", "group_name", "total_genomes", "positive_genomes", "prevalence_display", "mean_feature_burden_per_genome", "median_feature_burden_per_genome", "warning_flags"],
        max_rows=20,
    )
    geographic_feature_table_html = _html_table(
        top_geographic_features,
        ["database", "feature_id", "geo_level", "group_name", "total_genomes", "positive_genomes", "prevalence_display", "feature_rows", "mean_hits_per_positive_genome", "warning_flags"],
        max_rows=20,
    )
    geographic_figure_items = []
    for path in sorted((important_dir / "figures").glob("geographic_*_bar_*.svg"))[:6]:
        stem = path.stem
        geographic_figure_items.append(
            f"<div><h3>{html.escape(stem.replace('_', ' '))}</h3><img src='figures/{html.escape(path.name)}' alt='{html.escape(stem)}'>"
            f"<p><a href='figures/{html.escape(stem)}.png'>PNG</a> | <a href='figures/{html.escape(path.name)}'>SVG</a> | <a href='figures/{html.escape(stem)}.pdf'>PDF</a> | <a href='figures/{html.escape(stem)}.data.tsv'>Data TSV</a></p></div>"
        )
    for path in sorted((important_dir / "figures").glob("geographic_map_*.svg"))[:2]:
        stem = path.stem
        geographic_figure_items.append(
            f"<div><h3>{html.escape(stem.replace('_', ' '))}</h3><img src='figures/{html.escape(path.name)}' alt='{html.escape(stem)}'>"
            f"<p><a href='figures/{html.escape(stem)}.png'>PNG</a> | <a href='figures/{html.escape(path.name)}'>SVG</a> | <a href='figures/{html.escape(stem)}.pdf'>PDF</a> | <a href='figures/{html.escape(stem)}.data.tsv'>Data TSV</a></p></div>"
        )
    geographic_figures_html = "<div class='figure-row'>" + "".join(geographic_figure_items[:8]) + "</div>" if geographic_figure_items else "<p>No geographic figures were generated because no mappable country groups were available.</p>"
    variation_unique_features = len(variation_rows)
    variation_total_hits = sum(int(_float_or_none(row.get("total_hits", "")) or 0) for row in variation_rows)
    variation_high_features = sum(1 for row in variation_rows if row.get("variation_label") == "high_variation")
    variation_low_identity = sum(int(_float_or_none(row.get("low_identity_hits", "")) or 0) for row in variation_rows)
    variation_low_coverage = sum(int(_float_or_none(row.get("low_coverage_hits", "")) or 0) for row in variation_rows)
    variation_median_identity_values = [_float_or_none(row.get("median_identity", "")) for row in variation_rows]
    variation_median_coverage_values = [_float_or_none(row.get("median_coverage", "")) for row in variation_rows]
    variation_median_identity = _summary_stats([value for value in variation_median_identity_values if value is not None])["median"]
    variation_median_coverage = _summary_stats([value for value in variation_median_coverage_values if value is not None])["median"]
    variation_cards_html = (
        "<div class='cards'>"
        f"<div class='card'><span>Unique features analyzed</span><strong>{variation_unique_features}</strong></div>"
        f"<div class='card'><span>Total hits</span><strong>{variation_total_hits}</strong></div>"
        f"<div class='card'><span>Median identity</span><strong>{html.escape(variation_median_identity)}</strong></div>"
        f"<div class='card'><span>Median coverage</span><strong>{html.escape(variation_median_coverage)}</strong></div>"
        f"<div class='card'><span>High-variation features</span><strong>{variation_high_features}</strong></div>"
        f"<div class='card'><span>Low identity / coverage hits</span><strong>{variation_low_identity} / {variation_low_coverage}</strong></div>"
        "</div>"
    )
    variation_database_table_html = _html_table(
        variation_database_rows,
        ["database", "unique_features", "total_hits", "median_identity", "median_coverage", "high_variation_features", "low_identity_hits", "low_coverage_hits"],
        max_rows=20,
    )
    variation_figures = []
    for path in sorted((important_dir / "figures").glob("variation_*_top20.svg")):
        figure_name = path.name.replace(".svg", "")
        if figure_name.startswith("variation_identity_coverage_"):
            title = figure_name.replace("variation_identity_coverage_", "").replace("_top20", "") + " identity vs coverage"
        elif figure_name.startswith("variation_top_variable_"):
            title = figure_name.replace("variation_top_variable_", "").replace("_top20", "") + " top variable features"
        elif figure_name.startswith("variation_coverage_"):
            title = figure_name.replace("variation_coverage_", "").replace("_top20", "") + " coverage variation"
        elif figure_name.startswith("variation_identity_"):
            title = figure_name.replace("variation_identity_", "").replace("_top20", "") + " identity variation"
        else:
            title = figure_name
        stem = path.stem
        variation_figures.append(
            f"<div><h3>{html.escape(title)}</h3><img src='figures/{html.escape(path.name)}' alt='{html.escape(title)}'>"
            f"<p><a href='figures/{html.escape(stem)}.png'>PNG</a> | <a href='figures/{html.escape(path.name)}'>SVG</a> | <a href='figures/{html.escape(stem)}.pdf'>PDF</a> | <a href='figures/{html.escape(stem)}.data.tsv'>Data TSV</a></p></div>"
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
    metadata_summary = {row.get("metric", ""): row.get("value", "") for row in metadata_summary_rows}
    metadata_usable = metadata_summary.get("metadata_columns_usable", "0")
    metadata_sparse = metadata_summary.get("metadata_columns_sparse_or_biased", "0")
    metadata_excluded = metadata_summary.get("metadata_columns_excluded", "0")
    metadata_summary_html = (
        "<p>"
        f"Metadata association screening evaluated {html.escape(metadata_summary.get('metadata_columns_screened', '0'))} metadata columns, "
        f"{html.escape(metadata_summary.get('feature_enrichment_rows', '0'))} feature-by-group comparisons, and "
        f"{html.escape(metadata_summary.get('database_burden_associations', '0'))} database-burden comparisons. "
        f"Strong feature associations: {html.escape(metadata_summary.get('strong_feature_associations', '0'))}; "
        f"moderate feature associations: {html.escape(metadata_summary.get('moderate_feature_associations', '0'))}; "
        f"warning rows: {html.escape(metadata_summary.get('warning_rows', '0'))}."
        "</p>"
    )
    metadata_usability_cards = (
        "<div class='cards'>"
        f"<div class='card'><span>Usable metadata columns</span><strong>{html.escape(metadata_usable)}</strong></div>"
        f"<div class='card'><span>Sparse or biased columns</span><strong>{html.escape(metadata_sparse)}</strong></div>"
        f"<div class='card'><span>Columns excluded</span><strong>{html.escape(metadata_excluded)}</strong></div>"
        "</div>"
    )
    metadata_usability_table_html = _html_table(
        metadata_usability_rows,
        ["metadata_column", "non_missing_count", "missing_fraction", "unique_values", "largest_group", "largest_group_fraction", "eligible_for_testing", "recommended_use", "reason"],
        max_rows=20,
    )
    metadata_figure_items = []
    for key, title in [
        ("important_metadata_volcano_svg", "Metadata enrichment volcano"),
        ("important_metadata_enrichment_heatmap_svg", "Metadata enrichment heatmap"),
        ("important_metadata_burden_boxplot_svg", "Database burden by group"),
        ("important_metadata_category_enrichment_svg", "Category enrichment"),
    ]:
        figure_path = Path(metadata_association_outputs.get(key, ""))
        if not figure_path.exists():
            continue
        stem = figure_path.stem
        metadata_figure_items.append(
            f"<div><h3>{html.escape(title)}</h3><img src='figures/{html.escape(figure_path.name)}' alt='{html.escape(title)}'>"
            f"<p><a href='figures/{html.escape(stem)}.png'>PNG</a> | <a href='figures/{html.escape(figure_path.name)}'>SVG</a> | <a href='figures/{html.escape(stem)}.pdf'>PDF</a> | <a href='figures/{html.escape(stem)}.data.tsv'>Data TSV</a></p></div>"
        )
    metadata_figures_html = "<div class='figure-row'>" + "".join(metadata_figure_items) + "</div>" if metadata_figure_items else "<p>No metadata association figures were generated because metadata groups or feature rows were unavailable.</p>"
    lineage_summary = {row.get("metric", ""): row.get("value", "") for row in lineage_summary_rows}
    lineage_cards_html = (
        "<div class='cards'>"
        f"<div class='card'><span>Genomes with MLST ST</span><strong>{html.escape(lineage_summary.get('genomes_with_mlst_ST', '0'))}</strong></div>"
        f"<div class='card'><span>Unique STs</span><strong>{html.escape(lineage_summary.get('unique_STs', '0'))}</strong></div>"
        f"<div class='card'><span>Dominant ST</span><strong>{html.escape(lineage_summary.get('dominant_ST', '')) or 'NA'}</strong></div>"
        f"<div class='card'><span>ANI clusters</span><strong>{html.escape(lineage_summary.get('unique_ANI_clusters', '0'))}</strong></div>"
        f"<div class='card'><span>BioProjects</span><strong>{html.escape(lineage_summary.get('bioprojects_detected', '0'))}</strong></div>"
        f"<div class='card'><span>Lineage warnings</span><strong>{html.escape(lineage_summary.get('warning_rows', '0'))}</strong></div>"
        "</div>"
    )
    lineage_summary_html = (
        "<p>"
        f"Lineage reporting found {html.escape(lineage_summary.get('lineage_types_available', 'none'))} context, "
        f"{html.escape(lineage_summary.get('genomes_with_mlst_ST', '0'))} genome(s) with MLST ST calls, and "
        f"{html.escape(lineage_summary.get('genomes_with_ANI_cluster', '0'))} genome(s) with ANI cluster context. "
        f"The default view is {html.escape(lineage_summary.get('default_view', ''))}. "
        "Use this section to check whether feature, geography, temporal, or metadata patterns are concentrated in one ST, ANI cluster, or BioProject."
        "</p>"
    )
    lineage_written_html = "".join(
        f"<p><strong>{html.escape(row.get('section', '').replace('_', ' ').title())}:</strong> {html.escape(row.get('summary', ''))}</p>"
        for row in lineage_written_rows
    ) or lineage_summary_html
    lineage_figure_items = []
    for key, title in [
        ("important_lineage_distribution_svg", "Lineage distribution"),
        ("important_lineage_metadata_overlap_svg", "Metadata-lineage overlap"),
        ("important_lineage_database_burden_svg", "Database burden by lineage"),
        ("important_lineage_feature_heatmap_svg", "Feature prevalence by lineage"),
        ("important_lineage_feature_enrichment_svg", "Lineage feature enrichment"),
        ("important_lineage_feature_presence_svg", "Selected feature by lineage"),
        ("important_lineage_confounding_top_findings_svg", "Lineage-adjusted top findings"),
    ]:
        figure_path = Path(lineage_outputs.get(key, ""))
        if not figure_path.exists():
            continue
        stem = figure_path.stem
        lineage_figure_items.append(
            f"<div><h3>{html.escape(title)}</h3><img src='figures/{html.escape(figure_path.name)}' alt='{html.escape(title)}'>"
            f"<p><a href='figures/{html.escape(stem)}.png'>PNG</a> | <a href='figures/{html.escape(figure_path.name)}'>SVG</a> | <a href='figures/{html.escape(stem)}.pdf'>PDF</a> | <a href='figures/{html.escape(stem)}.data.tsv'>Data TSV</a></p></div>"
        )
    lineage_figures_html = "<div class='figure-row'>" + "".join(lineage_figure_items) + "</div>" if lineage_figure_items else "<p>No lineage figures were generated because lineage context was unavailable.</p>"
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
<a href="#metadata-associations">Metadata Associations</a>
<a href="#lineage">Lineage / Clonal Structure</a>
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
<section id="prevalence"><h2>Prevalence</h2><p>This section summarizes how frequently each detected feature appears across analyzed genomes. Prevalence uses positive genome count; feature rows can be higher when a feature appears multiple times in one genome.</p>
<div class="warning">Prevalence reflects the analyzed dataset, not global prevalence. This section is descriptive and does not test association or causality.</div>
{prevalence_cards_html}
{prevalence_written_html}
<iframe src="figures/prevalence_analysis.html" title="Feature prevalence interactive report"></iframe>
{prevalence_figures_html}
<h3>Database summary</h3>
{prevalence_database_table_html}
<h3>Top feature prevalence</h3>
{prevalence_table_html}
<div class="downloads"><a href="figures/prevalence_analysis.html">Open interactive prevalence report</a><a href="prevalence_tables.zip">Download prevalence tables ZIP</a><a href="prevalence_figures.zip">Download prevalence figures ZIP</a><a href="tables/feature_prevalence.tsv">Download full feature prevalence</a><a href="tables/feature_prevalence_top.tsv">Download top feature prevalence</a><a href="tables/prevalence_summary_by_database.tsv">Download database summary</a><a href="tables/prevalence_core_accessory_rare_summary.tsv">Download core/common/accessory/rare summary</a><a href="tables/prevalence_database_burden_by_sample.tsv">Download database burden by sample</a></div></section>
<section id="geography"><h2>Geographic Distribution</h2><div class="warning">Geographic patterns reflect the analyzed dataset only. They are not global prevalence estimates and can be affected by BioProject, lineage, country, and year sampling bias.</div>
{geographic_cards_html}
{geographic_summary_html}
<iframe src="figures/geographic_distribution.html" title="Geographic distribution interactive report"></iframe>
{geographic_figures_html}
<h3>Top database burden by country</h3>
{geographic_burden_table_html}
<h3>Top feature distributions by country</h3>
{geographic_feature_table_html}
<div class="downloads"><a href="figures/geographic_distribution.html">Open interactive geographic report</a><a href="figures/geographic_distribution_map.html">Open compatibility map</a><a href="geographic_tables.zip">Download geographic tables ZIP</a><a href="geographic_figures.zip">Download geographic figures ZIP</a><a href="tables/geographic_distribution_summary.tsv">Download geographic summary</a><a href="tables/geographic_database_burden.tsv">Download database burden table</a><a href="tables/geographic_feature_distribution.tsv">Download feature distribution table</a><a href="tables/geographic_warning_summary.tsv">Download warning summary</a></div></section>
<section id="variations"><h2>Variations</h2><p>Variation summaries use identity, coverage, alignment length, and hit-count values when available. Low identity, low coverage, high variation, and few-hit flags are review cues, not automatic failures.</p>
<div class="warning">A feature can be common but conserved, common and variable, or rare with unstable estimates. Use the complete hit-level table when reviewing low-confidence or partial hits.</div>
{variation_cards_html}
<iframe src="figures/variation_analysis.html" title="Feature variation interactive report"></iframe>
{variation_figures_html}
<h3>Variation by database</h3>
{variation_database_table_html}
<h3>Most variable features</h3>
{variation_table_html}
<div class="downloads"><a href="figures/variation_analysis.html">Open interactive variation report</a><a href="variation_figures.zip">Download variation figures ZIP</a><a href="key_tables/feature_variation_database_summary.tsv">Download variation database summary</a><a href="key_tables/feature_variation_summary.tsv">Download variation summary</a><a href="key_tables/feature_variation_hits.tsv">Download hit-level variation table</a></div></section>
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
<section id="metadata-associations"><h2>Metadata Associations</h2><div class="warning">Metadata associations are exploratory enrichment-style screens. They may reflect sampling, BioProject structure, lineage composition, geography, collection year, or missing metadata and should not be interpreted as causal.</div>
{metadata_summary_html}
<h3>Metadata usability</h3>
{metadata_usability_cards}
{metadata_usability_table_html}
<iframe src="figures/metadata_associations.html" title="Metadata associations interactive report"></iframe>
{metadata_figures_html}
<h3>Top feature associations</h3>
{metadata_feature_table_html}
<h3>Top database-burden associations</h3>
{metadata_burden_table_html}
<h3>Top multi-group burden tests</h3>
{metadata_omnibus_table_html}
<div class="downloads"><a href="figures/metadata_associations.html">Open interactive metadata association report</a><a href="metadata_association_tables.zip">Download all metadata association tables ZIP</a><a href="metadata_association_figures.zip">Download all metadata association figures ZIP</a><a href="tables/metadata_usability_summary.tsv">Download metadata usability summary</a><a href="tables/metadata_feature_enrichment.tsv">Download feature enrichment</a><a href="tables/metadata_burden_associations.tsv">Download burden associations</a><a href="tables/metadata_category_enrichment.tsv">Download category enrichment</a><a href="tables/metadata_burden_omnibus.tsv">Download burden omnibus tests</a><a href="tables/metadata_category_omnibus.tsv">Download category omnibus tests</a><a href="tables/metadata_association_summary.tsv">Download metadata association summary</a></div></section>
<section id="lineage"><h2>Lineage / Clonal Structure</h2><div class="warning">Lineage summaries are exploratory and do not replace phylogenetic analysis. Apparent metadata associations may reflect clonal structure, BioProject sampling, geography, or temporal sampling.</div>
{lineage_cards_html}
{lineage_summary_html}
{lineage_written_html}
<iframe src="figures/lineage_clonal_structure.html" title="Lineage and clonal structure interactive report"></iframe>
{lineage_figures_html}
<h3>Lineage distribution</h3>
{lineage_distribution_table_html}
<h3>Metadata-lineage overlap</h3>
{lineage_overlap_table_html}
<h3>Feature burden by lineage</h3>
{lineage_burden_table_html}
<h3>Feature enrichment by lineage</h3>
{lineage_enrichment_table_html}
<h3>Selected feature lineage report</h3>
{lineage_presence_table_html}
<div class="downloads"><a href="figures/lineage_clonal_structure.html">Open interactive lineage report</a><a href="lineage_tables.zip">Download lineage tables ZIP</a><a href="lineage_figures.zip">Download lineage figures ZIP</a><a href="tables/lineage_summary.tsv">Download sample lineage summary</a><a href="tables/lineage_distribution.tsv">Download lineage distribution</a><a href="tables/lineage_metadata_overlap.tsv">Download metadata-lineage overlap</a><a href="tables/lineage_feature_burden.tsv">Download lineage feature burden</a><a href="tables/lineage_feature_enrichment.tsv">Download lineage feature enrichment</a><a href="tables/lineage_adjusted_top_findings.tsv">Download lineage-adjusted top findings</a><a href="tables/lineage_feature_presence.tsv">Download selected feature lineage table</a><a href="tables/lineage_written_summaries.tsv">Download written summaries</a></div></section>
<section id="files"><h2>Important Files</h2><ul>
<li><a href="../basic/enriched_genome_dataset.csv">Enriched genome dataset CSV</a></li>
<li><a href="key_tables/qc_step_summary.tsv">QC step summary</a></li>
<li><a href="tables/feature_prevalence.tsv">Feature prevalence table</a></li>
<li><a href="tables/prevalence_summary_by_database.tsv">Prevalence database summary</a></li>
<li><a href="figures/prevalence_analysis.html">Interactive prevalence report</a></li>
<li><a href="figures/geographic_distribution.html">Interactive geographic report</a></li>
<li><a href="tables/geographic_database_burden.tsv">Geographic database burden table</a></li>
<li><a href="tables/geographic_feature_distribution.tsv">Geographic feature distribution table</a></li>
<li><a href="key_tables/feature_variation_summary.tsv">Feature variation summary</a></li>
<li><a href="key_tables/feature_variation_database_summary.tsv">Feature variation database summary</a></li>
<li><a href="figures/variation_analysis.html">Interactive variation report</a></li>
<li><a href="key_tables/temporal_trend_summary.tsv">Temporal trend summary</a></li>
<li><a href="tables/cooccurrence_pair_summary.tsv">Co-occurrence pair summary</a></li>
<li><a href="tables/genomic_context_evidence.tsv">Genomic context evidence</a></li>
<li><a href="tables/metadata_feature_enrichment.tsv">Metadata feature enrichment</a></li>
<li><a href="tables/metadata_burden_associations.tsv">Metadata burden associations</a></li>
<li><a href="tables/metadata_usability_summary.tsv">Metadata usability summary</a></li>
<li><a href="tables/metadata_burden_omnibus.tsv">Metadata burden omnibus tests</a></li>
<li><a href="tables/metadata_category_omnibus.tsv">Metadata category omnibus tests</a></li>
<li><a href="figures/lineage_clonal_structure.html">Interactive lineage / clonal structure report</a></li>
<li><a href="tables/lineage_distribution.tsv">Lineage distribution</a></li>
<li><a href="tables/lineage_metadata_overlap.tsv">Metadata-lineage overlap</a></li>
<li><a href="tables/lineage_feature_enrichment.tsv">Lineage feature enrichment</a></li>
<li><a href="tables/lineage_adjusted_top_findings.tsv">Lineage-adjusted top findings</a></li>
<li><a href="tables/lineage_written_summaries.tsv">Lineage written summaries</a></li>
<li><a href="../panr2_inputs/features/all_features.tsv">Complete standardized feature table</a></li>
<li><a href="../panr2_inputs/manifest/schema_validation_summary.txt">Feature-contract validation summary</a></li>
</ul></section>
<section id="warnings"><h2>Warnings And Limitations</h2><p>Association, geography, and co-occurrence summaries are exploratory. Confirm important findings with denominator checks, lineage context, BioProject balance, and independent datasets.</p></section>
</main></body></html>
""",
        encoding="utf-8",
    )
    return {
        "important_results_html": str(report_path),
        **geographic_outputs,
        **qc_outputs,
        **prevalence_outputs,
        **variation_outputs,
        **temporal_outputs,
        **cooccurrence_outputs,
        **metadata_association_outputs,
        **lineage_outputs,
    }


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
        metadata_association_outputs = write_important_metadata_association_outputs(sample_dir, out_dir, important_dir)
        lineage_outputs = write_important_lineage_outputs(sample_dir, out_dir, important_dir)
        outputs.update(write_important_results_report(sample_dir, out_dir, important_dir, geographic_outputs, qc_outputs, prevalence_outputs, variation_outputs, temporal_outputs, cooccurrence_outputs, metadata_association_outputs, lineage_outputs))
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
