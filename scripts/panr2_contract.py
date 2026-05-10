#!/usr/bin/env python3
"""Create and validate PanR2 contract feature tables from PanResistome outputs."""

from __future__ import annotations

import csv
import html
import itertools
import math
import re
from collections import Counter, defaultdict
from pathlib import Path


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


def discover_raw_feature_databases(sample_dir: Path) -> set[str]:
    raw_databases: set[str] = set()
    mlst_dirs = [
        sample_dir / "mlst",
        sample_dir / "tool_results" / "mlst",
    ]
    if any(path.exists() and any(child.is_file() for child in path.rglob("*")) for path in mlst_dirs):
        raw_databases.add("mlst")
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
                if any(path.name.endswith(suffix) for suffix in skip_suffixes):
                    continue
                if path.name.endswith(".features.tsv"):
                    continue
                paths.append(path)
    return paths


def discover_feature_rows(sample_dir: Path) -> list[dict[str, str]]:
    sample_map = load_sample_map(sample_dir)
    rows: list[dict[str, str]] = []
    seen_paths: set[Path] = set()

    for path in sorted(sample_dir.rglob("*_results.*")):
        if path.suffix.lower() not in {".tab", ".tsv", ".csv"}:
            continue
        if "panr2_inputs" in path.parts and "features" in path.parts:
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
) -> list[dict[str, str]]:
    annotated = []
    for finding in top_findings:
        support_samples = _supporting_samples_for_finding(finding, normalized_metadata_rows, presence)
        project_summary = _top_project_summary(support_samples, project_by_sample)
        flags = []
        if finding.get("warning"):
            flags.extend(flag.strip() for flag in finding["warning"].split(";") if flag.strip())
        if not finding.get("p_value") and not finding.get("q_value"):
            flags.append("screening_no_p_value")
        if 0 < len(support_samples) < 5:
            flags.append("low_sample_warning")
        if project_summary["warning"]:
            flags.append(project_summary["warning"])
        warning_flags = _flag_string(flags)
        out = dict(finding)
        out.update({
            "supporting_samples": str(len(support_samples)),
            "largest_bioproject": project_summary["largest_bioproject"],
            "largest_bioproject_fraction": project_summary["largest_bioproject_fraction"],
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
    top_findings = _annotate_top_findings(top_findings, normalized, presence, _bioproject_by_sample(metadata_rows))
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
    audit_rows = read_table(manifest_dir / "feature_completeness_audit.tsv")
    setup_body = "<h2>Database And Tool Setup</h2>"
    setup_body += _html_table(
        setup_rows,
        ["database_or_tool", "required_for_profile", "checked", "status", "setup_action", "version_or_path", "message"],
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
    cross_outputs = write_cross_database_outputs(
        all_features,
        metadata_rows,
        out_dir,
        max_features=max_features_network,
    )
    audit_outputs = write_feature_completeness_audit(sample_dir, out_dir)
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
    return {
        **written,
        **validation,
        **metadata_outputs,
        **feature_outputs,
        **matrix_outputs,
        **cross_outputs,
        **audit_outputs,
        **control_outputs,
        **report_outputs,
    }
