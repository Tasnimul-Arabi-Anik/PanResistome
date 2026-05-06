#!/usr/bin/env python3
"""Generate a FetchM2/PanResistome input TSV from NCBI Assembly E-utilities."""

from __future__ import annotations

import argparse
import csv
import json
import re
import time
import urllib.parse
import urllib.request
from datetime import date
from pathlib import Path


NCBI_COLUMNS = [
    "Assembly Accession",
    "Assembly Name",
    "Organism Name",
    "Organism Taxonomic ID",
    "ANI Check status",
    "Organism Infraspecific Names Breed",
    "Organism Infraspecific Names Strain",
    "Organism Infraspecific Names Cultivar",
    "Organism Infraspecific Names Ecotype",
    "Organism Infraspecific Names Isolate",
    "Organism Infraspecific Names Sex",
    "Annotation Name",
    "Assembly Stats Total Sequence Length",
    "Assembly Stats Total Number of Chromosomes",
    "Assembly Level",
    "Assembly Release Date",
    "WGS project accession",
    "Assembly Stats Contig N50",
    "Assembly Stats Scaffold N50",
    "Assembly Stats Number of Scaffolds",
    "Annotation BUSCO Complete",
    "Annotation BUSCO Single Copy",
    "Annotation BUSCO Duplicated",
    "Annotation BUSCO Fragmented",
    "Annotation BUSCO Missing",
    "Annotation BUSCO Lineage",
    "Assembly Submitter",
    "Assembly BioProject Accession",
    "Assembly BioSample Accession",
    "Annotation Count Gene Total",
    "Annotation Count Gene Protein-coding",
    "Annotation Count Gene Pseudogene",
    "Type Material Display Text",
    "CheckM marker set",
    "CheckM completeness",
    "CheckM contamination",
    "RefSeq Category",
    "Genome Representation",
    "FTP Path RefSeq",
    "FTP Path GenBank",
]


def fetch_json(url: str) -> dict:
    with urllib.request.urlopen(url, timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def ncbi_url(endpoint: str, params: dict[str, str]) -> str:
    return f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/{endpoint}?{urllib.parse.urlencode(params)}"


def stat_from_meta(meta: str, category: str) -> str:
    match = re.search(rf'<Stat category="{re.escape(category)}" sequence_tag="all">([^<]+)</Stat>', meta or "")
    return match.group(1) if match else ""


def first_bioproject(record: dict, key: str) -> str:
    projects = record.get(key) or []
    if projects:
        return str(projects[0].get("bioprojectaccn", "") or "")
    return ""


def record_to_row(record: dict) -> dict[str, str]:
    biosource = record.get("biosource") or {}
    infraspecies = {
        str(item.get("sub_type", "")).lower(): str(item.get("sub_value", "") or "")
        for item in biosource.get("infraspecieslist", []) or []
        if isinstance(item, dict)
    }
    busco = record.get("busco") or {}
    meta = record.get("meta", "")
    accession = str(record.get("assemblyaccession", "") or "")
    return {
        "Assembly Accession": accession,
        "Assembly Name": str(record.get("assemblyname", "") or ""),
        "Organism Name": str(record.get("speciesname", "") or record.get("organism", "") or ""),
        "Organism Taxonomic ID": str(record.get("taxid", "") or ""),
        "ANI Check status": "all",
        "Organism Infraspecific Names Breed": infraspecies.get("breed", ""),
        "Organism Infraspecific Names Strain": infraspecies.get("strain", ""),
        "Organism Infraspecific Names Cultivar": infraspecies.get("cultivar", ""),
        "Organism Infraspecific Names Ecotype": infraspecies.get("ecotype", ""),
        "Organism Infraspecific Names Isolate": infraspecies.get("isolate", str(biosource.get("isolate", "") or "")),
        "Organism Infraspecific Names Sex": str(biosource.get("sex", "") or ""),
        "Annotation Name": "",
        "Assembly Stats Total Sequence Length": stat_from_meta(meta, "total_length"),
        "Assembly Stats Total Number of Chromosomes": stat_from_meta(meta, "chromosome_count"),
        "Assembly Level": str(record.get("assemblystatus", "") or ""),
        "Assembly Release Date": str(record.get("asmreleasedate_refseq", "") or record.get("asmreleasedate_genbank", "") or "").split()[0].replace("/", "-"),
        "WGS project accession": str(record.get("wgs", "") or ""),
        "Assembly Stats Contig N50": str(record.get("contign50", "") or stat_from_meta(meta, "contig_n50")),
        "Assembly Stats Scaffold N50": str(record.get("scaffoldn50", "") or stat_from_meta(meta, "scaffold_n50")),
        "Assembly Stats Number of Scaffolds": stat_from_meta(meta, "scaffold_count"),
        "Annotation BUSCO Complete": str(busco.get("complete", "") or ""),
        "Annotation BUSCO Single Copy": str(busco.get("singlecopy", "") or ""),
        "Annotation BUSCO Duplicated": str(busco.get("duplicated", "") or ""),
        "Annotation BUSCO Fragmented": str(busco.get("fragmented", "") or ""),
        "Annotation BUSCO Missing": str(busco.get("missing", "") or ""),
        "Annotation BUSCO Lineage": str(busco.get("buscolineage", "") or ""),
        "Assembly Submitter": str(record.get("submitterorganization", "") or ""),
        "Assembly BioProject Accession": first_bioproject(record, "rs_bioprojects") or first_bioproject(record, "gb_bioprojects"),
        "Assembly BioSample Accession": str(record.get("biosampleaccn", "") or ""),
        "Annotation Count Gene Total": "",
        "Annotation Count Gene Protein-coding": "",
        "Annotation Count Gene Pseudogene": "",
        "Type Material Display Text": str(record.get("fromtype", "") or ""),
        "CheckM marker set": "",
        "CheckM completeness": "",
        "CheckM contamination": "",
        "RefSeq Category": str(record.get("refseq_category", "") or ""),
        "Genome Representation": "full" if "full-genome-representation" in (record.get("propertylist") or []) else "",
        "FTP Path RefSeq": str(record.get("ftppath_refseq", "") or ""),
        "FTP Path GenBank": str(record.get("ftppath_genbank", "") or ""),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Create an NCBI Assembly TSV input for PanResistome/FetchM2.")
    parser.add_argument("--organism", default="Delftia tsuruhatensis", help="NCBI organism query.")
    parser.add_argument("--outdir", required=True, type=Path, help="Validation input directory.")
    parser.add_argument("--max-records", type=int, default=0, help="Maximum records to write. 0 means all records returned by NCBI.")
    parser.add_argument("--sleep", type=float, default=0.34, help="Delay between NCBI requests.")
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    query = f'"{args.organism}"[Organism]'
    search = fetch_json(
        ncbi_url(
            "esearch.fcgi",
            {"db": "assembly", "term": query, "retmax": str(args.max_records or 10000), "retmode": "json"},
        )
    )
    ids = search.get("esearchresult", {}).get("idlist", [])
    if args.max_records:
        ids = ids[: args.max_records]
    if not ids:
        raise SystemExit(f"No NCBI Assembly records found for {args.organism!r}")

    records = []
    raw_summaries = {"organism": args.organism, "query": query, "generated_on": date.today().isoformat(), "records": []}
    for offset in range(0, len(ids), 100):
        batch = ids[offset : offset + 100]
        time.sleep(args.sleep)
        summary = fetch_json(
            ncbi_url(
                "esummary.fcgi",
                {"db": "assembly", "id": ",".join(batch), "retmode": "json"},
            )
        )
        result = summary.get("result", {})
        for uid in result.get("uids", []):
            record = result.get(uid, {})
            if record:
                raw_summaries["records"].append(record)
                records.append(record_to_row(record))

    records = sorted(records, key=lambda row: (not row["Assembly Accession"].startswith("GCF_"), row["Assembly Accession"]))
    tsv_path = args.outdir / "ncbi_dataset.tsv"
    with tsv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=NCBI_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)

    raw_path = args.outdir / "ncbi_assembly_esummary.json"
    raw_path.write_text(json.dumps(raw_summaries, indent=2, sort_keys=True), encoding="utf-8")

    readme_path = args.outdir / "README.md"
    readme_path.write_text(
        "\n".join(
            [
                "# Delftia tsuruhatensis Validation Input",
                "",
                f"Generated on: {date.today().isoformat()}",
                "",
                f"NCBI E-utilities query: `{query}`",
                f"Assembly records written: `{len(records)}`",
                "",
                "Files:",
                "- `ncbi_dataset.tsv`: FetchM2/PanResistome-compatible Assembly TSV.",
                "- `ncbi_assembly_esummary.json`: raw NCBI Assembly esummary JSON used to create the TSV.",
                "",
                "Recommended validation command:",
                "",
                "```bash",
                "nextflow run main.nf \\",
                f"  --input {tsv_path} \\",
                "  --outdir validation_runs/delftia_current \\",
                "  -profile conda,mamba \\",
                "  --analysis_profile comprehensive \\",
                "  --qc_filter true \\",
                "  --run_gtdbtk false \\",
                "  --run_quast true \\",
                "  --run_ani true \\",
                "  --run_mash true \\",
                "  --run_amrfinderplus true \\",
                "  --threads 8",
                "```",
                "",
                "The primary combined HTML output from a successful comprehensive run is:",
                "",
                "```text",
                "validation_runs/delftia_current/<organism>/report/index.html",
                "```",
                "",
            ]
        ),
        encoding="utf-8",
    )
    print(f"Wrote {len(records)} records to {tsv_path}")


if __name__ == "__main__":
    main()
