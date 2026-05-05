#!/usr/bin/env python3
"""Download NCBI assembly FASTA files from FetchM2 clean metadata.

This avoids FetchM2 0.1.3's threaded SQLite cache issue while preserving the
same sequence output files expected by PanResistome.
"""

from __future__ import annotations

import argparse
import gzip
import re
import shutil
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import requests


BASE_URL = "https://ftp.ncbi.nlm.nih.gov/genomes/all"


def normalize_text(value) -> str:
    return str(value or "").strip().lower()


def normalize_assembly_name(name: str) -> str:
    cleaned = str(name or "").strip()
    return cleaned.replace(" ", "_") if cleaned else "NA"


def build_parent_url(accession: str) -> str:
    prefix, digits = accession.split("_", 1)
    core = digits.split(".", 1)[0]
    return f"{BASE_URL}/{prefix}/{core[:3]}/{core[3:6]}/{core[6:9]}/{core[9:]}"


def resolve_assembly_directory(accession: str, name: str, session: requests.Session) -> str:
    parent_url = build_parent_url(accession)
    normalized_name = normalize_assembly_name(name)
    for candidate in [f"{accession}_{normalized_name}", f"{accession}_NA"]:
        response = session.get(f"{parent_url}/{candidate}/", timeout=30)
        if response.ok:
            return candidate
    response = session.get(parent_url, timeout=60)
    response.raise_for_status()
    matches = [item.rstrip("/") for item in re.findall(r'href="([^"]+/)"', response.text) if item.startswith(f"{accession}_")]
    if not matches:
        raise FileNotFoundError(f"No remote assembly directory found for {accession}")
    return matches[0]


def row_matches_filters(row: dict, filters: dict) -> bool:
    for field, values in {
        "Country": filters.get("country"),
        "Continent": filters.get("continent"),
        "Subcontinent": filters.get("subcontinent"),
        "Host_SD": filters.get("host"),
        "Host_Rank": filters.get("host_rank"),
        "Sample_Type_SD": filters.get("sample_type"),
        "Isolation_Source_SD": filters.get("isolation_source"),
        "Environment_Medium_SD": filters.get("environment_medium"),
    }.items():
        if values and normalize_text(row.get(field)) not in {normalize_text(value) for value in values}:
            return False
    year_from = filters.get("year_from")
    year_to = filters.get("year_to")
    if year_from is not None or year_to is not None:
        try:
            year = int(str(row.get("Collection_Year") or row.get("Collection Date") or "")[:4])
        except ValueError:
            return False
        if year_from is not None and year < int(year_from):
            return False
        if year_to is not None and year > int(year_to):
            return False
    return True


def select_rows(input_path: Path, filters: dict, max_genomes: int | None) -> list[dict]:
    df = pd.read_csv(input_path)
    rows = [row for row in df.fillna("").to_dict(orient="records") if row_matches_filters(row, filters)]
    if max_genomes is not None:
        rows = rows[:max_genomes]
    return rows


def download_one(row: dict, outdir: Path, retries: int, retry_delay: float, keep_gz: bool) -> tuple[str, str]:
    accession = str(row.get("Assembly Accession") or "").strip()
    name = str(row.get("Assembly Name") or "").strip()
    if not accession:
        return "", "missing accession"
    session = requests.Session()
    for attempt in range(1, retries + 1):
        try:
            directory = resolve_assembly_directory(accession, name, session)
            gz_name = f"{directory}_genomic.fna.gz"
            fna_name = f"{directory}_genomic.fna"
            gz_path = outdir / gz_name
            fna_path = outdir / fna_name
            if fna_path.exists() or gz_path.exists():
                return accession, "exists"
            url = f"{build_parent_url(accession)}/{directory}/{gz_name}"
            with session.get(url, stream=True, timeout=300) as response:
                response.raise_for_status()
                with gz_path.open("wb") as handle:
                    for chunk in response.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            handle.write(chunk)
            if not keep_gz:
                with gzip.open(gz_path, "rb") as source, fna_path.open("wb") as target:
                    shutil.copyfileobj(source, target)
                gz_path.unlink()
            return accession, "downloaded"
        except Exception as exc:
            if attempt >= retries:
                return accession, f"failed: {exc}"
            time.sleep(retry_delay * attempt)
    return accession, "failed"


def main() -> None:
    parser = argparse.ArgumentParser(description="Download FASTA assemblies from FetchM2 clean metadata.")
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--retries", type=int, default=3)
    parser.add_argument("--retry-delay", type=float, default=5.0)
    parser.add_argument("--max-genomes", type=int)
    parser.add_argument("--keep-gz", action="store_true")
    parser.add_argument("--host", nargs="+")
    parser.add_argument("--host-rank", nargs="+")
    parser.add_argument("--country", nargs="+")
    parser.add_argument("--continent", nargs="+")
    parser.add_argument("--subcontinent", nargs="+")
    parser.add_argument("--sample-type", nargs="+")
    parser.add_argument("--isolation-source", nargs="+")
    parser.add_argument("--environment-medium", nargs="+")
    parser.add_argument("--year-from", type=int)
    parser.add_argument("--year-to", type=int)
    args = parser.parse_args()

    args.outdir.mkdir(parents=True, exist_ok=True)
    filters = {
        "host": args.host,
        "host_rank": args.host_rank,
        "country": args.country,
        "continent": args.continent,
        "subcontinent": args.subcontinent,
        "sample_type": args.sample_type,
        "isolation_source": args.isolation_source,
        "environment_medium": args.environment_medium,
        "year_from": args.year_from,
        "year_to": args.year_to,
    }
    rows = select_rows(args.input, filters, args.max_genomes)
    results: list[tuple[str, str]] = []
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as executor:
        futures = [executor.submit(download_one, row, args.outdir, args.retries, args.retry_delay, args.keep_gz) for row in rows]
        for future in as_completed(futures):
            results.append(future.result())

    failed = [accession for accession, status in results if status.startswith("failed") or status == "missing accession"]
    (args.outdir / "failed_accessions.txt").write_text("\n".join(failed) + ("\n" if failed else ""), encoding="utf-8")
    pd.DataFrame(
        [{"assembly_accession": accession, "status": status} for accession, status in results],
        columns=["assembly_accession", "status"],
    ).to_csv(args.outdir / "sequence_download_summary.csv", index=False)
    print({
        "selected": len(rows),
        "downloaded": sum(1 for _, status in results if status == "downloaded"),
        "existing": sum(1 for _, status in results if status == "exists"),
        "failed": len(failed),
    })


if __name__ == "__main__":
    main()
