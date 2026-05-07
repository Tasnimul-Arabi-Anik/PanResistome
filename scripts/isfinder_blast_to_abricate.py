#!/usr/bin/env python3
"""Convert BLAST hits against an authorized ISfinder FASTA into PanR2 tables."""

import argparse
import csv
from pathlib import Path


ABRICATE_FIELDS = [
    "#FILE",
    "SEQUENCE",
    "START",
    "END",
    "GENE",
    "COVERAGE",
    "COVERAGE_MAP",
    "GAPS",
    "%COVERAGE",
    "%IDENTITY",
    "DATABASE",
    "ACCESSION",
    "PRODUCT",
    "RESISTANCE",
]


BLAST_FIELDS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "qlen",
    "slen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
]


def clean_feature_id(subject: str) -> str:
    subject = str(subject or "").strip()
    if not subject:
        return "unknown_is"
    token = subject.split()[0]
    for sep in ("|", ";", ","):
        if sep in token:
            token = token.split(sep)[0]
    return token or subject


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blast", required=True, help="BLAST outfmt 6 table.")
    parser.add_argument("--sample-id", required=True, help="FASTA/sample identifier.")
    parser.add_argument("--out-results", required=True, help="ABRicate-style results table.")
    parser.add_argument("--out-summary", required=True, help="ABRicate-style summary table.")
    parser.add_argument("--min-identity", type=float, default=90.0)
    parser.add_argument("--min-coverage", type=float, default=80.0)
    parser.add_argument("--database-version", default="")
    args = parser.parse_args()

    blast_path = Path(args.blast)
    out_results = Path(args.out_results)
    out_summary = Path(args.out_summary)
    out_results.parent.mkdir(parents=True, exist_ok=True)
    out_summary.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    if blast_path.exists():
        with blast_path.open(newline="", errors="ignore") as handle:
            reader = csv.DictReader(handle, fieldnames=BLAST_FIELDS, delimiter="\t")
            for row in reader:
                try:
                    identity = float(row.get("pident") or 0)
                    aln_length = float(row.get("length") or 0)
                    subject_length = float(row.get("slen") or 0)
                except ValueError:
                    continue
                coverage = 0.0 if subject_length <= 0 else 100.0 * aln_length / subject_length
                if identity < args.min_identity or coverage < args.min_coverage:
                    continue
                qstart = int(float(row.get("qstart") or 0))
                qend = int(float(row.get("qend") or 0))
                start, end = sorted((qstart, qend))
                feature_id = clean_feature_id(row.get("sseqid", ""))
                rows.append(
                    {
                        "#FILE": args.sample_id,
                        "SEQUENCE": row.get("qseqid", ""),
                        "START": str(start),
                        "END": str(end),
                        "GENE": feature_id,
                        "COVERAGE": f"{int(aln_length)}/{int(subject_length)}" if subject_length else str(int(aln_length)),
                        "COVERAGE_MAP": "",
                        "GAPS": "",
                        "%COVERAGE": f"{coverage:.2f}",
                        "%IDENTITY": f"{identity:.2f}",
                        "DATABASE": "isfinder",
                        "ACCESSION": row.get("sseqid", ""),
                        "PRODUCT": "insertion sequence",
                        "RESISTANCE": "",
                    }
                )

    rows.sort(key=lambda item: (item["#FILE"], item["SEQUENCE"], int(item["START"]), item["GENE"]))

    with out_results.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ABRICATE_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    unique_features = sorted({row["GENE"] for row in rows})
    with out_summary.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["#FILE", "NUM_FOUND", *unique_features])
        writer.writerow([args.sample_id, len(unique_features), *("1" for _ in unique_features)])


if __name__ == "__main__":
    main()
