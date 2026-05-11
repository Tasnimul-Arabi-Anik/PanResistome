#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


def accession(value):
    text = Path(str(value)).name
    for suffix in [".fna", ".fa", ".fasta", ".fas"]:
        if text.lower().endswith(suffix):
            text = text[:-len(suffix)]
            break
    return text


def main():
    parser = argparse.ArgumentParser(description="Summarize Mash distance output.")
    parser.add_argument("--sample-dir", required=True)
    parser.add_argument("--dist", required=True)
    parser.add_argument("--genomes-list", default="")
    args = parser.parse_args()
    sample_dir = Path(args.sample_dir)
    out_dir = sample_dir / "mash" / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)
    genomes = []
    genomes_list = Path(args.genomes_list) if args.genomes_list else None
    if genomes_list and genomes_list.exists():
        genomes = [line.strip() for line in genomes_list.read_text(encoding="utf-8").splitlines() if line.strip()]
    rows = []
    dist_path = Path(args.dist)
    if dist_path.exists():
        with dist_path.open() as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                rows.append({
                    "reference": accession(parts[0]),
                    "query": accession(parts[1]),
                    "mash_distance": parts[2],
                    "p_value": parts[3],
                    "matching_hashes": parts[4],
                })
    long_path = out_dir / "mash_distance_long.csv"
    with long_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["query", "reference", "mash_distance", "p_value", "matching_hashes"])
        writer.writeheader()
        writer.writerows(rows)
    best = {}
    for row in rows:
        if row["query"] == row["reference"]:
            continue
        distance = float(row["mash_distance"])
        current = best.get(row["query"])
        if current is None or distance < float(current["mash_distance"]):
            best[row["query"]] = row
    closest_path = out_dir / "closest_mash_neighbor.csv"
    with closest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["query", "reference", "mash_distance", "p_value", "matching_hashes"])
        writer.writeheader()
        writer.writerows(best.values())

    genome_count = len(genomes)
    nonself_rows = [row for row in rows if row["query"] != row["reference"]]
    if genome_count < 2:
        decision = "insufficient_genomes"
        status = "SKIPPED_INAPPLICABLE"
        message = "Mash requires at least two genomes; pairwise Mash screening was skipped."
    elif not nonself_rows:
        decision = "no_pairwise_hits"
        status = "WARNING_EMPTY"
        message = "Mash ran or was requested for at least two genomes but no non-self pairwise hits were found."
    else:
        decision = "summarized"
        status = "PASS"
        message = "Mash pairwise distances summarized."
    status_path = out_dir / "mash_run_status.tsv"
    with status_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["tool", "genome_count", "pair_rows", "nonself_pair_rows", "decision", "status", "message"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "tool": "mash",
                "genome_count": genome_count,
                "pair_rows": len(rows),
                "nonself_pair_rows": len(nonself_rows),
                "decision": decision,
                "status": status,
                "message": message,
            }
        )


if __name__ == "__main__":
    main()
