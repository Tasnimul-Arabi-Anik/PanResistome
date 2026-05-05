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
    args = parser.parse_args()
    sample_dir = Path(args.sample_dir)
    out_dir = sample_dir / "mash" / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)
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


if __name__ == "__main__":
    main()
