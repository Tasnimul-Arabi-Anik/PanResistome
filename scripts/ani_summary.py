#!/usr/bin/env python3
import argparse
import csv
import os
from pathlib import Path


def accession(value):
    text = Path(str(value)).name
    if text.endswith(".gz"):
        text = text[:-3]
    for suffix in [".fna", ".fa", ".fasta", ".fas"]:
        if text.lower().endswith(suffix):
            text = text[:-len(suffix)]
            break
    if text.endswith("_genomic"):
        text = text[:-8]
    return text


def read_pairs(path):
    rows = []
    if not path.exists():
        return rows
    with path.open(newline="") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                parts = line.split()
            if len(parts) < 3:
                continue
            try:
                ani = float(parts[2])
            except ValueError:
                continue
            rows.append({
                "query": accession(parts[0]),
                "reference": accession(parts[1]),
                "ani": ani,
                "fragments_mapped": parts[3] if len(parts) > 3 else "",
                "fragments_total": parts[4] if len(parts) > 4 else "",
            })
    return rows


def read_genomes(path):
    genomes = []
    if not path or not path.exists():
        return genomes
    with path.open() as handle:
        for line in handle:
            text = line.strip()
            if text:
                genomes.append(accession(text))
    return genomes


def union_find(items):
    parent = {item: item for item in items}

    def find(item):
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    def union(left, right):
        root_left = find(left)
        root_right = find(right)
        if root_left != root_right:
            parent[root_right] = root_left

    return parent, find, union


def main():
    parser = argparse.ArgumentParser(description="Summarize FastANI/skani pairwise ANI output for PanResistome.")
    parser.add_argument("--sample-dir", required=True)
    parser.add_argument("--pairs", required=True)
    parser.add_argument("--genomes-list")
    parser.add_argument("--tool", default="fastani")
    parser.add_argument("--duplicate-threshold", type=float, default=99.9)
    parser.add_argument("--species-threshold", type=float, default=95.0)
    args = parser.parse_args()

    sample_dir = Path(args.sample_dir)
    out_dir = sample_dir / "ani" / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = read_pairs(Path(args.pairs))
    listed_genomes = set(read_genomes(Path(args.genomes_list))) if args.genomes_list else set()
    genomes = sorted(listed_genomes | {row["query"] for row in rows} | {row["reference"] for row in rows})

    long_path = out_dir / "pairwise_ani_long.csv"
    with long_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["query", "reference", "ani", "fragments_mapped", "fragments_total", "tool"])
        writer.writeheader()
        for row in rows:
            out = dict(row)
            out["tool"] = args.tool
            writer.writerow(out)

    matrix = {genome: {other: "" for other in genomes} for genome in genomes}
    for genome in genomes:
        matrix[genome][genome] = "100.0"
    for row in rows:
        matrix[row["query"]][row["reference"]] = f"{row['ani']:.6g}"
    matrix_path = out_dir / "ani_matrix.csv"
    with matrix_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["genome"] + genomes)
        for genome in genomes:
            writer.writerow([genome] + [matrix[genome][other] for other in genomes])

    best = {}
    for row in rows:
        if row["query"] == row["reference"]:
            continue
        current = best.get(row["query"])
        if current is None or row["ani"] > current["closest_ani"]:
            best[row["query"]] = {"genome": row["query"], "closest_genome": row["reference"], "closest_ani": row["ani"]}
    closest_path = out_dir / "closest_genome.csv"
    with closest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["genome", "closest_genome", "closest_ani", "species_consistency_status"])
        writer.writeheader()
        for genome in genomes:
            row = best.get(genome, {"genome": genome, "closest_genome": "", "closest_ani": ""})
            ani = row["closest_ani"]
            row["species_consistency_status"] = "PASS" if ani != "" and float(ani) >= args.species_threshold else "WARN"
            writer.writerow(row)

    parent, find, union = union_find(genomes)
    for row in rows:
        if row["query"] != row["reference"] and row["ani"] >= args.duplicate_threshold:
            union(row["query"], row["reference"])
    clusters = {}
    for genome in genomes:
        clusters.setdefault(find(genome), []).append(genome)
    cluster_path = out_dir / "duplicate_clusters.csv"
    with cluster_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["ani_cluster", "representative", "genome", "cluster_size", "duplicate_threshold"])
        writer.writeheader()
        for idx, members in enumerate(sorted(clusters.values(), key=lambda m: (len(m), m[0]), reverse=True), start=1):
            representative = sorted(members)[0]
            for genome in sorted(members):
                writer.writerow({
                    "ani_cluster": f"ANI_CLUSTER_{idx:04d}",
                    "representative": representative,
                    "genome": genome,
                    "cluster_size": len(members),
                    "duplicate_threshold": args.duplicate_threshold,
                })

    outlier_path = out_dir / "ani_outliers.csv"
    with outlier_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["genome", "closest_genome", "closest_ani", "reason"])
        writer.writeheader()
        for genome in genomes:
            row = best.get(genome)
            if not row or row["closest_ani"] < args.species_threshold:
                writer.writerow({
                    "genome": genome,
                    "closest_genome": row["closest_genome"] if row else "",
                    "closest_ani": row["closest_ani"] if row else "",
                    "reason": f"closest_ani_below_{args.species_threshold:g}",
                })

    panr2_path = out_dir / "panr2_ani_summary.csv"
    with panr2_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "assembly_accession", "database", "feature_id", "feature_category", "presence", "identity", "coverage", "contig", "start", "end", "tool", "tool_version", "database_version"])
        writer.writeheader()
        cluster_by_genome = {}
        for row in csv.DictReader(cluster_path.open()):
            cluster_by_genome[row["genome"]] = row["ani_cluster"]
        for genome in genomes:
            writer.writerow({
                "sample_id": genome,
                "assembly_accession": genome,
                "database": "ani",
                "feature_id": cluster_by_genome.get(genome, "ANI_CLUSTER_UNASSIGNED"),
                "feature_category": "ani_cluster",
                "presence": 1,
                "identity": best.get(genome, {}).get("closest_ani", ""),
                "coverage": "",
                "contig": "",
                "start": "",
                "end": "",
                "tool": args.tool,
                "tool_version": "",
                "database_version": "",
            })


if __name__ == "__main__":
    main()
