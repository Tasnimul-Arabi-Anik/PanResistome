#!/usr/bin/env python3
import argparse
import csv
from pathlib import Path


METRICS = {
    "# contigs": "quast_num_contigs",
    "Total length": "quast_total_length",
    "Largest contig": "quast_largest_contig",
    "N50": "quast_n50",
    "GC (%)": "quast_gc_percent",
    "# N's per 100 kbp": "quast_ns_per_100kbp",
}


def accession(name):
    text = Path(str(name)).name
    for suffix in [".fna", ".fa", ".fasta", ".fas"]:
        if text.lower().endswith(suffix):
            return text[:-len(suffix)]
    return Path(text).stem


def main():
    parser = argparse.ArgumentParser(description="Convert QUAST report.tsv to PanResistome assembly QC tables.")
    parser.add_argument("--sample-dir", required=True)
    args = parser.parse_args()
    sample_dir = Path(args.sample_dir)
    report = sample_dir / "quast" / "report.tsv"
    out_dir = sample_dir / "quast" / "analysis"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "assembly_qc.csv"
    panr2_path = out_dir / "panr2_quast_summary.csv"

    rows = []
    if report.exists():
        with report.open(newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            header = next(reader, [])
            assemblies = header[1:]
            by_assembly = {assembly: {"assembly_accession": accession(assembly), "quast_assembly": assembly} for assembly in assemblies}
            for row in reader:
                if not row:
                    continue
                metric = row[0]
                field = METRICS.get(metric)
                if not field:
                    continue
                for assembly, value in zip(assemblies, row[1:]):
                    by_assembly[assembly][field] = value
            rows = list(by_assembly.values())

    fields = ["assembly_accession", "quast_assembly"] + sorted(METRICS.values())
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    with panr2_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["sample_id", "assembly_accession", "database", "feature_id", "feature_category", "presence", "identity", "coverage", "contig", "start", "end", "tool", "tool_version", "database_version"])
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "sample_id": row.get("assembly_accession", ""),
                "assembly_accession": row.get("assembly_accession", ""),
                "database": "assembly_qc",
                "feature_id": "QUAST_METRICS",
                "feature_category": "assembly_structure",
                "presence": 1,
                "identity": row.get("quast_n50", ""),
                "coverage": row.get("quast_total_length", ""),
                "contig": row.get("quast_num_contigs", ""),
                "start": "",
                "end": "",
                "tool": "quast",
                "tool_version": "",
                "database_version": "",
            })


if __name__ == "__main__":
    main()
