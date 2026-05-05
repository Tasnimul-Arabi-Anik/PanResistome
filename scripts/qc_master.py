#!/usr/bin/env python3
import argparse
import csv
import shutil
from pathlib import Path


def norm(value):
    return str(value or "").strip()


def key_candidates(row):
    values = [
        row.get("Assembly Accession", ""),
        row.get("assembly_accession", ""),
        row.get("sequence_accession", ""),
        row.get("sequence_file", ""),
    ]
    keys = set()
    for value in values:
        text = norm(value)
        if not text:
            continue
        stem = Path(text).stem
        keys.add(text)
        keys.add(stem)
        if stem.endswith("_genomic"):
            keys.add(stem[:-8])
    return {key.lower().replace(".", "_").replace("-", "_") for key in keys if key}


def read_csv(path):
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_list(path, rows, field):
    with path.open("w") as handle:
        for row in rows:
            value = row.get(field) or row.get("Assembly Accession") or row.get("sequence_accession") or row.get("sequence_file")
            if value:
                handle.write(str(value) + "\n")


def main():
    parser = argparse.ArgumentParser(description="Build PanResistome combined QC master report.")
    parser.add_argument("--sample-dir", required=True)
    parser.add_argument("--qc-filter", default="false")
    parser.add_argument("--representative-only", default="false")
    parser.add_argument("--max-contigs", type=float)
    parser.add_argument("--min-n50", type=float)
    parser.add_argument("--min-completeness", type=float)
    parser.add_argument("--max-contamination", type=float)
    parser.add_argument("--ani-species-threshold", type=float, default=95.0)
    args = parser.parse_args()

    sample_dir = Path(args.sample_dir)
    qc_dir = sample_dir / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir = sample_dir / "metadata_output"
    sequence_dir = sample_dir / "sequence"
    filtered_dir = sample_dir / "sequence_filtered"

    metadata_path = metadata_dir / "ncbi_enriched.csv"
    if not metadata_path.exists():
        metadata_path = metadata_dir / "ncbi_clean.csv"
    rows = read_csv(metadata_path)
    if not rows:
        rows = read_csv(sample_dir / "quast" / "analysis" / "assembly_qc.csv")

    quast = {row.get("assembly_accession", "").lower().replace(".", "_").replace("-", "_"): row for row in read_csv(sample_dir / "quast" / "analysis" / "assembly_qc.csv")}
    ani_closest = {row.get("genome", "").lower().replace(".", "_").replace("-", "_"): row for row in read_csv(sample_dir / "ani" / "analysis" / "closest_genome.csv")}
    ani_cluster = {}
    representatives = set()
    for row in read_csv(sample_dir / "ani" / "analysis" / "duplicate_clusters.csv"):
        key = row.get("genome", "").lower().replace(".", "_").replace("-", "_")
        ani_cluster[key] = row
        if row.get("genome") == row.get("representative"):
            representatives.add(key)

    output_rows = []
    for row in rows:
        enriched = dict(row)
        keys = key_candidates(row)
        qrow = next((quast[key] for key in keys if key in quast), {})
        arow = next((ani_closest[key] for key in keys if key in ani_closest), {})
        crow = next((ani_cluster[key] for key in keys if key in ani_cluster), {})
        for key, value in qrow.items():
            enriched.setdefault(key, value)
        if arow:
            enriched["ani_closest_genome"] = arow.get("closest_genome", "")
            enriched["ani_closest_ani"] = arow.get("closest_ani", "")
            enriched["ani_species_consistency_status"] = arow.get("species_consistency_status", "")
        if crow:
            enriched["ani_cluster"] = crow.get("ani_cluster", "")
            enriched["ani_cluster_representative"] = crow.get("representative", "")
            enriched["ani_cluster_size"] = crow.get("cluster_size", "")

        reasons = []
        warning_reasons = []
        for prefix in ["sequence", "checkm2", "gtdbtk"]:
            status = enriched.get(f"{prefix}_qc_status")
            if status == "FAIL":
                reasons.append(enriched.get(f"{prefix}_qc_fail_reasons") or f"{prefix.upper()}_QC_FAIL")

        def number(field):
            try:
                value = str(enriched.get(field, "")).replace(",", "")
                return float(value) if value != "" else None
            except ValueError:
                return None

        if args.max_contigs is not None:
            contigs = number("quast_num_contigs") or number("sequence_num_contigs")
            if contigs is not None and contigs > args.max_contigs:
                reasons.append(f"max_contigs:{contigs:g}>{args.max_contigs:g}")
        if args.min_n50 is not None:
            n50 = number("quast_n50") or number("sequence_n50")
            if n50 is not None and n50 < args.min_n50:
                reasons.append(f"min_n50:{n50:g}<{args.min_n50:g}")
        if args.min_completeness is not None:
            completeness = number("checkm2_completeness")
            if completeness is not None and completeness < args.min_completeness:
                reasons.append(f"min_completeness:{completeness:g}<{args.min_completeness:g}")
        if args.max_contamination is not None:
            contamination = number("checkm2_contamination")
            if contamination is not None and contamination > args.max_contamination:
                reasons.append(f"max_contamination:{contamination:g}>{args.max_contamination:g}")

        closest_ani = number("ani_closest_ani")
        if closest_ani is not None and closest_ani < args.ani_species_threshold:
            warning_reasons.append(f"ani_species_warning:{closest_ani:g}<{args.ani_species_threshold:g}")
        if args.representative_only.lower() in {"true", "1", "yes"} and crow:
            genome_key = crow.get("genome", "").lower().replace(".", "_").replace("-", "_")
            if genome_key not in representatives:
                reasons.append("non_representative_duplicate_cluster_member")

        enriched["qc_master_status"] = "FAIL" if reasons else "WARN" if warning_reasons else "PASS"
        enriched["qc_master_fail_reasons"] = ";".join(reasons)
        enriched["qc_master_warning_reasons"] = ";".join(warning_reasons)
        output_rows.append(enriched)

    fieldnames = []
    for row in output_rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    master_path = qc_dir / "qc_master_report.csv"
    with master_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(output_rows)

    pass_rows = [row for row in output_rows if row.get("qc_master_status") == "PASS"]
    fail_rows = [row for row in output_rows if row.get("qc_master_status") == "FAIL"]
    warn_rows = [row for row in output_rows if row.get("qc_master_status") == "WARN"]
    write_list(qc_dir / "qc_pass_samples.txt", pass_rows, "sequence_file")
    write_list(qc_dir / "qc_fail_samples.txt", fail_rows, "sequence_file")
    write_list(qc_dir / "qc_warning_samples.txt", warn_rows, "sequence_file")
    with (qc_dir / "excluded_for_panr2.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(fail_rows)

    if args.qc_filter.lower() in {"true", "1", "yes"}:
        filtered_dir.mkdir(parents=True, exist_ok=True)
        for old in filtered_dir.glob("*.fna"):
            old.unlink()
        for row in pass_rows:
            sequence_file = row.get("sequence_file")
            if not sequence_file:
                continue
            source = sequence_dir / sequence_file
            if source.exists():
                shutil.copy2(source, filtered_dir / sequence_file)
        clean_path = metadata_dir / "ncbi_clean.csv"
        if clean_path.exists() and not (metadata_dir / "ncbi_clean_unfiltered.csv").exists():
            shutil.copy2(clean_path, metadata_dir / "ncbi_clean_unfiltered.csv")
        if fieldnames:
            with clean_path.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(pass_rows)


if __name__ == "__main__":
    main()
