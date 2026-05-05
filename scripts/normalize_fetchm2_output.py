#!/usr/bin/env python3
"""Adapt FetchM2 output into the existing PanResistome sample layout."""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

import pandas as pd


COMPATIBILITY_COLUMNS = {
    "Geographic Location": ["Country"],
    "Collection Date": ["Collection_Year", "Collection Date"],
    "Host": ["Host_SD", "Host_Cleaned", "Host_Original"],
    "Isolation Source": ["Isolation_Source_SD", "Sample_Type_SD", "Environment_Medium_SD"],
}


def _clean_name(value: str) -> str:
    text = str(value or "").strip()
    if not text:
        return "fetchm2_dataset"
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", text).strip("_")
    return cleaned or "fetchm2_dataset"


def _derive_dataset_name(metadata: pd.DataFrame, fallback: str) -> str:
    for column in ["Organism Name", "Species", "Genus"]:
        if column not in metadata.columns:
            continue
        values = [
            str(value).strip()
            for value in metadata[column].dropna().astype(str).unique()
            if str(value).strip() and str(value).strip().lower() not in {"nan", "none", "0", "unknown"}
        ]
        if len(values) == 1:
            return _clean_name(values[0])
    return _clean_name(fallback)


def _fill_compatibility_columns(metadata: pd.DataFrame) -> pd.DataFrame:
    df = metadata.copy()
    for target, sources in COMPATIBILITY_COLUMNS.items():
        if target not in df.columns:
            df[target] = ""
        target_values = df[target].fillna("").astype(str).str.strip()
        for source in sources:
            if source not in df.columns:
                continue
            source_values = df[source].fillna("").astype(str).str.strip()
            mask = target_values.eq("") | target_values.str.lower().isin({"nan", "none", "0"})
            if mask.any():
                df.loc[mask, target] = source_values.loc[mask]
                target_values = df[target].fillna("").astype(str).str.strip()
    if "Collection Date" in df.columns:
        df["Collection Date"] = df["Collection Date"].astype(str).str.extract(r"(\d{4})", expand=False).fillna(df["Collection Date"].astype(str))
    if "Organism Taxonomic ID" in df.columns:
        if "TaxID" not in df.columns:
            df["TaxID"] = df["Organism Taxonomic ID"]
        else:
            missing_taxid = df["TaxID"].fillna("").astype(str).str.strip().eq("")
            df.loc[missing_taxid, "TaxID"] = df.loc[missing_taxid, "Organism Taxonomic ID"]
    if "Organism Name" in df.columns:
        organism = df["Organism Name"].fillna("").astype(str).str.strip()
        genus = organism.str.split().str[0].replace({"nan": "", "None": "", "0": ""})
        species = organism.str.split().str[:2].str.join(" ").replace({"nan": "", "None": "", "0": ""})
        if "Genus" not in df.columns:
            df["Genus"] = genus
        else:
            missing_genus = df["Genus"].fillna("").astype(str).str.strip().eq("")
            df.loc[missing_genus, "Genus"] = genus.loc[missing_genus]
        if "Species" not in df.columns:
            df["Species"] = species
        else:
            missing_species = df["Species"].fillna("").astype(str).str.strip().eq("")
            df.loc[missing_species, "Species"] = species.loc[missing_species]
    return df


def _move_child(root: Path, child_name: str, sample_dir: Path) -> None:
    source = root / child_name
    if not source.exists() or source.resolve() == sample_dir.resolve():
        return
    target = sample_dir / child_name
    if target.exists():
        return
    shutil.move(str(source), str(target))


def normalize_fetchm2_output(results_dir: Path, dataset_name: str | None = None) -> Path:
    clean_path = results_dir / "metadata_output" / "fetchm2_clean.csv"
    if not clean_path.exists():
        raise FileNotFoundError(f"FetchM2 clean metadata not found: {clean_path}")

    metadata = pd.read_csv(clean_path)
    sample_name = _derive_dataset_name(metadata, dataset_name or results_dir.name)
    sample_dir = results_dir / sample_name
    sample_dir.mkdir(parents=True, exist_ok=True)

    for child in ["metadata_output", "metadata_analysis", "audit", "sequence", "validation"]:
        _move_child(results_dir, child, sample_dir)

    metadata_dir = sample_dir / "metadata_output"
    metadata_dir.mkdir(parents=True, exist_ok=True)
    fetchm2_clean = metadata_dir / "fetchm2_clean.csv"
    if not fetchm2_clean.exists():
        raise FileNotFoundError(f"FetchM2 metadata was not moved into sample directory: {fetchm2_clean}")

    compat = _fill_compatibility_columns(pd.read_csv(fetchm2_clean))
    compat.to_csv(metadata_dir / "ncbi_clean.csv", index=False)
    compat.to_csv(metadata_dir / "fetchm2_clean_compat.csv", index=False)

    fetchm2_tsv = metadata_dir / "fetchm2_clean.tsv"
    if fetchm2_tsv.exists():
        compat.to_csv(metadata_dir / "ncbi_clean.tsv", sep="\t", index=False)

    manifest_dir = sample_dir / "metadata_output"
    (manifest_dir / "metadata_engine.txt").write_text("fetchm2\n", encoding="utf-8")
    return sample_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Normalize FetchM2 output into PanResistome sample directory layout.")
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--dataset-name", help="Optional PanResistome sample directory name.")
    args = parser.parse_args()
    sample_dir = normalize_fetchm2_output(args.results_dir, dataset_name=args.dataset_name)
    print(f"FetchM2 output normalized into {sample_dir}")


if __name__ == "__main__":
    main()
