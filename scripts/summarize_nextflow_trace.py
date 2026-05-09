#!/usr/bin/env python3
"""Summarize Nextflow trace output into compact runtime/resource tables."""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


TIME_UNITS = {
    "ms": 0.001,
    "s": 1.0,
    "sec": 1.0,
    "secs": 1.0,
    "m": 60.0,
    "min": 60.0,
    "mins": 60.0,
    "h": 3600.0,
    "hr": 3600.0,
    "hrs": 3600.0,
    "d": 86400.0,
}

MEM_UNITS = {
    "b": 1.0,
    "kb": 1024.0,
    "kib": 1024.0,
    "mb": 1024.0**2,
    "mib": 1024.0**2,
    "gb": 1024.0**3,
    "gib": 1024.0**3,
    "tb": 1024.0**4,
    "tib": 1024.0**4,
}


def parse_duration_seconds(value: str) -> float:
    text = str(value or "").strip().lower()
    if not text or text == "-":
        return 0.0
    total = 0.0
    for number, unit in re.findall(r"([0-9]*\.?[0-9]+)\s*([a-z]+)", text):
        total += float(number) * TIME_UNITS.get(unit, 0.0)
    if total:
        return total
    try:
        return float(text)
    except ValueError:
        return 0.0


def parse_memory_bytes(value: str) -> float:
    text = str(value or "").strip().lower().replace(",", "")
    if not text or text == "-":
        return 0.0
    match = re.match(r"^([0-9]*\.?[0-9]+)\s*([a-z]+)?$", text)
    if not match:
        return 0.0
    number = float(match.group(1))
    unit = match.group(2) or "b"
    return number * MEM_UNITS.get(unit, 1.0)


def process_name(task_name: str) -> str:
    return re.sub(r"\s+\(\d+\)$", "", str(task_name or "").strip())


def fmt_seconds(seconds: float) -> str:
    if seconds >= 3600:
        return f"{seconds / 3600:.2f}h"
    if seconds >= 60:
        return f"{seconds / 60:.2f}m"
    return f"{seconds:.2f}s"


def fmt_gib(bytes_value: float) -> str:
    return f"{bytes_value / (1024.0**3):.3f}"


def first_existing(row: dict[str, str], names: list[str]) -> str:
    for name in names:
        if name in row and row[name]:
            return row[name]
    return ""


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize a Nextflow trace TSV.")
    parser.add_argument("--trace", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--tasks-out")
    args = parser.parse_args()

    trace = Path(args.trace)
    out = Path(args.out)
    tasks_out = Path(args.tasks_out) if args.tasks_out else out.with_name(out.stem + "_tasks.tsv")
    out.parent.mkdir(parents=True, exist_ok=True)

    if not trace.exists() or trace.stat().st_size == 0:
        out.write_text("process\ttasks\tstatus_counts\ttotal_realtime\tmax_realtime\tmax_peak_rss_gib\tmax_peak_vmem_gib\tmean_cpu_percent\n", encoding="utf-8")
        tasks_out.write_text("process\ttask_name\tstatus\trealtime_seconds\tpeak_rss_gib\tpeak_vmem_gib\tcpu_percent\n", encoding="utf-8")
        return 0

    with trace.open(newline="", encoding="utf-8", errors="ignore") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    task_rows = []
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        task = first_existing(row, ["name", "task", "process"])
        proc = process_name(task)
        status = first_existing(row, ["status"])
        realtime_seconds = parse_duration_seconds(first_existing(row, ["realtime", "duration", "time"]))
        peak_rss = parse_memory_bytes(first_existing(row, ["peak_rss", "rss"]))
        peak_vmem = parse_memory_bytes(first_existing(row, ["peak_vmem", "vmem"]))
        cpu_text = first_existing(row, ["%cpu", "pcpu", "cpu"]).replace("%", "")
        try:
            cpu = float(cpu_text) if cpu_text else 0.0
        except ValueError:
            cpu = 0.0
        parsed = {
            "process": proc,
            "task_name": task,
            "status": status,
            "realtime_seconds": realtime_seconds,
            "peak_rss": peak_rss,
            "peak_vmem": peak_vmem,
            "cpu": cpu,
        }
        grouped[proc].append(parsed)
        task_rows.append({
            "process": proc,
            "task_name": task,
            "status": status,
            "realtime_seconds": f"{realtime_seconds:.3f}",
            "peak_rss_gib": fmt_gib(peak_rss),
            "peak_vmem_gib": fmt_gib(peak_vmem),
            "cpu_percent": f"{cpu:.2f}",
        })

    with tasks_out.open("w", newline="", encoding="utf-8") as handle:
        fields = ["process", "task_name", "status", "realtime_seconds", "peak_rss_gib", "peak_vmem_gib", "cpu_percent"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(task_rows, key=lambda row: (row["process"], row["task_name"])))

    summary_rows = []
    for proc, items in sorted(grouped.items()):
        status_counts: dict[str, int] = defaultdict(int)
        for item in items:
            status_counts[str(item["status"])] += 1
        total_realtime = sum(float(item["realtime_seconds"]) for item in items)
        max_realtime = max(float(item["realtime_seconds"]) for item in items)
        max_peak_rss = max(float(item["peak_rss"]) for item in items)
        max_peak_vmem = max(float(item["peak_vmem"]) for item in items)
        mean_cpu = sum(float(item["cpu"]) for item in items) / len(items)
        summary_rows.append({
            "process": proc,
            "tasks": str(len(items)),
            "status_counts": ";".join(f"{key}:{value}" for key, value in sorted(status_counts.items())),
            "total_realtime": fmt_seconds(total_realtime),
            "max_realtime": fmt_seconds(max_realtime),
            "max_peak_rss_gib": fmt_gib(max_peak_rss),
            "max_peak_vmem_gib": fmt_gib(max_peak_vmem),
            "mean_cpu_percent": f"{mean_cpu:.2f}",
        })

    with out.open("w", newline="", encoding="utf-8") as handle:
        fields = ["process", "tasks", "status_counts", "total_realtime", "max_realtime", "max_peak_rss_gib", "max_peak_vmem_gib", "mean_cpu_percent"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
