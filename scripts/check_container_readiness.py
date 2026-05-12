#!/usr/bin/env python3
"""Check whether an experimental container profile is ready to use."""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path


FIELDS = [
    "runtime",
    "runtime_available",
    "runtime_path",
    "image",
    "image_supplied",
    "database_paths_checked",
    "missing_database_paths",
    "status",
    "message",
]


RUNTIME_COMMANDS = {
    "docker": "docker",
    "apptainer": "apptainer",
    "singularity": "singularity",
}


def as_list(value: str) -> list[str]:
    return [item.strip() for item in str(value or "").split(",") if item.strip()]


def capture(command: list[str]) -> str:
    try:
        completed = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
            timeout=20,
        )
    except Exception as exc:  # pragma: no cover - defensive CLI helper
        return f"unavailable: {exc}"
    output = (completed.stdout or "").strip()
    return output.splitlines()[0] if output else "unavailable"


def write_report(path: Path, row: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime", choices=sorted(RUNTIME_COMMANDS), required=True)
    parser.add_argument("--image", default="", help="Container image passed to --container_image.")
    parser.add_argument(
        "--database-paths",
        default="",
        help="Comma-separated database paths that must exist on the host.",
    )
    parser.add_argument("--out", default="", help="Optional TSV report path.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    runtime_path = shutil.which(RUNTIME_COMMANDS[args.runtime])
    db_paths = as_list(args.database_paths)
    missing_paths = [path for path in db_paths if not Path(path).exists()]
    image = args.image.strip()

    if runtime_path and image and not missing_paths:
        status = "PASS"
        message = f"{args.runtime} is available: {capture([runtime_path, '--version'])}"
    elif not runtime_path:
        status = "FAIL_RUNTIME_MISSING"
        message = f"{args.runtime} executable was not found in PATH."
    elif not image:
        status = "FAIL_IMAGE_MISSING"
        message = "No container image was supplied. Pass --container_image to Nextflow."
    else:
        status = "FAIL_DATABASE_PATHS_MISSING"
        message = "One or more required database paths do not exist on the host."

    row = {
        "runtime": args.runtime,
        "runtime_available": str(runtime_path is not None).lower(),
        "runtime_path": runtime_path or "",
        "image": image,
        "image_supplied": str(bool(image)).lower(),
        "database_paths_checked": ",".join(db_paths),
        "missing_database_paths": ",".join(missing_paths),
        "status": status,
        "message": message,
    }
    if args.out:
        write_report(Path(args.out), row)
    print(f"status={status}")
    print(message)
    if args.out:
        print(f"report={args.out}")
    return 0 if status == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
