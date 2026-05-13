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
    "pull_test_requested",
    "pull_test_timeout_seconds",
    "pull_test_passed",
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


def image_exec_command(runtime: str, runtime_path: str, image: str) -> list[str]:
    if runtime == "docker":
        return [runtime_path, "run", "--rm", "--entrypoint", "true", image]
    return [runtime_path, "exec", image, "true"]


def run_image_pull_test(
    runtime: str,
    runtime_path: str,
    image: str,
    timeout_seconds: int,
) -> tuple[bool, str]:
    command = image_exec_command(runtime, runtime_path, image)
    try:
        completed = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
            timeout=timeout_seconds,
        )
    except Exception as exc:  # pragma: no cover - defensive CLI helper
        return False, f"container pull/exec test could not run: {exc}"
    output = (completed.stdout or "").strip()
    if completed.returncode == 0:
        return True, "container pull/exec test passed."
    detail = output.splitlines()[-1] if output else f"exit_status={completed.returncode}"
    return False, f"container pull/exec test failed: {detail}"


def write_report(path: Path, row: dict[str, str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runtime", choices=sorted(RUNTIME_COMMANDS), required=True)
    parser.add_argument("--image", default="", help="Container image passed to --container_image.")
    parser.add_argument(
        "--database-paths",
        default="",
        help="Comma-separated database paths that must exist on the host.",
    )
    parser.add_argument(
        "--pull-test",
        action="store_true",
        help="Run a small container pull/exec test with the requested runtime and image.",
    )
    parser.add_argument(
        "--pull-test-timeout",
        type=int,
        default=300,
        help=(
            "Seconds to allow for --pull-test. Large Singularity/Apptainer "
            "image conversion can require much longer than the default."
        ),
    )
    parser.add_argument("--out", default="", help="Optional TSV report path.")
    return parser.parse_args(argv)


def main() -> int:
    args = parse_args()
    runtime_path = shutil.which(RUNTIME_COMMANDS[args.runtime])
    db_paths = as_list(args.database_paths)
    missing_paths = [path for path in db_paths if not Path(path).exists()]
    image = args.image.strip()
    pull_test_passed = ""

    if runtime_path and image and not missing_paths and args.pull_test:
        pull_ok, pull_message = run_image_pull_test(
            args.runtime,
            runtime_path,
            image,
            max(1, int(args.pull_test_timeout)),
        )
        pull_test_passed = str(pull_ok).lower()
        if pull_ok:
            status = "PASS"
            message = f"{args.runtime} is available: {capture([runtime_path, '--version'])}; {pull_message}"
        else:
            status = "FAIL_PULL_TEST"
            message = pull_message
    elif runtime_path and image and not missing_paths:
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
        "pull_test_requested": str(bool(args.pull_test)).lower(),
        "pull_test_timeout_seconds": str(max(1, int(args.pull_test_timeout))),
        "pull_test_passed": pull_test_passed,
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
