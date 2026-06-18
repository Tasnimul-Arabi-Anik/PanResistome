#!/usr/bin/env python3
"""Conservative post-run cleanup for opt-in PanResistome runs."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


FIELDS = ["target", "path", "action", "status", "bytes_removed", "message"]


def _is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
        return True
    except ValueError:
        return False


def _dir_size(path: Path) -> int:
    total = 0
    if not path.exists():
        return total
    if path.is_file() or path.is_symlink():
        try:
            return path.stat().st_size
        except OSError:
            return 0
    for child in path.rglob("*"):
        try:
            if child.is_file() or child.is_symlink():
                total += child.stat().st_size
        except OSError:
            continue
    return total


def _safe_to_delete(target: Path, protected: list[Path]) -> tuple[bool, str]:
    target = target.resolve()
    protected_roots = {"", "/", str(Path.home().resolve()), "/tmp", "/var/tmp"}
    if str(target) in protected_roots:
        return False, "refusing to remove root, empty path, home directory, or shared temp root"
    for protected_path in protected:
        protected_path = protected_path.resolve()
        if target == protected_path:
            return False, f"target is protected path {protected_path}"
        if _is_relative_to(protected_path, target):
            return False, f"target contains protected path {protected_path}"
    return True, ""


def _remove_path(target: Path, label: str, protected: list[Path], dry_run: bool) -> dict[str, str]:
    row = {
        "target": label,
        "path": str(target),
        "action": "remove",
        "status": "SKIPPED",
        "bytes_removed": "0",
        "message": "",
    }
    if not target.exists():
        row["message"] = "path did not exist"
        return row
    safe, reason = _safe_to_delete(target, protected)
    if not safe:
        row["status"] = "SKIPPED_UNSAFE"
        row["message"] = reason
        return row
    bytes_before = _dir_size(target)
    if dry_run:
        row["status"] = "DRY_RUN"
        row["bytes_removed"] = str(bytes_before)
        row["message"] = "would remove path"
        return row
    try:
        if target.is_dir() and not target.is_symlink():
            shutil.rmtree(target)
        else:
            target.unlink()
        row["status"] = "PASS"
        row["bytes_removed"] = str(bytes_before)
        row["message"] = "removed"
    except PermissionError as error:
        row["status"] = "WARNING_PERMISSION_DENIED"
        row["message"] = f"permission denied: {error}"
    except OSError as error:
        row["status"] = "WARNING_FAILED"
        row["message"] = str(error)
    return row


def write_rows(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("\t".join(FIELDS) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(field, "")) for field in FIELDS) + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, required=True, help="Nextflow work directory to remove")
    parser.add_argument("--launch-dir", type=Path, required=True, help="Nextflow launch directory")
    parser.add_argument("--outdir", type=Path, required=True, help="Published output directory")
    parser.add_argument("--summary", type=Path, required=True, help="Cleanup audit TSV to write")
    parser.add_argument("--nextflow-session-id", default="", help="Optional Nextflow session id for cache cleanup")
    parser.add_argument("--dry-run", action="store_true", help="Write audit only; do not remove files")
    args = parser.parse_args(argv)

    launch_dir = args.launch_dir.resolve()
    outdir = args.outdir if args.outdir.is_absolute() else launch_dir / args.outdir
    protected = [launch_dir, outdir.resolve()]
    rows = [
        _remove_path(args.work_dir, "nextflow_work_dir", protected, args.dry_run),
    ]

    session_id = str(args.nextflow_session_id or "").strip()
    if session_id and "/" not in session_id and "\\" not in session_id:
        rows.append(
            _remove_path(
                launch_dir / ".nextflow" / "cache" / session_id,
                "nextflow_session_cache",
                protected,
                args.dry_run,
            )
        )

    for cache_path in [
        launch_dir / ".pytest_cache",
        launch_dir / "scripts" / "__pycache__",
        launch_dir / "tests" / "__pycache__",
    ]:
        rows.append(_remove_path(cache_path, "transient_cache", protected, args.dry_run))

    write_rows(args.summary, rows)
    failed = [row for row in rows if row["status"].startswith("WARNING")]
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
