import csv
import os
import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def write_executable(path: Path, text: str) -> None:
    path.write_text(textwrap.dedent(text).lstrip(), encoding="utf-8")
    path.chmod(0o755)


def read_status(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class OptionalParallelRunnerTests(unittest.TestCase):
    def test_mobsuite_parallel_runner_writes_per_sample_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bin_dir = root / "bin"
            bin_dir.mkdir()
            write_executable(
                bin_dir / "mob_recon",
                """
                #!/usr/bin/env python3
                import argparse
                from pathlib import Path
                parser = argparse.ArgumentParser()
                parser.add_argument("--infile")
                parser.add_argument("--outdir")
                parser.add_argument("--num_threads")
                parser.add_argument("--force", action="store_true")
                parser.add_argument("--database_directory", default="")
                args = parser.parse_args()
                out = Path(args.outdir)
                out.mkdir(parents=True, exist_ok=True)
                sample = Path(args.infile).stem
                (out / "contig_report.tsv").write_text(
                    "sample_id\\trep_type(s)\\tpredicted_mobility\\n"
                    f"{sample}\\tIncFIB\\tconjugative\\n",
                    encoding="utf-8",
                )
                """,
            )
            seq = root / "sequence"
            seq.mkdir()
            for sample in ["GCF_000001", "GCF_000002", "GCF_000003"]:
                (seq / f"{sample}.fna").write_text(">contig\nATGC\n", encoding="utf-8")
            status = root / "mobsuite_sample_status.tsv"
            env = os.environ.copy()
            env["PATH"] = f"{bin_dir}:{env['PATH']}"
            subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "run_mobsuite_parallel.py"),
                    "--sequence-dir",
                    str(seq),
                    "--raw-dir",
                    str(root / "raw"),
                    "--status-file",
                    str(status),
                    "--jobs",
                    "2",
                    "--threads-per-sample",
                    "1",
                ],
                env=env,
                check=True,
            )
            rows = read_status(status)
            self.assertEqual(len(rows), 3)
            self.assertEqual({row["status"] for row in rows}, {"PASS"})
            self.assertTrue((root / "raw" / "GCF_000001" / "contig_report.tsv").exists())

    def test_kleborate_parallel_runner_writes_per_sample_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            bin_dir = root / "bin"
            bin_dir.mkdir()
            write_executable(
                bin_dir / "kleborate",
                """
                #!/usr/bin/env python3
                import argparse
                from pathlib import Path
                parser = argparse.ArgumentParser()
                parser.add_argument("-a", "--assemblies", nargs="+")
                parser.add_argument("-o", "--output")
                parser.add_argument("--preset", default="")
                args = parser.parse_args()
                out = Path(args.output)
                out.mkdir(parents=True, exist_ok=True)
                sample = Path(args.assemblies[0]).stem
                (out / "kleborate.tsv").write_text(
                    "sample_id\\tST\\tyersiniabactin\\n"
                    f"{sample}\\tST11\\tybt_9\\n",
                    encoding="utf-8",
                )
                """,
            )
            seq = root / "sequence"
            seq.mkdir()
            for sample in ["GCF_000001", "GCF_000002"]:
                (seq / f"{sample}.fna").write_text(">contig\nATGC\n", encoding="utf-8")
            status = root / "kleborate_sample_status.tsv"
            env = os.environ.copy()
            env["PATH"] = f"{bin_dir}:{env['PATH']}"
            subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "run_kleborate_parallel.py"),
                    "--sequence-dir",
                    str(seq),
                    "--raw-dir",
                    str(root / "raw"),
                    "--status-file",
                    str(status),
                    "--jobs",
                    "2",
                ],
                env=env,
                check=True,
            )
            rows = read_status(status)
            self.assertEqual(len(rows), 2)
            self.assertEqual({row["status"] for row in rows}, {"PASS"})
            self.assertTrue((root / "raw" / "GCF_000001" / "kleborate.tsv").exists())


if __name__ == "__main__":
    unittest.main()
