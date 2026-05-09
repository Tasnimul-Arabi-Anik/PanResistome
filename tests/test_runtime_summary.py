import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]


class RuntimeSummaryTests(unittest.TestCase):
    def test_summarizes_nextflow_trace_by_process(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            trace = tmp / "trace.tsv"
            summary = tmp / "summary.tsv"
            tasks = tmp / "tasks.tsv"
            trace.write_text(
                "name\tstatus\trealtime\tpeak_rss\tpeak_vmem\t%cpu\n"
                "AMRFINDERPLUS_ANALYSIS\tCOMPLETED\t2m 30s\t1.5 GB\t2 GB\t180%\n"
                "AMRFINDERPLUS_ANALYSIS (2)\tCOMPLETED\t30s\t512 MB\t1 GB\t90%\n"
                "CHECKM2_QC\tFAILED\t1h\t4 GB\t5 GB\t75%\n",
                encoding="utf-8",
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "summarize_nextflow_trace.py"),
                    "--trace",
                    str(trace),
                    "--out",
                    str(summary),
                    "--tasks-out",
                    str(tasks),
                ],
                check=True,
            )

            summary_df = pd.read_csv(summary, sep="\t", dtype=str)
            tasks_df = pd.read_csv(tasks, sep="\t", dtype=str)
            amr = summary_df.loc[summary_df["process"] == "AMRFINDERPLUS_ANALYSIS"].iloc[0]
            checkm2 = summary_df.loc[summary_df["process"] == "CHECKM2_QC"].iloc[0]

            self.assertEqual(amr["tasks"], "2")
            self.assertEqual(amr["status_counts"], "COMPLETED:2")
            self.assertEqual(amr["total_realtime"], "3.00m")
            self.assertEqual(amr["max_realtime"], "2.50m")
            self.assertEqual(checkm2["status_counts"], "FAILED:1")
            self.assertEqual(checkm2["max_peak_rss_gib"], "4.000")
            self.assertEqual(set(tasks_df["process"]), {"AMRFINDERPLUS_ANALYSIS", "CHECKM2_QC"})


class AMRFinderPlusParallelRunnerTests(unittest.TestCase):
    def test_runs_fake_amrfinder_per_sample_and_writes_status(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bin_dir = tmp / "bin"
            sequence_dir = tmp / "sequence"
            raw_dir = tmp / "raw"
            status_file = tmp / "tables" / "amrfinderplus_sample_status.tsv"
            bin_dir.mkdir()
            sequence_dir.mkdir()
            (sequence_dir / "GCF_000000001.1.fna").write_text(">c1\nATGC\n", encoding="utf-8")
            (sequence_dir / "GCF_000000002.1.fna").write_text(">c1\nATGC\n", encoding="utf-8")
            fake_amrfinder = bin_dir / "amrfinder"
            fake_amrfinder.write_text(
                "#!/usr/bin/env bash\n"
                "if [ \"$1\" = \"-u\" ]; then exit 0; fi\n"
                "out=\"\"\n"
                "while [ \"$#\" -gt 0 ]; do\n"
                "  if [ \"$1\" = \"-o\" ]; then out=\"$2\"; shift 2; continue; fi\n"
                "  shift\n"
                "done\n"
                "printf 'Gene symbol\\tClass\\nfoo\\tbeta-lactam\\n' > \"$out\"\n",
                encoding="utf-8",
            )
            fake_amrfinder.chmod(0o755)
            env = os.environ.copy()
            env["PATH"] = f"{bin_dir}:{env.get('PATH', '')}"

            subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "run_amrfinderplus_parallel.py"),
                    "--sequence-dir",
                    str(sequence_dir),
                    "--raw-dir",
                    str(raw_dir),
                    "--status-file",
                    str(status_file),
                    "--jobs",
                    "2",
                    "--threads-per-sample",
                    "1",
                    "--update-db",
                    "true",
                ],
                check=True,
                env=env,
            )

            status = pd.read_csv(status_file, sep="\t")
            self.assertEqual(set(status["status"]), {"PASS"})
            self.assertEqual(len(status), 2)
            self.assertTrue((raw_dir / "GCF_000000001.1.tsv").exists())
            self.assertTrue((raw_dir / "GCF_000000002.1.tsv").exists())
            self.assertTrue((raw_dir / "amrfinder_update.log").exists())


if __name__ == "__main__":
    unittest.main()
