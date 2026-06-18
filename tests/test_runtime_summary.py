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


class CleanupAfterRunTests(unittest.TestCase):
    def test_removes_work_dir_session_cache_and_transient_caches(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            launch = tmp / "launch"
            outdir = tmp / "results"
            work_dir = tmp / "work"
            summary = outdir / "pipeline_cleanup_summary.tsv"
            session_cache = launch / ".nextflow" / "cache" / "session-123"
            scripts_cache = launch / "scripts" / "__pycache__"
            tests_cache = launch / "tests" / "__pycache__"
            pytest_cache = launch / ".pytest_cache"

            for directory in [outdir, work_dir, session_cache, scripts_cache, tests_cache, pytest_cache]:
                directory.mkdir(parents=True)
                (directory / "marker.txt").write_text("temporary\n", encoding="utf-8")

            subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "cleanup_after_run.py"),
                    "--work-dir",
                    str(work_dir),
                    "--launch-dir",
                    str(launch),
                    "--outdir",
                    str(outdir),
                    "--summary",
                    str(summary),
                    "--nextflow-session-id",
                    "session-123",
                ],
                check=True,
            )

            self.assertFalse(work_dir.exists())
            self.assertFalse(session_cache.exists())
            self.assertFalse(scripts_cache.exists())
            self.assertFalse(tests_cache.exists())
            self.assertFalse(pytest_cache.exists())
            self.assertTrue(outdir.exists())

            rows = pd.read_csv(summary, sep="\t", dtype=str)
            self.assertIn("nextflow_work_dir", set(rows["target"]))
            self.assertIn("nextflow_session_cache", set(rows["target"]))
            self.assertEqual(set(rows["status"]), {"PASS"})
            self.assertGreater(int(rows.loc[rows["target"] == "nextflow_work_dir", "bytes_removed"].iloc[0]), 0)

    def test_refuses_to_remove_launch_or_output_directory(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            launch = tmp / "launch"
            outdir = tmp / "results"
            summary = outdir / "pipeline_cleanup_summary.tsv"
            launch.mkdir()
            outdir.mkdir()
            (launch / "keep.txt").write_text("do not delete\n", encoding="utf-8")
            (outdir / "keep.txt").write_text("published output\n", encoding="utf-8")

            subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "cleanup_after_run.py"),
                    "--work-dir",
                    str(launch),
                    "--launch-dir",
                    str(launch),
                    "--outdir",
                    str(outdir),
                    "--summary",
                    str(summary),
                ],
                check=True,
            )

            self.assertTrue((launch / "keep.txt").exists())
            self.assertTrue((outdir / "keep.txt").exists())
            rows = pd.read_csv(summary, sep="\t", dtype=str)
            work_row = rows.loc[rows["target"] == "nextflow_work_dir"].iloc[0]
            self.assertEqual(work_row["status"], "SKIPPED_UNSAFE")
            self.assertIn("protected path", work_row["message"])


class AniSummaryTests(unittest.TestCase):
    def test_empty_pairs_with_genome_list_keeps_samples(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            sample_dir = tmp / "sample"
            ani_dir = sample_dir / "ani"
            ani_dir.mkdir(parents=True)
            pairs = ani_dir / "fastani_pairs.tsv"
            genomes = ani_dir / "genomes.list"
            pairs.write_text("query\treference\tani\tfragments_mapped\tfragments_total\n", encoding="utf-8")
            genomes.write_text(
                f"{tmp}/GCF_000000001.1_ASM1_genomic.fna\n{tmp}/GCF_000000002.1_ASM2_genomic.fna\n",
                encoding="utf-8",
            )

            subprocess.run(
                [
                    sys.executable,
                    str(REPO_ROOT / "scripts" / "ani_summary.py"),
                    "--sample-dir",
                    str(sample_dir),
                    "--pairs",
                    str(pairs),
                    "--genomes-list",
                    str(genomes),
                    "--tool",
                    "fastani",
                ],
                check=True,
            )

            panr2 = pd.read_csv(sample_dir / "ani" / "analysis" / "panr2_ani_summary.csv")
            closest = pd.read_csv(sample_dir / "ani" / "analysis" / "closest_genome.csv")
            clusters = pd.read_csv(sample_dir / "ani" / "analysis" / "duplicate_clusters.csv")
            self.assertEqual(len(panr2), 2)
            self.assertEqual(set(panr2["feature_id"]), {"ANI_CLUSTER_0001", "ANI_CLUSTER_0002"})
            self.assertEqual(set(closest["species_consistency_status"]), {"WARN"})
            self.assertEqual(set(clusters["cluster_size"]), {1})


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

    def test_reuses_existing_sample_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bin_dir = tmp / "bin"
            sequence_dir = tmp / "sequence"
            raw_dir = tmp / "raw"
            status_file = tmp / "tables" / "amrfinderplus_sample_status.tsv"
            bin_dir.mkdir()
            sequence_dir.mkdir()
            raw_dir.mkdir()
            (sequence_dir / "GCF_000000001.1.fna").write_text(">c1\nATGC\n", encoding="utf-8")
            existing = raw_dir / "GCF_000000001.1.tsv"
            existing.write_text("Gene symbol\tClass\nexisting\tbeta-lactam\n", encoding="utf-8")
            fake_amrfinder = bin_dir / "amrfinder"
            fake_amrfinder.write_text(
                "#!/usr/bin/env bash\n"
                "if [ \"$1\" = \"-u\" ]; then exit 0; fi\n"
                "exit 99\n",
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
                    "1",
                    "--threads-per-sample",
                    "1",
                    "--update-db",
                    "true",
                ],
                check=True,
                env=env,
            )

            status = pd.read_csv(status_file, sep="\t")
            self.assertEqual(status.loc[0, "status"], "PASS")
            self.assertEqual(status.loc[0, "message"], "existing_output_reused")
            self.assertEqual(existing.read_text(encoding="utf-8").splitlines()[1], "existing\tbeta-lactam")


if __name__ == "__main__":
    unittest.main()
