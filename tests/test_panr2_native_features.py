import importlib.util
from pathlib import Path


def load_native_feature_module():
    module_path = Path(__file__).resolve().parents[1] / "scripts" / "run_panr2_native_features.py"
    spec = importlib.util.spec_from_file_location("run_panr2_native_features", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_mlst_status_count_ignores_placeholder_rows(tmp_path):
    native = load_native_feature_module()
    raw = tmp_path / "mlst.tsv"
    raw.write_text(
        "/data/GCF_000001.1.fna\t-\t-\n"
        "/data/GCF_000002.1.fna\t-\t-\n",
        encoding="utf-8",
    )

    assert native.count_mlst_feature_rows(raw) == (0, 0)


def test_mlst_status_count_detects_st_and_alleles(tmp_path):
    native = load_native_feature_module()
    raw = tmp_path / "mlst.tsv"
    raw.write_text(
        "/data/GCF_000001.1.fna\tecoli\t42\tadk(1)\tfumC(2)\n"
        "/data/GCF_000002.1.fna\tecoli\t42\tadk(1)\tfumC(3)\n",
        encoding="utf-8",
    )

    assert native.count_mlst_feature_rows(raw) == (6, 4)


def test_worker_count_bounds_threads_to_tasks():
    native = load_native_feature_module()

    assert native.worker_count(16, 3) == 3
    assert native.worker_count(4, 20) == 4
    assert native.worker_count(0, 20) == 1
    assert native.worker_count(16, 0) == 1


def test_native_runner_audit_writer(tmp_path):
    native = load_native_feature_module()
    audit_path = tmp_path / "native_runner_merge_audit.tsv"
    row = native.audit_row(
        "abricate",
        "parallel",
        6,
        6,
        2,
        2,
        0,
        10,
        4,
        "PASS",
        "ABRicate completed",
    )

    native.write_audit(audit_path, [row])

    text = audit_path.read_text(encoding="utf-8")
    assert "expected_raw_tables\tobserved_raw_tables" in text
    assert "abricate\tparallel\t6\t6\t2\t2\t0\t10\t4\tPASS" in text
