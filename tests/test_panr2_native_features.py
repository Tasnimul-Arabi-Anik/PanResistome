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
