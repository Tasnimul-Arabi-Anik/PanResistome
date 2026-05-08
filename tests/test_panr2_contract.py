import importlib.util
from pathlib import Path


def load_contract_module():
    module_path = Path(__file__).resolve().parents[1] / "scripts" / "panr2_contract.py"
    spec = importlib.util.spec_from_file_location("panr2_contract", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_headerless_mlst_output_exports_features(tmp_path):
    contract = load_contract_module()
    sample_dir = tmp_path / "sample"
    raw_dir = sample_dir / "tool_results" / "mlst" / "raw"
    raw_dir.mkdir(parents=True)
    (raw_dir / "mlst.tsv").write_text(
        "/data/GCF_000001.1_sample.fna\tecoli\t42\tadk(1)\tfumC(2)\n",
        encoding="utf-8",
    )

    written = contract.write_feature_tables(sample_dir, sample_dir / "panr2_inputs")

    assert "mlst" in written
    rows = contract.read_table(Path(written["mlst"]))
    assert {row["feature_id"] for row in rows} == {"ST_42", "adk_1", "fumC_2"}
    assert all(row["database"] == "mlst" for row in rows)


def test_headerless_mlst_placeholders_write_empty_feature_table(tmp_path):
    contract = load_contract_module()
    sample_dir = tmp_path / "sample"
    raw_dir = sample_dir / "tool_results" / "mlst" / "raw"
    raw_dir.mkdir(parents=True)
    (raw_dir / "mlst.tsv").write_text(
        "/data/GCF_000001.1_sample.fna\t-\t-\n"
        "/data/GCF_000002.1_sample.fna\t-\t-\n",
        encoding="utf-8",
    )

    written = contract.write_feature_tables(sample_dir, sample_dir / "panr2_inputs")

    assert "mlst" in written
    rows = contract.read_table(Path(written["mlst"]))
    assert rows == []
    header = contract.read_header(Path(written["mlst"]))
    assert "sample_id" in header
    assert "feature_id" in header
