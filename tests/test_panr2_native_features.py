import csv
import importlib.util
from pathlib import Path


def load_native_feature_module():
    module_path = Path(__file__).resolve().parents[1] / "scripts" / "run_panr2_native_features.py"
    spec = importlib.util.spec_from_file_location("run_panr2_native_features", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def load_panr2_contract_module():
    module_path = Path(__file__).resolve().parents[1] / "scripts" / "panr2_contract.py"
    spec = importlib.util.spec_from_file_location("panr2_contract", module_path)
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


def test_abricate_command_includes_external_datadir():
    native = load_native_feature_module()

    command = native.abricate_command("abricate", Path("/shared/abricate_db"), ["--db", "ncbi", "sample.fna"])

    assert command == ["abricate", "--datadir", "/shared/abricate_db", "--db", "ncbi", "sample.fna"]


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


def test_mobileelementfinder_failure_writes_header_status_audit_and_warning(tmp_path, capsys):
    native = load_native_feature_module()
    sample_dir = tmp_path / "sample"
    status_path = sample_dir / "panr2_native_feature_runners" / "module_status.tsv"
    audit_path = sample_dir / "panr2_native_feature_runners" / "native_runner_merge_audit.tsv"
    rows = []
    audit_rows = []

    status = native.write_mobileelementfinder_failure_outputs(
        sample_dir=sample_dir,
        status_path=status_path,
        audit_path=audit_path,
        rows=rows,
        audit_rows=audit_rows,
        runner_mode="parallel",
        started="2026-05-16T00:00:00Z",
        sample_count=2,
        allow_failure=True,
        error=RuntimeError("parser failed"),
    )

    assert status == "WARNING_FAILED"
    header_only = sample_dir / "tool_results" / "mobileelementfinder" / "panr2_inputs" / "mobileelementfinder_results.tab"
    assert header_only.read_text(encoding="utf-8").count("\n") == 1
    with status_path.open(newline="", encoding="utf-8") as handle:
        status_rows = list(csv.DictReader(handle, delimiter="\t"))
    with audit_path.open(newline="", encoding="utf-8") as handle:
        audit_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert status_rows[0]["status"] == "WARNING_FAILED"
    assert status_rows[0]["samples_failed"] == "2"
    assert audit_rows[0]["status"] == "WARNING_FAILED"
    assert "header-only PanR2-compatible output was written" in capsys.readouterr().err


def test_header_only_native_plasmidfinder_exports_feature_table(tmp_path):
    contract = load_panr2_contract_module()
    sample_dir = tmp_path / "Acinetobacter_pittii"
    plasmidfinder_dir = sample_dir / "tool_results" / "abricate" / "plasmidfinder"
    plasmidfinder_dir.mkdir(parents=True)
    (sample_dir / "metadata_output").mkdir(parents=True)
    (sample_dir / "metadata_output" / "ncbi_clean.csv").write_text(
        "Assembly Accession,Organism Name\nGCF_000000001.1,Acinetobacter pittii\n",
        encoding="utf-8",
    )
    (plasmidfinder_dir / "plasmidfinder_results.tab").write_text(
        "#FILE\tSEQUENCE\tSTART\tEND\tGENE\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n",
        encoding="utf-8",
    )

    written = contract.write_feature_tables(sample_dir, sample_dir / "panr2_inputs")

    plasmid_features = sample_dir / "panr2_inputs" / "features" / "plasmidfinder.features.tsv"
    assert written["plasmidfinder"] == str(plasmid_features)
    assert plasmid_features.exists()
    assert plasmid_features.read_text(encoding="utf-8").count("\n") == 1
