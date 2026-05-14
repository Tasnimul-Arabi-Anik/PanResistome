import importlib.util
import json
import sys
from argparse import Namespace
from pathlib import Path


def load_export_module():
    scripts_dir = Path(__file__).resolve().parents[1] / "scripts"
    if str(scripts_dir) not in sys.path:
        sys.path.insert(0, str(scripts_dir))
    module_path = scripts_dir / "export_panr2_inputs.py"
    spec = importlib.util.spec_from_file_location("export_panr2_inputs", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_reproducibility_manifest_records_core_run_context(tmp_path):
    exporter = load_export_module()
    sample_dir = tmp_path / "sample"
    out = sample_dir / "panr2_inputs"
    features = out / "features"
    manifest = out / "manifest"
    features.mkdir(parents=True)
    manifest.mkdir(parents=True)
    feature_text = (
        "sample_id\tassembly_accession\tdatabase\tfeature_id\tfeature_category\tpresence\tidentity\tcoverage\tcontig\tstart\tend\ttool\ttool_version\tdatabase_version\n"
        "s1\tGCF_000001.1\tamr\tblaX\tbeta-lactam\t1\t99\t100\tcontig1\t1\t100\tabricate\t1.0\tncbi\n"
    )
    (features / "amr.features.tsv").write_text(feature_text, encoding="utf-8")
    (features / "all_features.tsv").write_text(feature_text, encoding="utf-8")
    (manifest / "report_controls.tsv").write_text(
        "setting\tvalue\tmessage\nlarge_dataset\ttrue\tLarge-dataset safeguards enabled.\n",
        encoding="utf-8",
    )
    (manifest / "database_setup_status.tsv").write_text(
        "database_or_tool\tstatus\tmessage\nabricate_ncbi\tPASS\tok\n",
        encoding="utf-8",
    )
    (manifest / "module_status_summary.tsv").write_text(
        "module\tstatus\tmessage\namr\tPASS\tok\n",
        encoding="utf-8",
    )

    args = Namespace(
        repo_dir=str(tmp_path / "not-a-repo"),
        pipeline_outdir=str(tmp_path / "results"),
        pipeline_version="0.4.1-dev",
        run_command="nextflow run main.nf -profile docker,large",
        profile_stack="docker,large",
        launch_dir=str(tmp_path),
        container_engine="",
        container_image="ghcr.io/tasnimul-arabi-anik/panresistome@sha256:abc123",
        container_digest="",
        git_commit="",
        git_tag="",
        nextflow_session_id="session-1",
        nextflow_run_name="test-run",
    )

    path = exporter.write_reproducibility_manifest(sample_dir, out, args)

    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["pipeline"]["version"] == "0.4.1-dev"
    assert payload["pipeline"]["contract_version"] == "1.0"
    assert payload["execution"]["profile_stack"] == "docker,large"
    assert payload["container"]["engine"] == "docker"
    assert payload["container"]["digest"] == "sha256:abc123"
    assert payload["features"]["row_counts_by_table"]["amr"] == 1
    assert payload["features"]["row_counts_by_table"]["all_features"] == 1
    assert payload["report_controls"]["large_dataset"] == "true"
