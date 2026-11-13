import hashlib
import json
import subprocess
import sys
from pathlib import Path

from metadata.create_json_batch import initialize_worker, process_file


GENERATOR = '''
from pathlib import Path
import hashlib
def build_ummg(path, mod03_root, collection_shortname, collection_version, project=None, pge_version=None):
    source = Path(path)
    content = source.read_bytes()
    if content == b"broken":
        raise ValueError("unreadable granule")
    return {
        "GranuleUR": source.name,
        "CollectionReference": {"ShortName": collection_shortname, "Version": collection_version},
        "DataGranule": {"ArchiveAndDistributionInformation": [{
            "Name": source.name, "Size": len(content) / (1024.0 * 1024.0), "SizeUnit": "MB",
            "Checksum": {"Algorithm": "MD5", "Value": hashlib.md5(content).hexdigest()}
        }]}
    }, "matched_MOD03.hdf"
'''


def setup_worker(tmp_path):
    script = tmp_path / "generator.py"
    script.write_text(GENERATOR)
    initialize_worker(str(script), {
        "generator_mtime_ns": script.stat().st_mtime_ns, "overwrite": False,
        "mod03_root": str(tmp_path), "collection_shortname": "MODMISR_L2_CP",
        "collection_version": "V00", "project": None, "pge_version": None,
    })
    return script


def test_writes_adjacent_metadata_and_resumes_after_source_changes(tmp_path):
    setup_worker(tmp_path)
    source = tmp_path / "MODMISR_L2_CP_example_MM.nc"
    source.write_bytes(b"first granule")
    assert process_file(str(source))[0] == "written"
    sidecar = Path(str(source) + ".cmr.json")
    assert sidecar.exists()
    assert process_file(str(source))[0] == "skipped"
    source.write_bytes(b"changed granule")
    assert process_file(str(source))[0] == "written"
    data = json.loads(sidecar.read_text())
    assert data["DataGranule"]["ArchiveAndDistributionInformation"][0]["Checksum"]["Value"] == hashlib.md5(source.read_bytes()).hexdigest()


def test_repairs_invalid_json_and_preserves_previous_sidecar_on_failure(tmp_path):
    setup_worker(tmp_path)
    source = tmp_path / "MODMISR_L2_CP_example_MM.nc"
    source.write_bytes(b"valid")
    sidecar = Path(str(source) + ".cmr.json")
    sidecar.write_text('{"unfinished":')
    assert process_file(str(source))[0] == "written"
    previous = sidecar.read_bytes()
    source.write_bytes(b"broken")
    assert process_file(str(source))[0] == "failed"
    assert sidecar.read_bytes() == previous
    assert not list(tmp_path.glob(".*.tmp.*"))


def test_batch_finishes_other_files_and_exits_nonzero_when_one_fails(tmp_path):
    script = setup_worker(tmp_path)
    inputs = tmp_path / "2011"
    inputs.mkdir()
    (inputs / "MODMISR_L2_CP_good_MM.nc").write_bytes(b"valid")
    (inputs / "MODMISR_L2_CP_bad_MM.nc").write_bytes(b"broken")
    batch = Path(__file__).resolve().parents[1] / "metadata/create_json_batch.py"
    logs = tmp_path / "logs"
    result = subprocess.run([sys.executable, str(batch), "--input-dir", str(inputs), "--metadata-script", str(script),
                             "--workers", "2", "--log-dir", str(logs)], capture_output=True, text=True)
    assert result.returncode == 1, result.stdout + result.stderr
    summary = json.loads((logs / "summary.json").read_text())
    assert summary["files_selected"] == 2
    assert summary["written"] == summary["failed"] == 1
    assert (inputs / "MODMISR_L2_CP_good_MM.nc.cmr.json").exists()
    assert not (inputs / "MODMISR_L2_CP_bad_MM.nc.cmr.json").exists()
    assert "unreadable granule" in (logs / "failures.tsv").read_text()
