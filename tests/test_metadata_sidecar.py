import hashlib
import json
import logging
import os
import stat
from pathlib import Path
from types import SimpleNamespace

import pytest

from src.data_processors import metadata_sidecar


GENERATOR = '''
import hashlib
import os
from pathlib import Path
with Path(__file__).with_suffix(".imports").open("a") as stream:
    stream.write(str(os.getpid()) + "\\n")
def build_ummg(path, mod03_root, collection_shortname, collection_version):
    source = Path(path)
    content = source.read_bytes()
    return {
        "GranuleUR": source.name,
        "CollectionReference": {"ShortName": collection_shortname, "Version": collection_version},
        "DataGranule": {"ArchiveAndDistributionInformation": [{
            "Name": source.name, "Size": len(content) / (1024.0 * 1024.0), "SizeUnit": "MB",
            "Checksum": {"Algorithm": "MD5", "Value": hashlib.md5(content).hexdigest()}
        }]}
    }, "matched_MOD03.hdf"
'''


@pytest.fixture
def enabled(tmp_path, monkeypatch):
    script = tmp_path / "generator.py"
    script.write_text(GENERATOR)
    monkeypatch.setenv("CREATE_JSON_METADATA", "1")
    monkeypatch.setenv("METADATA_SCRIPT", str(script))
    monkeypatch.setenv("MOD03_ROOT", str(tmp_path / "mod03"))
    monkeypatch.delenv("COLLECTION_SHORTNAME", raising=False)
    monkeypatch.delenv("COLLECTION_VERSION", raising=False)
    metadata_sidecar._load_generator.cache_clear()
    yield script
    metadata_sidecar._load_generator.cache_clear()


def test_disabled_without_environment_does_not_open_source_or_generator(tmp_path, monkeypatch):
    monkeypatch.delenv("CREATE_JSON_METADATA", raising=False)
    monkeypatch.setenv("METADATA_SCRIPT", str(tmp_path / "missing_generator.py"))
    assert metadata_sidecar.create_json_sidecar(tmp_path / "missing.nc") is None
    assert list(tmp_path.iterdir()) == []


def test_fresh_sidecars_have_correct_checksums_and_import_generator_once(tmp_path, enabled):
    for i in range(2):
        source = tmp_path / f"granule_{i}.nc"
        source.write_bytes(f"complete NetCDF bytes {i}".encode())
        source.chmod(0o640)
        signature = source.stat()
        sidecar = metadata_sidecar.create_json_sidecar(source)
        assert sidecar == Path(str(source) + ".cmr.json")
        metadata = json.loads(sidecar.read_text())
        archive = metadata["DataGranule"]["ArchiveAndDistributionInformation"][0]
        assert archive["Checksum"]["Value"] == hashlib.md5(source.read_bytes()).hexdigest()
        assert metadata["CollectionReference"] == {"ShortName": "MODMISR_L2_CP", "Version": "V00"}
        assert stat.S_IMODE(sidecar.stat().st_mode) == 0o640
        assert source.stat().st_mtime_ns == signature.st_mtime_ns
    assert enabled.with_suffix(".imports").read_text().splitlines() == [str(os.getpid())]


@pytest.mark.parametrize("failure", ["generation", "serialization", "replacement"])
def test_metadata_failure_keeps_netcdf_and_previous_sidecar(tmp_path, enabled, monkeypatch, caplog, failure):
    source = tmp_path / "good.nc"
    content = b"successful scientific output"
    source.write_bytes(content)
    sidecar = Path(str(source) + ".cmr.json")
    previous = '{"previous": "metadata"}\n'
    sidecar.write_text(previous)
    if failure == "generation":
        enabled.write_text(GENERATOR + '\ndef build_ummg(*args):\n    raise OSError("metadata input unavailable")\n')
    elif failure == "serialization":
        enabled.write_text(GENERATOR + '\ndef build_ummg(*args):\n    return {"bad": object()}, "MOD03.hdf"\n')
    else:
        def fail_replace(*args):
            raise OSError("metadata publication failed")
        monkeypatch.setattr(metadata_sidecar.os, "replace", fail_replace)
    assert metadata_sidecar.create_json_sidecar(source) is None
    assert source.read_bytes() == content
    assert sidecar.read_text() == previous
    assert not list(tmp_path.glob(".*.tmp.*"))
    assert "NetCDF remains saved" in caplog.text


def test_does_not_publish_metadata_if_source_changes(tmp_path, enabled, caplog):
    enabled.write_text(GENERATOR + '''
def build_ummg(path, *args):
    with open(path, "ab") as stream:
        stream.write(b"changed by another writer")
    return {"GranuleUR": Path(path).name}, "MOD03.hdf"
''')
    source = tmp_path / "changing.nc"
    source.write_bytes(b"original")
    assert metadata_sidecar.create_json_sidecar(source) is None
    assert not Path(str(source) + ".cmr.json").exists()
    assert "changed while generating JSON" in caplog.text


def test_replaces_stale_metadata_when_netcdf_is_rewritten(tmp_path, enabled):
    source = tmp_path / "rewritten.nc"
    source.write_bytes(b"old data")
    sidecar = metadata_sidecar.create_json_sidecar(source)
    old_metadata = sidecar.read_text()
    source.write_bytes(b"new scientific output")
    assert metadata_sidecar.create_json_sidecar(source) == sidecar
    metadata = json.loads(sidecar.read_text())
    assert sidecar.read_text() != old_metadata
    assert metadata["DataGranule"]["ArchiveAndDistributionInformation"][0]["Checksum"]["Value"] == hashlib.md5(source.read_bytes()).hexdigest()


@pytest.mark.parametrize("has_mod06,has_tc,flag", [(True, True, "MM"), (True, False, "MO"), (False, True, "MI")])
def test_production_save_creates_real_metadata_after_netcdf_closes(tmp_path, enabled, monkeypatch, caplog, has_mod06, has_tc, flag):
    import numpy as np
    from netCDF4 import Dataset
    from src.data_processors import process_main

    script = Path(__file__).resolve().parents[1] / "metadata/creat_json.py"
    monkeypatch.setenv("METADATA_SCRIPT", str(script))
    processor = process_main.MainProcessor.__new__(process_main.MainProcessor)
    processor.log = logging.getLogger("test.inline_metadata")
    processor.input_files = ["MOD03.A2015161.0555.hdf"]
    processor.mm_processor = SimpleNamespace(orbit="O082320")
    processor.has_mod06 = has_mod06
    processor.has_tc = has_tc
    processor.met_date = np.array([2015, 6, 10, 0])
    processor.timestamp = "0555"
    shape = (12, 14)
    geo = {
        "lat": np.linspace(30, 31, 168, dtype=np.float32).reshape(shape),
        "lon": np.linspace(-100, -98, 168, dtype=np.float32).reshape(shape),
        "vza": np.zeros(shape, dtype=np.float32),
    }
    source = tmp_path / f"MODMISR_L2_CP_V00_061_A2015161.0555_0610_O082320_{flag}.nc"
    checked = []
    real_helper = metadata_sidecar.create_json_sidecar
    def check_closed_then_create(path, logger=None):
        before = Path(path).read_bytes()
        # Reopening the actual saved file here verifies synchronous completion.
        with Dataset(path) as dataset:
            assert "MM_CloudOpticalDepth" in dataset.variables
        result = real_helper(path, logger=logger)
        assert result is not None
        assert Path(path).read_bytes() == before
        checked.append(hashlib.md5(before).hexdigest())
        return result
    monkeypatch.setattr(process_main, "create_json_sidecar", check_closed_then_create)
    processor.save_pixels(
        {"opt": np.ones(shape, dtype=np.float32), "cth": np.full(shape, 1000, dtype=np.int16)},
        None, geo, np.full(shape, 1000, dtype=np.int16),
        {"surface_pressure": np.full(shape, 100000.0, dtype=np.float32)}, str(source),
    )
    metadata = json.loads(Path(str(source) + ".cmr.json").read_text())
    assert checked == [metadata["DataGranule"]["ArchiveAndDistributionInformation"][0]["Checksum"]["Value"]]
    assert metadata["GranuleUR"] == source.name
    bounds = metadata["SpatialExtent"]["HorizontalSpatialDomain"]["Geometry"]["BoundingRectangles"][0]
    assert bounds == {"WestBoundingCoordinate": -100.0, "SouthBoundingCoordinate": 30.0,
                      "EastBoundingCoordinate": -98.0, "NorthBoundingCoordinate": 31.0}
    with Dataset(source) as dataset:
        assert metadata["DataGranule"]["ProductionDateTime"] == dataset.getncattr("date_created")
    assert "fallback sensing time/day-night" in caplog.text
