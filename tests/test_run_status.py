import json
from pathlib import Path

from src.misc.run_status import RunStatusWriter, expected_output_path, work_item_record
from src.data_readers import mod_read
from scripts.build_mm_rerun_list import audit_run


def _files(tmp_path):
    names = [
        "MOD021KM.A2015214.1245.061.example.hdf",
        "MOD06_L2.A2015214.1245.061.example.hdf",
        "MOD03.A2015214.1245.061.example.hdf",
        "MISR_AM1_TC_CLOUD_P062_O083096_F01_0001.hdf",
        "MISR_AM1_AGP_P062_F01_24.hdf",
        "era5_single_2015_08_02.nc",
        "era5_profile_2015_08_02.nc",
        "era5_single_2015_08_03.nc",
        "era5_profile_2015_08_03.nc",
    ]
    paths = []
    for name in names:
        path = tmp_path / name
        path.touch()
        paths.append(str(path))
    return paths


def test_expected_output_path_matches_production_name(tmp_path):
    files = _files(tmp_path)

    output = expected_output_path(files, 83096, str(tmp_path / "output"))

    assert Path(output).name == (
        "MODMISR_L2_CP_V00_061_A2015214.1245_0802_O083096_MM.nc"
    )


def test_run_status_writer_records_plan_and_events_per_rank(tmp_path):
    files = _files(tmp_path)
    output_dir = tmp_path / "output"
    item = (123, files, 83096)
    writer = RunStatusWriter(str(output_dir), "test_run", rank=7)

    writer.write_plan([item], str(output_dir))
    base = work_item_record(*item, str(output_dir))
    writer.event(base, "STARTED")
    writer.event(base, "FAILED_NO_OUTPUT", message="MOD03_READ_FAILED")
    writer.close()

    plan = [json.loads(line) for line in writer.plan_path.read_text().splitlines()]
    events = [json.loads(line) for line in writer.events_path.read_text().splitlines()]
    assert plan[0]["gid"] == 123
    assert plan[0]["modis_id"] == ".A2015214.1245."
    assert [event["status"] for event in events] == [
        "STARTED",
        "FAILED_NO_OUTPUT",
    ]


def test_hdf_open_retries_transient_failures(monkeypatch):
    attempts = []
    sentinel = object()

    def fake_sd(path, mode):
        attempts.append((path, mode))
        if len(attempts) < 3:
            raise OSError("temporary filesystem error")
        return sentinel

    monkeypatch.setattr(mod_read, "SD", fake_sd)
    monkeypatch.setattr(mod_read.time, "sleep", lambda _seconds: None)

    result = mod_read._open_sd_with_retry("example.hdf", attempts=3)

    assert result is sentinel
    assert len(attempts) == 3


def test_run_audit_selects_failed_and_not_started_but_not_no_overlap(tmp_path):
    files = _files(tmp_path)
    output_dir = tmp_path / "output"
    writer = RunStatusWriter(str(output_dir), "audit_run", rank=0)
    items = [(1, files, 83096), (2, files, 83096), (3, files, 83096)]
    writer.write_plan(items, str(output_dir))

    failed = work_item_record(*items[0], str(output_dir), check_input_exists=False)
    skipped = work_item_record(*items[1], str(output_dir), check_input_exists=False)
    writer.event(failed, "STARTED")
    writer.event(failed, "FAILED_NO_OUTPUT")
    writer.event(skipped, "STARTED")
    writer.event(skipped, "SKIPPED_NO_OVERLAP")
    writer.close()

    details = audit_run(writer.run_dir, required_decoder="")

    assert {detail["gid"] for detail in details} == {1, 3}
    assert {detail["reason"] for detail in details} == {
        "FAILED_NO_OUTPUT",
        "NOT_STARTED",
    }

