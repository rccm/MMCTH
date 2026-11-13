import json
import sqlite3
import subprocess
import sys
from pathlib import Path

from scripts.audit_archive_completion import matching_items, product_key


def test_filename_layouts_and_bracketed_orbits():
    expected = ("A2015214.1245", 83096, "MM")
    assert product_key("MODMISR_L2_CP_V00_061_A2015214.1245_0802_O083096_MM.nc") == expected
    assert product_key("MODMISR_L2_CP_061_A2015214.1245_0802_O083096_MM_V00.nc") == expected
    assert product_key("MODMISR_L2_CP_061_A2015214.1245_0802_O[83096]_MO_V00.nc") == (*expected[:2], "MO")
    assert product_key("MM_CO2FailureDiagnostics_A2015214.1245_O083096.nc") is None


def test_orbit_alias_is_only_allowed_without_a_tc_input():
    item = (("A2015214.1245", 83000, "MM"), "wrong.nc", 10, "")
    index = {"A2015214.1245": [item]}
    assert not matching_items(index, ("A2015214.1245", 83096, "MM"), True)
    mo = (("A2015214.1245", 83000, "MO"), "mo.nc", 10, "")
    index["A2015214.1245"].append(mo)
    assert matching_items(index, ("A2015214.1245", 83096, "MO"), False) == [mo]


def test_full_audit_does_not_hide_missing_groups_or_recount_resolved_failures(tmp_path):
    root = tmp_path / "products"
    year = root / "2005"
    year.mkdir(parents=True)
    legacy = tmp_path / "legacy"
    legacy.mkdir()
    database = tmp_path / "inputs.sqlite"
    report = tmp_path / "report"
    names = {}
    with sqlite3.connect(database) as conn:
        conn.execute("CREATE TABLE files (file_id INTEGER PRIMARY KEY, group_id INTEGER, file_path TEXT, orbit_number INTEGER, date TEXT, modis_id TEXT)")
        for gid in range(1, 9):
            stamp = f"A2005001.{gid:04d}"
            files = [f"/inputs/MOD021KM.{stamp}.hdf", f"/inputs/MOD06_L2.{stamp}.hdf", f"/inputs/MOD03.{stamp}.hdf",
                     f"/inputs/MISR_TC_CLOUD_P001_O{gid:06d}.hdf", "/inputs/AGP.hdf", "/inputs/single.nc",
                     "/inputs/profile.nc", "/inputs/single_next.nc", "/inputs/profile_next.nc"]
            names[gid] = year / f"MODMISR_L2_CP_V00_061_{stamp}_0101_O{gid:06d}_MM.nc"
            for file in files:
                conn.execute("INSERT INTO files (group_id,file_path,orbit_number,date,modis_id) VALUES (?,?,?,?,?)",
                             (gid, file, gid, "2005-01-01 00:00:00", f".{stamp}."))
    # Two present products, one expected skip, one failed, one unlogged missing,
    # one stored elsewhere, one partial fallback, and one invalid stored file.
    for gid in [1, 3]:
        names[gid].write_bytes(b"\x89HDF\r\n\x1a\n")
    (legacy / names[6].name).write_bytes(b"CDF\x01test")
    names[7].with_name(names[7].name.replace("_MM.nc", "_MO.nc")).write_bytes(b"CDF\x01test")
    names[8].write_text("an error page, not a product")
    old = year / "_run_logs" / "z_old"
    new = year / "_run_logs" / "a_new"
    old.mkdir(parents=True)
    new.mkdir()
    def event(gid, status, timestamp):
        return {"gid": gid, "status": status, "expected_output": str(names[gid]), "timestamp": timestamp}
    (old / "events_rank_0000.jsonl").write_text("\n".join(json.dumps(event(g, s, "2026-01-01T00:00:00+00:00")) for g, s in
        [(2, "FAILED_NO_OUTPUT"), (3, "FAILED_NO_OUTPUT"), (4, "FAILED_NO_OUTPUT")]) + "\n")
    (new / "events_rank_0000.jsonl").write_text(json.dumps(event(2, "SKIPPED_NO_OVERLAP", "2026-02-01T00:00:00+00:00")) + "\n")
    before = database.read_bytes()
    script = Path(__file__).resolve().parents[1] / "scripts" / "audit_archive_completion.py"
    result = subprocess.run([sys.executable, str(script), "--database", str(database), "--output-root", str(root),
                             "--supplemental-root", str(legacy), "--start-date", "2005-01-01", "--end-date", "2005-12-31",
                             "--report-dir", str(report), "--workers", "1"], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    summary = json.loads((report / "summary.json").read_text())
    assert summary["totals"]["database_groups"] == 8
    assert summary["totals"]["PRESENT"] == 2
    assert summary["expected_no_overlap_groups"] == 1
    assert summary["totals"]["CONFIRMED_FAILED"] == 1
    assert summary["totals"]["MISSING_WITHOUT_FAILURE_RECORD"] == 1
    assert summary["totals"]["PARTIAL_PRODUCT_ONLY"] == 1
    assert summary["totals"]["INVALID_STORED_OUTPUT"] == 1
    assert summary["stored_outside_production_groups"] == 1
    assert summary["unresolved_groups"] == 4
    ids = (report / "rerun_candidate_ids.txt").read_text().splitlines()
    assert ids == [f".A2005001.{g:04d}." for g in [4, 5, 7, 8]]
    assert database.read_bytes() == before
