#!/usr/bin/env python3
"""Repair reviewed missing TCCloud slots without rebuilding ordered input groups.

Defaults to a read-only dry run. --apply backs up SQLite, guards all original
rows, changes only file_path in existing fourth input slots, and validates
the whole catalog against its backup before committing one transaction.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from contextlib import closing
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import sqlite3


DEFAULT_DATABASE = Path("/data/keeling/a/gzhao1/f/Database/inputfiles_all.sqlite")
DEFAULT_AUDIT = Path(
    "/data/gdi/f/gzhao1/mmcth_validate/climate/"
    "pilot1_misr_availability_2010_2019/source_gap_audit"
)
TC_NAME = re.compile(r"^MISR_AM1_TC_CLOUD_P(\d{3})_O(\d{6})_F01_0001\.hdf$")
FILE_COLUMNS = (
    "file_id", "group_id", "file_path", "orbit_number", "date", "modis_id", "path_number"
)


def connect(database: Path, writable: bool = False) -> sqlite3.Connection:
    connection = sqlite3.connect(
        database.absolute().as_uri() + ("?mode=rw" if writable else "?mode=ro"),
        uri=True, timeout=60,
    )
    if not writable:
        connection.execute("PRAGMA query_only=ON")
    return connection


def fetch_module():
    path = Path(__file__).resolve().parents[1] / "misc" / "fetch_methods.py"
    spec = importlib.util.spec_from_file_location("repair_fetch_methods", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def inspect_groups(connection, candidates):
    """Check every group, including unavailable candidates, using stable row IDs."""
    plans = []
    unavailable = []
    already = []
    before_groups = {}
    paths_checked = set()
    group_ids = set()
    for candidate in candidates:
        group_id = int(candidate["catalog_group_id"])
        if group_id in group_ids:
            raise ValueError(f"Duplicate candidate group {group_id}")
        group_ids.add(group_id)
        rows = connection.execute(
            "SELECT file_id,group_id,file_path,orbit_number,date,modis_id,path_number "
            "FROM files WHERE group_id=? ORDER BY file_id", (group_id,)
        ).fetchall()
        if len(rows) != 9:
            raise ValueError(f"Expected nine inputs for group {group_id}; found {len(rows)}")
        for row in rows:
            if (row[3], row[4][:10], row[5], row[6]) != (
                int(candidate["orbit_number"]), candidate["date"],
                candidate["modis_id"], int(candidate["path_number"]),
            ):
                raise ValueError(f"Candidate metadata changed for group {group_id}")
        for index, column in ((0, "catalog_modis021km_path"),
                              (1, "catalog_modis06_path"), (2, "catalog_modis03_path")):
            if rows[index][2] != candidate[column]:
                raise ValueError(f"MODIS slot {index} changed for group {group_id}")
        if int(candidate["legacy_matching_file_count"]) != 1:
            if rows[3][2] != "X":
                raise ValueError(f"Unmatched TC slot changed for group {group_id}")
            unavailable.append(candidate)
            continue
        target_path = candidate["legacy_matching_file_paths"]
        match = TC_NAME.fullmatch(Path(target_path).name)
        if not match or (int(match[2]), int(match[1])) != (rows[3][3], rows[3][6]):
            raise ValueError(f"TCCloud filename/version/orbit/path mismatch: {target_path}")
        if target_path not in paths_checked:
            if not Path(target_path).is_file() or Path(target_path).stat().st_size <= 0:
                raise ValueError(f"Unavailable TCCloud file: {target_path}")
            paths_checked.add(target_path)
        if rows[3][2] == target_path:
            already.append(candidate)
            continue
        if rows[3][2] != "X":
            raise ValueError(f"Refusing to replace populated TC slot for group {group_id}")
        before_groups[group_id] = rows
        plans.append({
            "file_id": rows[3][0], "group_id": group_id, "date": candidate["date"],
            "modis_id": candidate["modis_id"], "orbit_number": rows[3][3],
            "path_number": rows[3][6], "old_file_path": "X", "new_file_path": target_path,
        })
    return plans, unavailable, already, before_groups


def table_counts(connection, schema="main"):
    return {name: connection.execute(f"SELECT count(*) FROM {schema}.{name}").fetchone()[0]
            for name in ("files", "groups")}


def validate_against_backup(connection, backup_path, plans):
    """Check all catalog rows and schema, allowing only the planned path edits."""
    connection.execute("ATTACH DATABASE ? AS before_repair",
                       (backup_path.absolute().as_uri() + "?mode=ro",))
    connection.execute(
        "CREATE TEMP TABLE repair_expected(file_id INTEGER PRIMARY KEY,new_path TEXT NOT NULL)"
    )
    connection.executemany("INSERT INTO repair_expected VALUES (?,?)",
                           [(row["file_id"], row["new_file_path"]) for row in plans])
    before = table_counts(connection, "before_repair")
    after = table_counts(connection)
    if before != after:
        raise ValueError(f"Catalog row counts changed: {before} -> {after}")
    changed = connection.execute("""
        SELECT count(*) FROM main.files f
        LEFT JOIN before_repair.files b ON f.file_id=b.file_id
        LEFT JOIN repair_expected e ON f.file_id=e.file_id
        WHERE b.file_id IS NULL OR f.group_id IS NOT b.group_id
           OR f.orbit_number IS NOT b.orbit_number OR f.date IS NOT b.date
           OR f.modis_id IS NOT b.modis_id OR f.path_number IS NOT b.path_number
           OR f.file_path IS NOT CASE WHEN e.file_id IS NULL THEN b.file_path ELSE e.new_path END
    """).fetchone()[0]
    if changed:
        raise ValueError(f"Unexpected files-table differences: {changed}")
    changed_groups = connection.execute("""
        SELECT count(*) FROM main.groups g
        LEFT JOIN before_repair.groups b ON g.group_id=b.group_id
        WHERE b.group_id IS NULL OR g.group_name IS NOT b.group_name
           OR g.group_datetime IS NOT b.group_datetime
    """).fetchone()[0]
    if changed_groups:
        raise ValueError(f"Unexpected groups-table differences: {changed_groups}")
    for query in (
        "SELECT type,name,tbl_name,sql FROM {schema}.sqlite_master ORDER BY type,name",
        "SELECT name,seq FROM {schema}.sqlite_sequence ORDER BY name",
    ):
        if connection.execute(query.format(schema="main")).fetchall() != connection.execute(
                query.format(schema="before_repair")).fetchall():
            raise ValueError("Catalog schema or SQLite sequences changed")
    integrity = connection.execute("PRAGMA main.integrity_check").fetchall()
    if integrity != [("ok",)]:
        raise ValueError(f"SQLite integrity check failed: {integrity[:10]}")
    return {"table_counts": after, "all_files_rows_match_expected": True,
            "all_groups_rows_unchanged": True, "schema_and_sequences_unchanged": True,
            "sqlite_integrity_check": "ok"}


def validate_lookup(connection, plans, before_groups):
    fetch = fetch_module()
    for row in plans:
        group_id = row["group_id"]
        expected = list(before_groups[group_id])
        replacement = list(expected[3])
        replacement[2] = row["new_file_path"]
        expected[3] = tuple(replacement)
        actual = connection.execute(
            "SELECT file_id,group_id,file_path,orbit_number,date,modis_id,path_number "
            "FROM files WHERE group_id=? ORDER BY file_id", (group_id,)
        ).fetchall()
        if actual != expected:
            raise ValueError(f"Group/order validation failed for {group_id}")
        paths, orbit = fetch.fetch_files_by_group(connection, group_id)
        if paths != [item[2] for item in expected] or orbit != row["orbit_number"]:
            raise ValueError(f"Production lookup validation failed for {group_id}")
    return {"groups_validated_with_production_lookup": len(plans),
            "nine_input_slots_and_file_ids_preserved": True}


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, default=DEFAULT_DATABASE)
    parser.add_argument("--candidates", type=Path,
                        default=DEFAULT_AUDIT / "legacy_tccloud_orbit_path_match_candidates.csv")
    parser.add_argument("--audit-dir", type=Path)
    parser.add_argument("--apply", action="store_true", help="Back up and apply one transaction.")
    args = parser.parse_args()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    audit_dir = args.audit_dir or DEFAULT_AUDIT / f"database_repair_{stamp}"
    audit_dir.mkdir(parents=True, exist_ok=False)
    candidates = list(csv.DictReader(args.candidates.open()))
    with closing(connect(args.database)) as readonly:
        plans, unavailable, already, before_groups = inspect_groups(readonly, candidates)
    write_csv(audit_dir / "planned_updates.csv", plans)
    write_csv(audit_dir / "unmatched_groups.csv", unavailable)
    by_date = Counter(row["date"] for row in plans)
    matched_dates = set(by_date) | {row["date"] for row in already}
    partial_unmatched = [row for row in unavailable if row["date"] in matched_dates]
    summary = {
        "created_utc": stamp, "database": str(args.database),
        "database_real_path": str(args.database.resolve()),
        "candidate_csv": str(args.candidates),
        "candidate_csv_sha256": hashlib.sha256(args.candidates.read_bytes()).hexdigest(),
        "mode": "apply" if args.apply else "dry_run", "committed": False,
        "planned_updates": len(plans), "already_repaired": len(already),
        "affected_dates": len(by_date), "unique_legacy_files": len({row["new_file_path"] for row in plans}),
        "updates_by_year": dict(Counter(row["date"][:4] for row in plans)),
        "updates_by_date": dict(sorted(by_date.items())),
        "unmatched_groups_on_repairable_dates": len(partial_unmatched),
        "unmatched_groups_by_repairable_date": dict(Counter(row["date"] for row in partial_unmatched)),
        "unmatched_groups_on_other_dates": len(unavailable) - len(partial_unmatched),
    }
    summary_path = audit_dir / "repair_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Validated {len(plans):,} updates on {len(by_date)} dates; audit: {audit_dir}", flush=True)
    if not args.apply or not plans:
        print(json.dumps(summary, indent=2), flush=True)
        return
    backup_path = args.database.with_name(args.database.name + f".backup_tccloud_{stamp}")
    if backup_path.exists():
        raise FileExistsError(backup_path)
    connection = connect(args.database, writable=True)
    try:
        connection.execute("PRAGMA synchronous=FULL")
        # A reserved writer lock excludes concurrent catalog mutations while
        # allowing the separate read connection to take a SQLite snapshot.
        connection.execute("BEGIN IMMEDIATE")
        live_plans, _, _, live_before = inspect_groups(connection, candidates)
        if live_plans != plans or live_before != before_groups:
            raise ValueError("Catalog changed between preflight and the writer lock")
        print(f"Creating consistent SQLite backup: {backup_path}", flush=True)
        with closing(connect(args.database)) as source, closing(sqlite3.connect(backup_path)) as target:
            source.backup(target, pages=4096)
        backup_path.chmod(0o444)
        summary["backup"] = str(backup_path)
        summary["backup_size_bytes"] = backup_path.stat().st_size
        summary_path.write_text(json.dumps(summary, indent=2) + "\n")
        print("Backup complete; updating the fourth input slot in place.", flush=True)
        for row in plans:
            cursor = connection.execute(
                "UPDATE files SET file_path=? WHERE file_id=? AND group_id=? AND file_path='X' "
                "AND orbit_number=? AND path_number=? AND date=? AND modis_id=?",
                (row["new_file_path"], row["file_id"], row["group_id"], row["orbit_number"],
                 row["path_number"], row["date"] + " 00:00:00", row["modis_id"]),
            )
            if cursor.rowcount != 1:
                raise ValueError(f"Guarded update failed for group {row['group_id']}")
        if connection.total_changes != len(plans):
            raise ValueError("Unexpected SQLite modification count")
        print("Checking every catalog row against the backup and SQLite integrity.", flush=True)
        validation = validate_against_backup(connection, backup_path, plans)
        validation.update(validate_lookup(connection, plans, before_groups))
        connection.commit()
        summary.update(committed=True, updated_rows=len(plans), validation=validation)
        # Reopen after commit and exercise the production lookup again.
        with closing(connect(args.database)) as final:
            summary["post_commit_lookup_validation"] = validate_lookup(final, plans, before_groups)
        write_csv(audit_dir / "updated_rows.csv", plans)
        for year in sorted({row["date"][:4] for row in plans}):
            ids = sorted({row["modis_id"] for row in plans if row["date"].startswith(year)})
            (audit_dir / f"repaired_modis_ids_{year}.txt").write_text("\n".join(ids) + "\n")
        summary_path.write_text(json.dumps(summary, indent=2) + "\n")
        print(json.dumps(summary, indent=2), flush=True)
    except BaseException as exc:
        if connection.in_transaction:
            connection.rollback()
        summary["error"] = f"{type(exc).__name__}: {exc}"
        summary_path.write_text(json.dumps(summary, indent=2) + "\n")
        raise
    finally:
        connection.close()


if __name__ == "__main__":
    main()
