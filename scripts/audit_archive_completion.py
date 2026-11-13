#!/usr/bin/env python3
"""Compare a read-only MMCTH input catalog with actual production products.

Check every product's size and sample NetCDF/HDF5 signatures. Reconcile absent products
with the newest matching per-group event, rather than requiring saved job logs.
Accept both historical filename layouts; keep alternate/test outputs separate.
This is a completion audit, not an exhaustive NetCDF or scientific-data audit.
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import os
import re
import sqlite3
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo

STAMP = re.compile(r"A\d{7}\.\d{4}")
PRODUCT = re.compile(r"(A\d{7}\.\d{4})_\d{4}_O\[?(\d+)\]?_(MM|MI|MO|NW)(?:_V\d+)?\.nc$", re.I)
ORBIT = re.compile(r"O(\d{6})")
FIELDS = ["year", "date", "gid", "modis_id", "product_type", "status", "last_status",
          "last_run", "last_event_at", "message", "recorded_status", "recorded_product_type",
          "recorded_expected_output", "expected_output", "stored_matches",
          "other_product_types", "other_stored_paths", "missing_required_inputs", "input_check",
          "modis_file", "mod06_file", "mod03_file", "misr_file", "agp_file",
          "era5_single", "era5_profile", "era5_single_next", "era5_profile_next"]


def product_key(name):
    match = PRODUCT.search(name)
    if match is None:
        return None
    stamp, orbit, flag = match.groups()
    return stamp, int(orbit), flag.upper()


def scan_directory(spec):
    """One worker per production year (or supplemental root). No shared writes."""
    directory, primary, first_year, last_year = spec
    items = []
    signature_checks = 0
    unrecognized = []
    if not os.path.isdir(directory):
        return directory, primary, items, unrecognized
    # Deliberately exclude debug, plots, temp, and run manifests.
    with os.scandir(directory) as entries:
        for entry in entries:
            if not entry.name.endswith(".nc") or not entry.name.startswith("MODMISR_"):
                continue
            key = product_key(entry.name)
            if key is None:
                unrecognized.append(entry.path)
                continue
            if not first_year <= int(key[0][1:5]) <= last_year:
                continue
            problem = ""
            checked = False
            try:
                stat = entry.stat()
                size = stat.st_size
                if size == 0:
                    problem = "EMPTY_OUTPUT"
                elif size < 65536 or signature_checks < 8:
                    # Shared storage advertises a 1 MiB block size; Python's
                    # default buffering would read that much for eight bytes.
                    with open(entry.path, "rb", buffering=0) as stream:
                        signature = stream.read(8)
                    checked = True
                    signature_checks += 1
                    if signature != b"\x89HDF\r\n\x1a\n" and signature[:4] not in (b"CDF\x01", b"CDF\x02", b"CDF\x05"):
                        problem = "INVALID_OUTPUT_SIGNATURE"
            except OSError as exc:
                size = -1
                problem = f"UNREADABLE_OUTPUT:{type(exc).__name__}:{exc}"
            items.append((key, entry.path, size, problem, checked))
    return directory, primary, items, unrecognized


def expected_record(rows, output_root):
    gid, _path, orbit, date, catalog_id = rows[0]
    files = [r[1] for r in rows]
    if len(files) != 9 or "MOD021KM" not in os.path.basename(files[0]):
        raise ValueError(f"Invalid database group {gid}: {len(files)} rows")
    match = STAMP.search(files[0])
    if not match:
        raise ValueError(f"Missing MODIS timestamp in group {gid}")
    stamp = match.group()
    has_mod06 = os.path.basename(files[1]) != "X"
    has_tc = os.path.basename(files[3]) != "X"
    flag = "MM" if has_mod06 and has_tc else "MI" if not has_mod06 else "MO"
    tc_orbit = ORBIT.search(files[3])
    orbit = int(tc_orbit.group(1)) if tc_orbit else int(orbit) if orbit is not None else None
    year = int(stamp[1:5])
    key = (stamp, orbit, flag)
    expected = f"{output_root}/{year}/MODMISR_L2_CP_V00_061_{stamp}_{date[5:7]}{date[8:10]}_O{orbit:06d}_{flag}.nc" if orbit is not None else ""
    record = {"year": year, "date": date[:10], "gid": gid, "modis_id": f".{stamp}.",
              "product_type": flag, "expected_output": expected}
    for name, path in zip(FIELDS[-9:], files):
        record[name] = path
    return key, record, has_tc, has_mod06


def matching_items(index, key, has_tc):
    # Date-range runs historically used the first orbit in the chunk for MO
    # names. Without a TC file, match MODIS timestamp + product type instead.
    if not has_tc:
        return [item for item in index.get(key[0], []) if item[0][2] == key[2]]
    return [item for item in index.get(key[0], []) if item[0] == key]


def scan_run_records(spec):
    year, root, wanted = spec
    latest = {}
    malformed = []
    runs = []
    run_root = Path(root) / str(year) / "_run_logs"
    if not run_root.exists():
        return year, latest, runs, malformed
    for run in sorted(run_root.iterdir()):
        if not run.is_dir():
            continue
        runs.append(run.name)
        for path in sorted(run.glob("events_rank_*.jsonl")):
            with path.open() as stream:
                for number, line in enumerate(stream, 1):
                    if not line.strip():
                        continue
                    try:
                        row = json.loads(line)
                    except json.JSONDecodeError:
                        malformed.append(f"{path}:{number}")
                        continue
                    gid = int(row["gid"])
                    if gid not in wanted:
                        continue
                    order = (row.get("timestamp", ""), row.get("status") != "STARTED")
                    old = latest.get(gid)
                    if old is None or order >= old[0]:
                        latest[gid] = (order, {**row, "run": run.name})
    # Only consult plans if no event survived for an affected group.
    without_event = wanted - latest.keys()
    if without_event:
        for run in sorted(run_root.iterdir()):
            for path in sorted(run.glob("planned_rank_*.jsonl")):
                with path.open() as stream:
                    for number, line in enumerate(stream, 1):
                        if not line.strip():
                            continue
                        try:
                            row = json.loads(line)
                        except json.JSONDecodeError:
                            malformed.append(f"{path}:{number}")
                            continue
                        gid = int(row["gid"])
                        if gid not in without_event:
                            continue
                        order = (row.get("planned_at", ""), False)
                        if gid not in latest or order >= latest[gid][0]:
                            latest[gid] = (order, {**row, "status": "NOT_STARTED", "run": run.name})
    return year, {gid: item[1] for gid, item in latest.items()}, runs, malformed


def write_tsv(path, rows, fields=FIELDS):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def check_inputs(row):
    roles = list(FIELDS[-9:])
    required = {0, 2, 5, 6, 7, 8}
    if row["product_type"] in {"MM", "MO"}:
        required.add(1)
    if row["product_type"] in {"MM", "MI"}:
        required.update({3, 4})
    row["missing_required_inputs"] = ";".join(row[roles[i]] for i in sorted(required) if not os.path.isfile(row[roles[i]]))
    row["input_check"] = "REQUIRED_INPUTS_EXIST" if not row["missing_required_inputs"] else "REQUIRED_INPUTS_MISSING"
    for product, role in [("MOD06", "mod06_file"), ("MOD03", "mod03_file")]:
        if product + "_READ_FAILED" in row.get("message", "") and os.path.isfile(row[role]):
            try:
                from pyhdf.SD import SD, SDC
                dataset = SD(row[role], SDC.READ)
                dataset.end()
                row["input_check"] += f";{product}_HDF_OPEN_OK_NOW"
            except Exception as exc:
                row["input_check"] += f";{product}_HDF_OPEN_FAILED:{type(exc).__name__}:{exc}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--database", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--supplemental-root", type=Path, action="append", default=[])
    parser.add_argument("--start-date", default="2000-03-01")
    parser.add_argument("--end-date", default="2022-03-31")
    parser.add_argument("--report-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    datetime.strptime(args.start_date, "%Y-%m-%d")
    datetime.strptime(args.end_date, "%Y-%m-%d")
    root = str(args.output_root.resolve())
    dest = args.report_dir.resolve()
    dest.mkdir(parents=True, exist_ok=True)
    started = datetime.now(ZoneInfo("America/Chicago")).isoformat(timespec="seconds")
    db_before = args.database.stat()
    first_year, last_year = int(args.start_date[:4]), int(args.end_date[:4])
    primary, supplemental = defaultdict(list), defaultdict(list)
    inventory_stats = []
    bad_outputs = []
    unknown_names = []
    specs = [(f"{root}/{y}", True, first_year, last_year) for y in range(first_year, last_year + 1)]
    specs.append((root, True, first_year, last_year))
    specs.extend((str(p.resolve()), False, first_year, last_year) for p in args.supplemental_root)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(scan_directory, spec) for spec in specs]
        for future in as_completed(futures):
            directory, is_primary, items, unknown = future.result()
            index = primary if is_primary else supplemental
            for item in items:
                index[item[0][0]].append(item)
                if item[3]:
                    bad_outputs.append({"path": item[1], "bytes": item[2], "problem": item[3], "primary": is_primary})
            unknown_names.extend(unknown)
            inventory_stats.append({"directory": directory, "primary": is_primary,
                                    "products": len(items), "bad_signatures_or_empty": sum(bool(i[3]) for i in items),
                                    "signature_checks": sum(bool(i[4]) for i in items),
                                    "unrecognized_names": len(unknown)})
            print("INVENTORY", directory, len(items), "products;", sum(bool(i[3]) for i in items), "invalid", flush=True)
    write_tsv(dest / "invalid_stored_outputs.tsv", bad_outputs, ["path", "bytes", "problem", "primary"])
    (dest / "unrecognized_output_names.txt").write_text("".join(p + "\n" for p in unknown_names))
    missing = []
    yearly = defaultdict(Counter)
    database_errors = []
    group_count = 0
    print("Comparing database groups with stored products", flush=True)
    uri = args.database.resolve().as_uri() + "?mode=ro"
    with sqlite3.connect(uri, uri=True) as conn:
        query = "SELECT group_id,file_path,orbit_number,date,modis_id FROM files WHERE date>=? AND date<=? ORDER BY group_id,file_id"
        cursor = conn.execute(query, (args.start_date + " 00:00:00", args.end_date + " 23:59:59"))
        for gid, grouped in itertools.groupby(cursor, key=lambda r: r[0]):
            rows = list(grouped)
            group_count += 1
            year = int(rows[0][3][:4])
            yearly[year]["database_groups"] += 1
            try:
                key, row, has_tc, has_mod06 = expected_record(rows, root)
            except ValueError as exc:
                yearly[year]["DATABASE_GROUP_INVALID"] += 1
                database_errors.append({"gid": gid, "year": year, "problem": str(exc)})
                continue
            matches = matching_items(primary, key, has_tc)
            valid = [i for i in matches if not i[3]]
            if valid:
                yearly[year]["PRESENT"] += 1
                yearly[year]["PRESENT_" + key[2]] += 1
                if not any(i[0] == key for i in valid):
                    yearly[year]["MO_ORBIT_NAME_ALIAS"] += 1
            else:
                row["_key"] = key
                row["_has_tc"] = has_tc
                row["_has_mod06"] = has_mod06
                row["_invalid"] = [i[3] for i in matches]
                row["stored_matches"] = ";".join(i[1] for i in matching_items(supplemental, key, has_tc) if not i[3])
                others = [i for i in primary.get(key[0], []) if not i[3] and i[0][2] != key[2]
                          and (i[0][2] == "MO" or i[0][1] == key[1])]
                row["other_product_types"] = ";".join(sorted({i[0][2] for i in others}))
                row["other_stored_paths"] = ";".join(i[1] for i in others)
                missing.append(row)
            if group_count % 100000 == 0:
                print("DATABASE", group_count, "groups checked;", len(missing), "absent/invalid candidates", flush=True)
    print("DATABASE_COMPLETE", group_count, "groups;", len(missing), "candidates", flush=True)
    # Save candidates before interpreting histories, so the filesystem result is
    # independently reviewable and retained if reconciliation is interrupted.
    write_tsv(dest / "missing_before_run_reconciliation.tsv", missing)
    wanted_by_year = defaultdict(set)
    for row in missing:
        wanted_by_year[row["year"]].add(row["gid"])
    latest, manifests, malformed = {}, {}, []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(scan_run_records, (year, root, wanted)) for year, wanted in wanted_by_year.items()]
        for future in as_completed(futures):
            year, events, runs, errors = future.result()
            latest.update(events)
            manifests[year] = runs
            malformed.extend(errors)
            print("MANIFESTS", year, len(runs), "runs;", len(events), "candidates with records", flush=True)
    unresolved, skips, elsewhere = [], [], []
    for row in missing:
        event = latest.get(row["gid"], {})
        # A repaired catalog may retain its group ID but change the pairing.
        # Do not suppress a missing current pair using a different old pair.
        ekey = product_key(os.path.basename(event.get("expected_output") or ""))
        compatible = ekey is not None and (ekey == row["_key"] or
                      (not row["_has_tc"] and ekey[0] == row["_key"][0] and ekey[2] == row["_key"][2]))
        row.update(last_status=event.get("status", "NO_RUN_RECORD") if compatible else "NO_MATCHING_RUN_RECORD",
                   last_run=event.get("run", ""), last_event_at=event.get("timestamp", event.get("planned_at", "")),
                   recorded_status=event.get("status", "NO_RUN_RECORD"), recorded_product_type=ekey[2] if ekey else "",
                   recorded_expected_output=event.get("expected_output", ""),
                   message=event.get("message", ""))
        status = row["last_status"]
        if not row["_has_tc"] and not row["_has_mod06"]:
            row["status"] = "NO_CLOUD_INPUTS_IN_CATALOG"
        elif row["_invalid"]:
            row["status"] = "INVALID_STORED_OUTPUT"
        elif status == "SKIPPED_NO_OVERLAP":
            row["status"] = "EXPECTED_NO_OVERLAP"
        elif row["stored_matches"]:
            row["status"] = "STORED_OUTSIDE_PRODUCTION"
        elif row["product_type"] == "MM" and row["other_stored_paths"]:
            row["status"] = "PARTIAL_PRODUCT_ONLY"
        elif status.startswith("FAILED") or status == "SKIPPED_VALUE_ERROR":
            row["status"] = "CONFIRMED_FAILED"
        elif status == "SAVED":
            row["status"] = "SAVED_OUTPUT_NOW_MISSING"
        elif status in {"STARTED", "NOT_STARTED"}:
            row["status"] = "INCOMPLETE_RUN"
        else:
            row["status"] = "MISSING_WITHOUT_FAILURE_RECORD"
        yearly[row["year"]][row["status"]] += 1
        if row["status"] == "EXPECTED_NO_OVERLAP":
            skips.append(row)
        elif row["status"] == "STORED_OUTSIDE_PRODUCTION":
            elsewhere.append(row)
        else:
            check_inputs(row)
            unresolved.append(row)
    for name, rows in [("unresolved_groups.tsv", unresolved), ("expected_no_overlap.tsv", skips),
                       ("stored_outside_production.tsv", elsewhere),
                       ("production_gaps.tsv", sorted(unresolved + elsewhere, key=lambda r: r["gid"]))]:
        write_tsv(dest / name, rows)
    write_tsv(dest / "database_group_problems.tsv", database_errors, ["gid", "year", "problem"])
    (dest / "malformed_manifest_records.txt").write_text("".join(p + "\n" for p in malformed))
    def id_list(path, rows):
        ids = sorted({r["modis_id"] for r in rows})
        path.write_text("".join(i + "\n" for i in ids))
        return len(ids)
    rerunnable = [r for r in unresolved if r["status"] != "NO_CLOUD_INPUTS_IN_CATALOG"]
    unique_rerun = id_list(dest / "rerun_candidate_ids.txt", rerunnable)
    id_list(dest / "production_gap_ids.txt", unresolved + elsewhere)
    id_list(dest / "stored_outside_production_ids.txt", elsewhere)
    id_list(dest / "partial_product_ids.txt", [r for r in unresolved if r["status"] == "PARTIAL_PRODUCT_ONLY"])
    id_list(dest / "unexplained_missing_ids.txt", [r for r in unresolved if r["status"] == "MISSING_WITHOUT_FAILURE_RECORD"])
    missing_inputs = sorted({path for row in unresolved for path in row["missing_required_inputs"].split(";") if path})
    (dest / "missing_required_input_paths.txt").write_text("".join(path + "\n" for path in missing_inputs))
    input_read_errors = []
    for row in unresolved:
        for product, role in [("MOD06", "mod06_file"), ("MOD03", "mod03_file")]:
            if product + "_READ_FAILED" in row["message"]:
                input_read_errors.append({"modis_id": row["modis_id"], "input_product": product,
                                          "input_path": row[role], "current_check": row["input_check"]})
    write_tsv(dest / "input_read_errors.tsv", input_read_errors, ["modis_id", "input_product", "input_path", "current_check"])
    for year in range(first_year, last_year + 1):
        rows = [r for r in rerunnable if r["year"] == year]
        if rows:
            id_list(dest / f"{year}_rerun_candidate_ids.txt", rows)
    id_list(dest / "confirmed_failure_ids.txt", [r for r in unresolved if r["status"] == "CONFIRMED_FAILED"])
    totals = sum(yearly.values(), Counter())
    db_after = args.database.stat()
    summary = {"started_at": started, "finished_at": datetime.now(ZoneInfo("America/Chicago")).isoformat(timespec="seconds"),
               "database": str(args.database), "start_date": args.start_date, "end_date": args.end_date,
               "database_unchanged_during_audit": (db_before.st_size, db_before.st_mtime_ns) == (db_after.st_size, db_after.st_mtime_ns),
               "output_root": root, "totals": dict(totals), "years": {str(y): dict(c) for y, c in sorted(yearly.items())},
               "inventory": inventory_stats, "unresolved_groups": len(unresolved), "unique_rerun_candidate_ids": unique_rerun,
               "production_gap_groups": len(unresolved) + len(elsewhere), "missing_required_input_paths": missing_inputs,
               "expected_no_overlap_groups": len(skips), "stored_outside_production_groups": len(elsewhere),
               "invalid_stored_outputs": len(bad_outputs), "database_group_problems": len(database_errors),
               "malformed_manifest_records": len(malformed), "manifests_used_for_missing_groups": manifests,
               "validation": "All database groups; sizes of all stored production products; signature checks on all products smaller than 64 KiB and eight sampled products per directory; latest compatible events for missing groups. No exhaustive file-content, NetCDF metadata, pixel-array, decoder-version, or scientific-validity audit."}
    (dest / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    yearly_fields = ["year", "database_groups", "PRESENT", "EXPECTED_NO_OVERLAP", "CONFIRMED_FAILED", "SAVED_OUTPUT_NOW_MISSING",
                     "INCOMPLETE_RUN", "PARTIAL_PRODUCT_ONLY", "MISSING_WITHOUT_FAILURE_RECORD", "INVALID_STORED_OUTPUT", "STORED_OUTSIDE_PRODUCTION", "NO_CLOUD_INPUTS_IN_CATALOG"]
    write_tsv(dest / "yearly_summary.tsv", [{"year": y, **c} for y, c in sorted(yearly.items())], yearly_fields)
    lines = ["# MMCTH archive completion audit", "", f"Period: {args.start_date} through {args.end_date}, inclusive.",
             f"Finished: {summary['finished_at']}", "", f"Database: `{args.database}` (opened read-only).", "",
             f"Compared **{group_count:,} database groups** with products actually stored in `{root}`.", "",
             f"- Present production products: **{totals['PRESENT']:,} groups**.",
             f"- Expected no-overlap skips supported by the latest matching records: **{len(skips):,} groups**.",
             f"- Unresolved missing/invalid groups: **{len(unresolved):,}** ({unique_rerun:,} unique rerun-candidate MODIS IDs).",
             f"- Matching products stored outside production: **{len(elsewhere):,} groups**.", "",
             "| Year | Database groups | Present | No overlap | Unresolved | Elsewhere |", "|---|---:|---:|---:|---:|---:|"]
    for y, c in sorted(yearly.items()):
        n = c["database_groups"] - c["PRESENT"] - c["EXPECTED_NO_OVERLAP"] - c["STORED_OUTSIDE_PRODUCTION"]
        lines.append(f"| {y} | {c['database_groups']:,} | {c['PRESENT']:,} | {c['EXPECTED_NO_OVERLAP']:,} | {n:,} | {c['STORED_OUTSIDE_PRODUCTION']:,} |")
    lines += ["", "## Unresolved reasons", ""]
    for name, count in sorted(Counter(r["status"] for r in unresolved).items()):
        lines.append(f"- {name}: {count:,}")
    lines += ["", "`unresolved_groups.tsv` includes exact input paths, expected output paths, newest run status, messages, and current required-input checks.",
              "`production_gaps.tsv` and `production_gap_ids.txt` include all missing production products, including those with supplemental copies.",
              "`rerun_candidate_ids.txt` and yearly lists include recoverable missing/invalid products. `confirmed_failure_ids.txt` contains only failures documented by a matching run record.",
              "`partial_product_ids.txt` identifies expected MM products for which only MO or MI products were saved. `unexplained_missing_ids.txt` contains missing current pairings without a compatible failure or no-overlap record.",
              "`stored_outside_production.tsv` identifies products found only in supplemental output directories; these may need restoring rather than reprocessing.", "",
              "## Input checks", "",
              f"Missing required source paths: {len(missing_inputs)} unique files across {sum(bool(r['missing_required_inputs']) for r in unresolved)} groups. See `missing_required_input_paths.txt`.",
              f"Recorded MOD03/MOD06 read errors: {len(input_read_errors)} groups. Their input paths and current HDF open results are in `input_read_errors.tsv`.", "",
              "## Matching and limits", "",
              "Both historical product filename layouts are recognized. MM and MI require the database TC orbit. MO names without a TC input are matched by MODIS timestamp and product type because historical date-range jobs used a chunk-level orbit in their filenames.",
              "Only a compatible latest SKIPPED_NO_OVERLAP event excludes an absent output. A missing file without such a record remains a candidate; the audit does not assume processing failed or spatial overlap exists when that evidence is unavailable.",
              summary["validation"], "", f"Database unchanged during audit: {summary['database_unchanged_during_audit']}.",
              "Source database and production files were not modified. No production jobs were submitted.", ""]
    (dest / "summary.md").write_text("\n".join(lines))
    print("AUDIT_COMPLETE", json.dumps({"totals": dict(totals), "unresolved": len(unresolved), "unique_rerun_ids": unique_rerun, "report": str(dest / 'summary.md')}), flush=True)


if __name__ == "__main__":
    main()
