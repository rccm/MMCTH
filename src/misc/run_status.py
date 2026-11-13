"""Small, append-only run manifests for long MPI production jobs."""

from __future__ import annotations

import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path


TERMINAL_STATUSES = {
    "SAVED",
    "SKIPPED_NO_OVERLAP",
    "FAILED_NO_OUTPUT",
    "FAILED_EXCEPTION",
    "SKIPPED_VALUE_ERROR",
}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def make_run_id() -> str:
    """Return a filesystem-safe ID; call on rank 0 and broadcast for manual MPI."""
    array_job = os.environ.get("SLURM_ARRAY_JOB_ID")
    array_task = os.environ.get("SLURM_ARRAY_TASK_ID")
    if array_job:
        suffix = f"_{array_task}" if array_task is not None else ""
        return f"slurm_{array_job}{suffix}"

    slurm_job = os.environ.get("SLURM_JOB_ID")
    if slurm_job:
        return f"slurm_{slurm_job}"

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    return f"manual_{stamp}_{os.getpid()}"


def _is_real_input(path: str, check_exists: bool = True) -> bool:
    if os.path.basename(str(path)) == "X":
        return False
    return os.path.exists(path) if check_exists else True


def expected_output_path(
    files,
    orbit_number,
    output_dir: str,
    *,
    check_input_exists: bool = True,
) -> str | None:
    """Reproduce MainProcessor's production output name without opening inputs."""
    if not files:
        return None

    modis_match = re.search(r"(A\d{7}\.\d{4})", str(files[0]))
    if modis_match is None:
        return None

    modis_id = modis_match.group(1)
    date = datetime.strptime(modis_id[1:8], "%Y%j")
    mmdd = date.strftime("%m%d")

    has_mod06 = len(files) > 1 and _is_real_input(files[1], check_input_exists)
    has_tc = len(files) > 3 and _is_real_input(files[3], check_input_exists)
    if has_mod06 and has_tc:
        file_flag = "MM"
    elif not has_mod06:
        file_flag = "MI"
    elif not has_tc:
        file_flag = "MO"
    else:
        file_flag = "NW"

    tc_orbit = re.search(r"(O\d{6})", str(files[3])) if len(files) > 3 else None
    if tc_orbit is not None:
        orbit = tc_orbit.group(1)
    elif orbit_number is not None:
        orbit = "O" + str(orbit_number).replace("O", "").zfill(6)
    else:
        return None

    filename = f"MODMISR_L2_CP_V00_061_{modis_id}_{mmdd}_{orbit}_{file_flag}.nc"
    return str(Path(output_dir).expanduser().resolve() / filename)


def work_item_record(
    gid,
    files,
    orbit_number,
    output_dir: str,
    *,
    check_input_exists: bool = True,
) -> dict:
    modis_name = os.path.basename(str(files[0])) if files else ""
    modis_match = re.search(r"(\.A\d{7}\.\d{4}\.)", modis_name)
    tc_name = os.path.basename(str(files[3])) if len(files) > 3 else ""
    return {
        "gid": int(gid) if isinstance(gid, int) or str(gid).isdigit() else str(gid),
        "modis_id": modis_match.group(1) if modis_match else "",
        "modis_file": modis_name,
        "misr_file": tc_name,
        "orbit": orbit_number,
        "expected_output": expected_output_path(
            files,
            orbit_number,
            output_dir,
            check_input_exists=check_input_exists,
        ),
    }


class RunStatusWriter:
    """Write a per-rank plan and line-buffered events without shared-file races."""

    def __init__(self, output_dir: str, run_id: str, rank: int):
        self.rank = int(rank)
        self.run_dir = Path(output_dir).expanduser().resolve() / "_run_logs" / run_id
        self.run_dir.mkdir(parents=True, exist_ok=True)
        self.plan_path = self.run_dir / f"planned_rank_{self.rank:04d}.jsonl"
        self.events_path = self.run_dir / f"events_rank_{self.rank:04d}.jsonl"
        self._events = self.events_path.open("a", encoding="utf-8", buffering=1)

    def write_plan(self, work_items, output_dir: str) -> None:
        temp_path = self.plan_path.with_name(
            f".{self.plan_path.name}.tmp.{os.getpid()}"
        )
        with temp_path.open("w", encoding="utf-8") as stream:
            for gid, files, orbit_number in work_items:
                record = work_item_record(
                    gid,
                    files,
                    orbit_number,
                    output_dir,
                    check_input_exists=False,
                )
                record.update({"rank": self.rank, "planned_at": utc_now()})
                stream.write(json.dumps(record, sort_keys=True) + "\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temp_path, self.plan_path)

    def event(self, base_record: dict, status: str, **extra) -> None:
        record = dict(base_record)
        record.update(
            {
                "rank": self.rank,
                "status": status,
                "timestamp": utc_now(),
            }
        )
        record.update(extra)
        self._events.write(json.dumps(record, sort_keys=True) + "\n")
        self._events.flush()

    def close(self) -> None:
        if not self._events.closed:
            self._events.close()
