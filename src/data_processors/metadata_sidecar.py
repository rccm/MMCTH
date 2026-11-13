"""Optional CMR JSON creation immediately after a production NetCDF is closed."""

import importlib.util
import json
import logging
import os
import tempfile
from functools import lru_cache
from pathlib import Path


@lru_cache(maxsize=1)
def _load_generator(script):
    """Load the existing generator once in each producing MPI process."""
    spec = importlib.util.spec_from_file_location("mmcth_inline_metadata", script)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load metadata generator: {script}")
    generator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(generator)
    return generator


def create_json_sidecar(output_path, logger=None):
    """Write adjacent JSON when enabled; metadata failures leave NetCDF intact.

    CREATE_JSON_METADATA defaults to off outside the chunk submission script.
    Each newly saved NetCDF gets fresh metadata, even if an older sidecar exists.
    No additional processes are started and the existing generator is unchanged.
    """
    if os.environ.get("CREATE_JSON_METADATA", "0").strip().lower() not in {
        "1", "true", "yes", "on",
    }:
        return None

    log = logger or logging.getLogger(__name__)
    temporary = None
    try:
        source = Path(output_path).resolve()
        default_script = Path(__file__).resolve().parents[2] / "metadata/creat_json.py"
        script = Path(os.environ.get("METADATA_SCRIPT") or default_script).resolve()
        generator = _load_generator(str(script))
        before = source.stat()
        metadata, mod03 = generator.build_ummg(
            str(source),
            os.environ.get("MOD03_ROOT") or "/data/gdi/satellite/TerraDataArchive/MODIS/MOD03",
            os.environ.get("COLLECTION_SHORTNAME") or "MODMISR_L2_CP",
            os.environ.get("COLLECTION_VERSION") or "V00",
        )
        after = source.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise RuntimeError("Source NetCDF changed while generating JSON; retry metadata")

        sidecar = source.with_name(source.name + ".cmr.json")
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=source.parent,
            prefix=f".{sidecar.name}.tmp.", delete=False,
        ) as stream:
            temporary = Path(stream.name)
            json.dump(metadata, stream, indent=2)
            stream.write("\n")
            os.fchmod(stream.fileno(), before.st_mode & 0o666)
        os.replace(temporary, sidecar)
        temporary = None
        log.info("Created JSON metadata: %s", sidecar)
        if mod03 is None:
            log.warning("JSON metadata used fallback sensing time/day-night for %s", source)
        return sidecar
    except Exception:
        log.exception("JSON metadata failed for %s; NetCDF remains saved", output_path)
        return None
    finally:
        if temporary is not None:
            try:
                temporary.unlink(missing_ok=True)
            except OSError:
                log.warning("Could not remove temporary JSON: %s", temporary, exc_info=True)
