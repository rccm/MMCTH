#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Quick‑and‑clean processing of multilayer MODIS–MISR cloud pixels.

* Reads a single combined **MODMISR_L2_CP** NetCDF file (path set in
  ``INPUT_FILE`` below).
* Selects multilayer pixels (MISR–MODIS CTH difference > 1 km, ice‑phase
  MODIS cloud, valid BT/radiance).
* Converts MISR cloud‑top height (CTH) → cloud‑top pressure (CTP) with the
  supplied Fortran routine.
* Runs the CO₂‑slicing retrieval for each selected pixel via
  ``src.pixel.modis`` bindings.
* If ``DEBUG_PLOT`` is ``True`` a quick‑look WV/T profile for the **first**
  pixel is saved for sanity‑checking.

Run the script directly – no command‑line arguments are needed.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import xarray as xr

# -----------------------------------------------------------------------------
# Matplotlib – headless backend *before* importing pyplot
# -----------------------------------------------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402  pylint: disable=wrong-import-position

# -----------------------------------------------------------------------------
# Local Fortran/Cython wrappers
# -----------------------------------------------------------------------------
from misr_cth_to_pressure import height_to_log_pressure  # converts CTH → CTP
from src.pixel import modis  # CO₂‑slicing algorithm

# -----------------------------------------------------------------------------
# USER SETTINGS – adjust as needed
# -----------------------------------------------------------------------------
INPUT_FILE = Path(
    "/data/keeling/a/gzhao1/f/mmcth/output/"
    "MODMISR_L2_CP_A2016074.0330_O086367_debug.nc"
)
DEBUG_PLOT = False  # True ➜ save WV/T profiles for the first pixel

# -----------------------------------------------------------------------------
# STANDARD ATMOSPHERIC PROFILES (ECMWF 101‑level)
# -----------------------------------------------------------------------------
PSTD = np.array([...], dtype=np.float64)  # truncated for brevity – identical to ECMWF table
TSTD = np.array([...], dtype=np.float32)  # truncated – standard temperature profile
WSTD = np.array([...], dtype=np.float32)  # truncated – standard specific humidity
PLEV = PSTD  # alias for readability

# (NOTE: Full numeric arrays are unchanged – see previous revision.)

# -----------------------------------------------------------------------------
# Helper – quick‑look plotting
# -----------------------------------------------------------------------------

def _plot_profile(wprof: np.ndarray, tprof: np.ndarray, psfc: float,
                  misr_ctp: float, modis_ctp: float, vza: float, stem: str) -> None:
    """Save water‑vapour & temperature profiles on log‑pressure axes."""

    def _mark(ax, y: float, label: str) -> None:
        ax.axhline(y, ls="--", lw=0.8, color="k")
        ax.text(ax.get_xlim()[1] * 0.97, y * 1.03,
                f"{label}: {y:.0f} hPa", ha="right", va="bottom", fontsize=8)

    # Water‑vapour profile ----------------------------------------------------
    fig, ax = plt.subplots(figsize=(5, 4))
    valid = wprof > 0
    ax.plot(wprof[valid], PLEV[valid], label="ERA5+Std")
    ax.plot(WSTD, PLEV, label="Std")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.invert_yaxis()
    ax.set_xlabel("q  [g kg⁻¹]"); ax.set_ylabel("Pressure [hPa]")
    ax.grid(True, ls=":", which="both")
    _mark(ax, psfc, "Sfc"); _mark(ax, misr_ctp, "MISR"); _mark(ax, modis_ctp, "MODIS")
    ax.legend(); fig.tight_layout(); fig.savefig(f"wv_{stem}.png", dpi=300); plt.close(fig)

    # Temperature profile ----------------------------------------------------
    fig, ax = plt.subplots(figsize=(5, 4))
    valid = tprof > 0
    ax.plot(tprof[valid], PSTD[valid], label="ERA5")
    ax.plot(TSTD, PSTD, label="Std")
    ax.set_yscale("log"); ax.invert_yaxis()
    ax.set_xlabel("T  [K]"); ax.set_ylabel("Pressure [hPa]")
    ax.grid(True, ls=":", which="both")
    _mark(ax, psfc, "Sfc"); _mark(ax, misr_ctp, "MISR"); _mark(ax, modis_ctp, "MODIS")
    ax.legend(); fig.tight_layout(); fig.savefig(f"tmp_{stem}.png", dpi=300); plt.close(fig)

# -----------------------------------------------------------------------------
# Pixel selection & processing
# -----------------------------------------------------------------------------

def _select_pixels(ds: xr.Dataset) -> tuple[np.ndarray, np.ndarray]:
    """Return indices (iy, ix) of pixels satisfying multilayer criteria."""
    misr_cth = ds["misr_cth"].where(ds["misr_cth"] > -888).values
    mod06 = ds["mod06_product"].values
    mod_cth, mod_ctp, mod_cphase = mod06[3], mod06[0], mod06[2]

    height_diff = mod_cth - misr_cth
    bt_rad = ds["bands_BT"].values
    bt_ok = np.all(bt_rad[:5] > 0, axis=0)  # first 5 indexes are BT

    mask = (height_diff > 1000) & (mod_cphase == 3) & bt_ok & (misr_cth > 0)
    return np.where(mask)


def main() -> None:
    print(f"Opening {INPUT_FILE}")
    ds = xr.open_dataset(INPUT_FILE)

    iy, ix = _select_pixels(ds)
    if iy.size == 0:
        print("No multilayer pixels found – exiting.")
        return
    print(f"Selected {iy.size} multilayer pixels")

    # Gather arrays ----------------------------------------------------------
    sp_prof = ds["mixingratio_exp"].values[:, iy, ix]
    t_prof  = ds["temperature_exp"].values[:, iy, ix]
    z_prof  = ds["geopotential_exp"].values[:, iy, ix]

    misr_cth = ds["misr_cth"].values[iy, ix]
    misr_ctp = height_to_log_pressure(z_prof, misr_cth / 1000.0,
                                      z_prof.shape[1], z_prof.shape[2])

    bt_rad = ds["bands_BT"].values
    trad = bt_rad[5:, iy, ix]        # radiance (not BT) bands for retrieval

    landsea = ds["landsea"].values[iy, ix]
    sst     = ds["sst"].values[iy, ix]
    skint   = ds["skint"].values[iy, ix]
    surftmp = np.where(landsea == 1, skint, sst)

    surfp = ds["surface_pressure"].values[iy, ix] / 100.0  # Pa→hPa
    mslp  = ds["msp"].values[iy, ix] / 100.0
    vza   = ds["vza"].values[iy, ix]
    lat   = ds["latitude"].values[iy, ix]
    lon   = ds["longitude"].values[iy, ix]

    # Example date (year, month, day, hour)
    met_date = np.array([2014, 12, 12, 0], dtype="int32")

    print("Running CO₂‑slicing …")
    ctp, cee, copt, pmask = modis.process_selected_pixels(
        sp_prof,
        t_prof,
        surfp,
        mslp,
        surftmp,
        vza,
        trad,
        lat,
        lon,
        landsea,
        misr_ctp,
        met_date,
        iy.size,
    )

    valid_ctp = np.count_nonzero(ctp > 0)
    print(f"Finished. Valid CTPs retrieved: {valid_ctp}")

    # Optional quick‑look ----------------------------------------------------
    if DEBUG_PLOT:
        print("Saving debug WV/T profile for first pixel …")
        _plot_profile(
            sp_prof[:, 0], t_prof[:, 0], surfp[0],
            misr_ctp[0], ds["mod06_product"].values[0, iy[0], ix[0]],
            vza[0], stem="pix0"
        )
        print("Profiles saved (wv_pix0.png, tmp_pix0.png)")


if __name__ == "__main__":
    main()
