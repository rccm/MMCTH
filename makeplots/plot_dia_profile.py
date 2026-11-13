#!/usr/bin/env python3
"""
Plot vertical diagnostic meteorological profiles for one MM diagnostic pixel.

Selection options:

1. By swath row/column:

   python plot_diag_profile.py file.nc \
     --swath-row 1200 \
     --swath-col 300 \
     -o profile.png

2. By diagnostic-pixel index:

   python plot_diag_profile.py file.nc \
     --pixel-index 123 \
     -o profile.png

Plot layout:

  lower x-axis: TemperatureProfile, K, red
  upper x-axis: WaterVaporProfile, g kg-1, blue
  left y-axis: pressure, hPa, log scale
  right y-axis: height, km

Horizontal lines shown:

  SurfacePressure
  MISR CTP / CTH
  MODIS CTP / CTH
  OneLayer_CTP
  FinalTwoLayer_CloudTopPressure, or TwoLayer_CandidateCTP if no final two-layer solution exists
  CO2_TropopauseHeight / diag_ltrp
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional

import h5py
import numpy as np
import matplotlib.pyplot as plt


MISSING_VALUE = -9999.0
MISSING_THRESHOLD = -900.0


def clean_array(a: np.ndarray) -> np.ndarray:
    """Convert to float and replace invalid/missing sentinels with NaN."""
    out = np.asarray(a, dtype=float)
    out[~np.isfinite(out)] = np.nan
    out[np.isclose(out, MISSING_VALUE) | (out <= MISSING_THRESHOLD)] = np.nan
    return out


def clean_scalar(x) -> Optional[float]:
    """Convert scalar to float; return None if invalid or missing."""
    try:
        v = float(x)
    except Exception:
        return None

    if not np.isfinite(v):
        return None

    if np.isclose(v, MISSING_VALUE) or v <= MISSING_THRESHOLD:
        return None

    return v


def scalar_at(f: h5py.File, name: str, i: int) -> Optional[float]:
    """Read scalar value at diagnostic pixel index i."""
    if name not in f:
        return None

    try:
        return clean_scalar(f[name][i])
    except Exception:
        return None


def array_at(f: h5py.File, name: str, i: int) -> Optional[np.ndarray]:
    """Read 1-D profile at diagnostic pixel index i."""
    if name not in f:
        return None

    try:
        return clean_array(f[name][i, :])
    except Exception:
        return None


def find_diagnostic_pixel_by_row_col(
    f: h5py.File,
    swath_row: int,
    swath_col: int,
) -> int:
    """Find diagnostic_pixel index matching swath_row and swath_col."""
    rows = np.asarray(f["swath_row"][:])
    cols = np.asarray(f["swath_col"][:])

    matches = np.where((rows == swath_row) & (cols == swath_col))[0]

    if matches.size == 0:
        raise ValueError(
            f"No diagnostic pixel found for swath_row={swath_row}, "
            f"swath_col={swath_col}."
        )

    if matches.size > 1:
        print(
            f"Warning: found {matches.size} matching diagnostic pixels. "
            f"Using the first one: diagnostic_pixel index {matches[0]}"
        )

    return int(matches[0])


def choose_diagnostic_pixel(
    f: h5py.File,
    pixel_index: Optional[int],
    swath_row: Optional[int],
    swath_col: Optional[int],
) -> int:
    """
    Choose diagnostic-pixel index either by direct pixel index or by swath row/col.
    """
    n_pix = f["diagnostic_pixel"].shape[0]

    if pixel_index is not None:
        if pixel_index < 0 or pixel_index >= n_pix:
            raise ValueError(
                f"--pixel-index {pixel_index} is outside valid range "
                f"0 to {n_pix - 1}."
            )
        return int(pixel_index)

    if swath_row is not None and swath_col is not None:
        return find_diagnostic_pixel_by_row_col(f, swath_row, swath_col)

    raise ValueError(
        "You must provide either --pixel-index N, "
        "or both --swath-row ROW and --swath-col COL."
    )


def pressure_from_height(
    target_height_km: Optional[float],
    pressure_hpa: np.ndarray,
    height_km: np.ndarray,
) -> Optional[float]:
    """Convert height in km to pressure in hPa using local HeightProfile."""
    if target_height_km is None or not np.isfinite(target_height_km):
        return None

    valid = (
        np.isfinite(height_km)
        & np.isfinite(pressure_hpa)
        & (pressure_hpa > 0)
    )

    if valid.sum() < 2:
        return None

    h = height_km[valid]
    p = pressure_hpa[valid]

    order = np.argsort(h)
    h = h[order]
    p = p[order]

    h_unique, unique_idx = np.unique(h, return_index=True)
    p_unique = p[unique_idx]

    if target_height_km < np.nanmin(h_unique) or target_height_km > np.nanmax(h_unique):
        return None

    return float(np.interp(target_height_km, h_unique, p_unique))


def height_from_pressure(
    target_pressure_hpa: Optional[float],
    pressure_hpa: np.ndarray,
    height_km: np.ndarray,
) -> Optional[float]:
    """Convert pressure in hPa to height in km using local HeightProfile."""
    if target_pressure_hpa is None or not np.isfinite(target_pressure_hpa):
        return None

    valid = (
        np.isfinite(height_km)
        & np.isfinite(pressure_hpa)
        & (pressure_hpa > 0)
    )

    if valid.sum() < 2:
        return None

    p = pressure_hpa[valid]
    h = height_km[valid]

    order = np.argsort(p)
    p = p[order]
    h = h[order]

    p_unique, unique_idx = np.unique(p, return_index=True)
    h_unique = h[unique_idx]

    if target_pressure_hpa < np.nanmin(p_unique) or target_pressure_hpa > np.nanmax(p_unique):
        return None

    return float(np.interp(target_pressure_hpa, p_unique, h_unique))


def pressure_from_fortran_level(
    level_index: Optional[float],
    pressure_hpa: np.ndarray,
) -> Optional[float]:
    """Read pressure at a Fortran 1-based level index."""
    if level_index is None or not np.isfinite(level_index):
        return None

    rounded = int(round(level_index))
    if not np.isclose(level_index, rounded):
        return None

    zero_based = rounded - 1
    if zero_based < 0 or zero_based >= pressure_hpa.size:
        return None

    pressure = float(pressure_hpa[zero_based])
    if not np.isfinite(pressure) or pressure <= 0:
        return None

    return pressure


def height_m_to_km(height_m: Optional[float]) -> Optional[float]:
    """Convert a cloud-top height scalar in meters to kilometers."""
    if height_m is None or not np.isfinite(height_m):
        return None

    return height_m / 1000.0


def add_horizontal_marker_inside(
    ax,
    y_pressure: Optional[float],
    label: str,
    color: str = "0.25",
    linestyle: str = "--",
    linewidth: float = 1.0,
    xfrac: float = 0.04,
):
    """Add horizontal line with label inside plot."""
    if y_pressure is None or not np.isfinite(y_pressure) or y_pressure <= 0:
        return

    y0, y1 = ax.get_ylim()
    p_min = min(y0, y1)
    p_max = max(y0, y1)

    if y_pressure < p_min or y_pressure > p_max:
        return

    ax.axhline(
        y_pressure,
        color=color,
        linestyle=linestyle,
        linewidth=linewidth,
        alpha=0.85,
        zorder=1,
    )

    ax.text(
        xfrac,
        y_pressure,
        label,
        transform=ax.get_yaxis_transform(),
        va="center",
        ha="left",
        fontsize=8,
        color=color,
        bbox=dict(
            facecolor="white",
            edgecolor="none",
            alpha=0.70,
            pad=1.5,
        ),
        zorder=5,
    )


def collect_required_markers(
    f: h5py.File,
    i: int,
    pressure_hpa: np.ndarray,
    height_km_profile: np.ndarray,
):
    """
    Collect only requested horizontal line markers.

    Returns list of:

        label, pressure_hpa
    """
    markers = []

    # Surface pressure.
    sp = scalar_at(f, "SurfacePressure", i)
    if sp is not None:
        markers.append((f"Surface P: {sp:.1f} hPa", sp, "black", "-", 1.3))

    # MISR CTP and CTH.
    misr_ctp = scalar_at(f, "MISR_CloudTopPressure_Used", i)
    if misr_ctp is None:
        misr_ctp = scalar_at(f, "MISR_CloudTopPressure_Input", i)

    misr_cth_m = scalar_at(f, "MISR_CloudTopHeight", i)
    misr_cth_km = None if misr_cth_m is None else misr_cth_m / 1000.0

    if misr_ctp is None and misr_cth_km is not None:
        misr_ctp = pressure_from_height(misr_cth_km, pressure_hpa, height_km_profile)

    if misr_ctp is not None:
        if misr_cth_km is not None:
            label = f"MISR: {misr_ctp:.1f} hPa, {misr_cth_km:.2f} km"
        else:
            label = f"MISR CTP: {misr_ctp:.1f} hPa"
        markers.append((label, misr_ctp, "0.20", "--", 1.0))

    # MODIS CTP and CTH.
    modis_ctp = scalar_at(f, "MODIS_CloudTopPressure", i)
    modis_cth_m = scalar_at(f, "MODIS_CloudTopHeight", i)
    modis_cth_km = None if modis_cth_m is None else modis_cth_m / 1000.0

    if modis_ctp is None and modis_cth_km is not None:
        modis_ctp = pressure_from_height(modis_cth_km, pressure_hpa, height_km_profile)

    if modis_ctp is not None:
        if modis_cth_km is not None:
            label = f"MODIS: {modis_ctp:.1f} hPa, {modis_cth_km:.2f} km"
        else:
            label = f"MODIS CTP: {modis_ctp:.1f} hPa"
        markers.append((label, modis_ctp, "0.35", "--", 1.0))

    # One-layer CTP only.
    one_ctp = scalar_at(f, "OneLayer_CTP", i)
    if one_ctp is not None:
        markers.append((f"One-layer CTP: {one_ctp:.1f} hPa", one_ctp, "0.45", "-.", 1.0))

    # Final accepted two-layer result only.  Older diagnostic files do not have
    # FinalTwoLayer_* fields, so fall back only when the final MM flag says the
    # two-layer solution was accepted.
    final_two_ctp = scalar_at(f, "FinalTwoLayer_CloudTopPressure", i)
    final_two_cth_km = height_m_to_km(scalar_at(f, "FinalTwoLayer_CloudTopHeight", i))

    if final_two_ctp is None:
        mm_flag = scalar_at(f, "MM_Flag", i)
        if mm_flag is not None and int(round(mm_flag)) in (8, 9):
            final_two_ctp = scalar_at(f, "TwoLayer_CandidateCTP", i)
            final_two_cth_km = scalar_at(f, "TwoLayer_CandidateCTH", i)

    if final_two_ctp is not None:
        if final_two_cth_km is not None:
            label = f"Final two-layer: {final_two_ctp:.1f} hPa, {final_two_cth_km:.2f} km"
        else:
            label = f"Final two-layer CTP: {final_two_ctp:.1f} hPa"
        markers.append((label, final_two_ctp, "0.55", "-.", 1.0))
    else:
        candidate_ctp = scalar_at(f, "TwoLayer_CandidateCTP", i)
        candidate_cth_km = scalar_at(f, "TwoLayer_CandidateCTH", i)
        if candidate_ctp is not None:
            if candidate_cth_km is not None:
                label = f"Two-layer candidate: {candidate_ctp:.1f} hPa, {candidate_cth_km:.2f} km"
            else:
                label = f"Two-layer candidate CTP: {candidate_ctp:.1f} hPa"
            markers.append((label, candidate_ctp, "0.55", "-.", 1.0))

    # CO2 Fortran-selected tropopause level.
    co2_ltrp = scalar_at(f, "diag_ltrp", i)
    co2_tropopause_height = scalar_at(f, "CO2_TropopauseHeight", i)
    co2_tropopause_pressure = pressure_from_fortran_level(co2_ltrp, pressure_hpa)

    if co2_tropopause_pressure is None and co2_tropopause_height is not None:
        co2_tropopause_pressure = pressure_from_height(
            co2_tropopause_height,
            pressure_hpa,
            height_km_profile,
        )

    if co2_tropopause_pressure is not None:
        label = f"tropopause:{co2_tropopause_pressure:.1f}hPa"
        if co2_tropopause_height is not None:
            label += f", {co2_tropopause_height:.2f} km"
        if co2_ltrp is not None:
            label += f", L={int(round(co2_ltrp))}"
        markers.append((label, co2_tropopause_pressure, "purple", ":", 1.2))

    return markers


def set_exact_temp_limits(ax, temperature: np.ndarray):
    """Set lower x-axis exactly to min/max temperature, excluding NaN/-9999."""
    temperature[temperature< 200] =np.nan
    valid_temp = temperature[np.isfinite(temperature)]

    if valid_temp.size == 0:
        return

    tmin = float(np.nanmin(valid_temp))
    tmax = float(np.nanmax(valid_temp))

    if np.isclose(tmin, tmax):
        ax.set_xlim(tmin - 1.0, tmax + 1.0)
    else:
        ax.set_xlim(tmin, tmax)


def set_exact_wv_limits(ax, water_vapor: np.ndarray):
    """Set upper x-axis exactly to min/max water vapor, excluding NaN/-9999."""
    valid_wv = water_vapor[np.isfinite(water_vapor)]

    if valid_wv.size == 0:
        return

    wmin = float(np.nanmin(valid_wv))
    wmax = float(np.nanmax(valid_wv))

    if np.isclose(wmin, wmax):
        pad = 0.1 if np.isclose(wmin, 0.0) else abs(wmin) * 0.05
        ax.set_xlim(wmin - pad, wmax + pad)
    else:
        ax.set_xlim(wmin, wmax)


def make_plot(
    ncfile: Path,
    pixel_index: Optional[int],
    swath_row: Optional[int],
    swath_col: Optional[int],
    output: Optional[Path],
    show: bool,
):
    with h5py.File(ncfile, "r") as f:
        i = choose_diagnostic_pixel(
            f=f,
            pixel_index=pixel_index,
            swath_row=swath_row,
            swath_col=swath_col,
        )

        pressure_hpa = clean_array(f["level"][:])

        temperature = array_at(f, "TemperatureProfile_Padded", i)
        water_vapor = array_at(f, "WaterVaporProfile_Padded", i)
        height_km = array_at(f, "HeightProfile", i)

        if temperature is None:
            raise ValueError("TemperatureProfile not found in file.")
        if water_vapor is None:
            raise ValueError("WaterVaporProfile not found in file.")
        if height_km is None:
            raise ValueError("HeightProfile not found in file.")

        # Metadata for title.
        row = scalar_at(f, "swath_row", i)
        col = scalar_at(f, "swath_col", i)
        flat_index = scalar_at(f, "flat_index_fortran0", i)
        original_row = scalar_at(f, "original_modis_row", i)
        original_col = scalar_at(f, "original_modis_col", i)
        lat = scalar_at(f, "Latitude", i)
        lon = scalar_at(f, "Longitude", i)
        mm_flag = scalar_at(f, "MM_Flag", i)
        mm_reason = scalar_at(f, "MM_FailureReason", i)

        markers = collect_required_markers(
            f=f,
            i=i,
            pressure_hpa=pressure_hpa,
            height_km_profile=height_km,
        )

    valid_y = (
        np.isfinite(pressure_hpa)
        & (pressure_hpa > 0)
        & (
            np.isfinite(temperature)
            | np.isfinite(water_vapor)
            | np.isfinite(height_km)
        )
    )
    pressure_hpa = pressure_hpa[valid_y]
    temperature = temperature[valid_y]
    water_vapor = water_vapor[valid_y]
    height_km = height_km[valid_y]

    fig, ax_temp = plt.subplots(figsize=(8.2, 8.8))

    # Lower x-axis: temperature.
    ax_temp.plot(
        temperature,
        pressure_hpa,
        color="red",
        marker="o",
        markersize=2.5,
        linewidth=1.5,
        label="Temperature",
        zorder=3,
    )

    ax_temp.set_xlabel("Temperature profile (K)", color="red")
    ax_temp.tick_params(axis="x", colors="red")
    ax_temp.spines["bottom"].set_color("red")

    set_exact_temp_limits(ax_temp, temperature)

    # Left y-axis: pressure.
    ax_temp.set_ylabel("Pressure (hPa)")
    ax_temp.set_yscale("log")

    # Limit top of y-axis to 10^1 hPa.
    p_bottom = np.nanmax(pressure_hpa)
    surface_pressures = [
        p for _, p, _, _, _ in markers
        if p is not None and np.isfinite(p)
    ]
    if surface_pressures:
        p_bottom = max(p_bottom, max(surface_pressures))

    ax_temp.set_ylim(p_bottom, 200.0)

    ax_temp.grid(True, which="both", linewidth=0.4, alpha=0.4)

    # Upper x-axis: water vapor.
    ax_wv = ax_temp.twiny()

    ax_wv.plot(
        water_vapor,
        pressure_hpa,
        color="blue",
        marker="s",
        markersize=2.3,
        linewidth=1.3,
        label="Water vapor",
        zorder=2,
    )

    ax_wv.set_xlabel("Water vapor profile (g kg$^{-1}$)", color="blue")
    ax_wv.tick_params(axis="x", colors="blue")
    ax_wv.spines["top"].set_color("blue")

    set_exact_wv_limits(ax_wv, water_vapor)

    # Right y-axis: height labels.
    ax_height = ax_temp.twinx()
    ax_height.set_yscale("log")
    ax_height.set_ylim(ax_temp.get_ylim())

    tick_pressures = np.array(
        [1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10],
        dtype=float,
    )

    p_min = 10.0
    p_max = p_bottom
    tick_pressures = tick_pressures[
        (tick_pressures >= p_min) & (tick_pressures <= p_max)
    ]

    tick_heights = [
        height_from_pressure(p, pressure_hpa, height_km)
        for p in tick_pressures
    ]

    tick_labels = [
        "" if h is None else f"{h:.1f}"
        for h in tick_heights
    ]

    ax_height.set_yticks(tick_pressures)
    ax_height.set_yticklabels(tick_labels)
    ax_height.set_ylabel("Height from HeightProfile (km)")

    # Requested horizontal lines only.
    x_positions = [0.04, 0.32, 0.44, 0.56, 0.68, 0.80]

    for n, (label, p, color, linestyle, linewidth) in enumerate(markers):
        add_horizontal_marker_inside(
            ax_temp,
            p,
            label,
            color=color,
            linestyle=linestyle,
            linewidth=linewidth,
            xfrac=x_positions[n % len(x_positions)],
        )

    # Title.
    title_lines = [
        f"{ncfile.name}",
        f"diagnostic_pixel index={i}",
    ]

    meta_parts = []

    if row is not None and col is not None:
        meta_parts.append(f"swath_row={int(row)}, swath_col={int(col)}")

    # if flat_index is not None:
    #     meta_parts.append(f"flat_index={int(flat_index)}")

    # if original_row is not None and original_col is not None:
    #     original_differs = (
    #         row is None
    #         or col is None
    #         or int(original_row) != int(row)
    #         or int(original_col) != int(col)
    #     )
    #     if original_differs:
    #         meta_parts.append(
    #             f"original_modis_row={int(original_row)}, "
    #             f"original_modis_col={int(original_col)}"
    #         )

    if lat is not None and lon is not None:
        meta_parts.append(f"lat={lat:.3f}, lon={lon:.3f}")

    if mm_flag is not None:
        meta_parts.append(f"MM_Flag={int(mm_flag)}")

    if mm_reason is not None:
        meta_parts.append(f"MM_FailureReason={int(mm_reason)}")

    if meta_parts:
        title_lines.append(", ".join(meta_parts))

    ax_temp.set_title("\n".join(title_lines), fontsize=10)

    fig.tight_layout()

    if output is None:
        output = Path(f"diag_profile_pixel{i}.png")

    fig.savefig(output, dpi=220, bbox_inches="tight")
    print(f"Saved plot to: {output}")

    if show:
        plt.show()

    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(
        description="Plot vertical diagnostic profile for one MM diagnostic pixel."
    )

    parser.add_argument(
        "ncfile",
        type=Path,
        help="Input MM_CO2FailureDiagnostics NetCDF/HDF5 file.",
    )

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument(
        "--pixel-index",
        type=int,
        default=None,
        help="Diagnostic-pixel index to plot, zero-based.",
    )

    group.add_argument(
        "--swath-row",
        type=int,
        default=None,
        help="Swath row. Must be used together with --swath-col.",
    )

    parser.add_argument(
        "--swath-col",
        type=int,
        default=None,
        help="Swath column. Must be used together with --swath-row.",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=None,
        help="Output PNG filename.",
    )

    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the figure interactively after saving.",
    )

    args = parser.parse_args()

    if args.swath_row is not None and args.swath_col is None:
        parser.error("--swath-row requires --swath-col.")

    if args.pixel_index is None and args.swath_row is None:
        parser.error("Provide either --pixel-index or --swath-row/--swath-col.")

    make_plot(
        ncfile=args.ncfile,
        pixel_index=args.pixel_index,
        swath_row=args.swath_row,
        swath_col=args.swath_col,
        output=args.output,
        show=args.show,
    )


if __name__ == "__main__":
    main()
