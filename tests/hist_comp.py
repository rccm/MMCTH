#!/usr/bin/env python3
import numpy as np
import xarray as xr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FILE_A = "/data/gdi/f/gzhao1/mmcth/output/MODMISR_L2_CP_061_A2015176.1650_0625_O082545_MM_V00.nc"
FILE_B = "/data/gdi/f/gzhao1/mmcth/ds_output/debug/MODMISR_L2_CP_061_A2015176.1650_0625_O082545_MM_V00.nc"

FLAG_VAR = "MM_Flag"
HGT_VAR  = "MM_CloudTopHeight"

BINS = 120
OUT_PNG = f"{HGT_VAR}_compare_steps_maskFromA_{FLAG_VAR}_gt5_lt8.png"

def main():
    dsA = xr.open_dataset(FILE_A)
    dsB = xr.open_dataset(FILE_B)

    flagA = dsA[FLAG_VAR].values
    hgtA  = dsA[HGT_VAR].values
    hgtB  = dsB[HGT_VAR].values

    # Must align 1:1 to apply FILE_A-based mask to both
    if hgtA.shape != hgtB.shape or flagA.shape != hgtA.shape:
        raise RuntimeError(
            f"Shape mismatch:\n"
            f"  flagA {flagA.shape}\n"
            f"  hgtA  {hgtA.shape}\n"
            f"  hgtB  {hgtB.shape}\n"
            "These must match to apply FILE_A-based mask to both files."
        )

    # Mask defined ONLY from FILE_A MM_Flag: 5 < flag < 8 (i.e., 6 or 7)
    # Also require valid heights in BOTH files so you're comparing like-for-like points
    mask = (
        (flagA > 5) & (flagA < 8) &
        np.isfinite(hgtA) & (hgtA > 0) &
        np.isfinite(hgtB) & (hgtB > 0)
    )

    xA = np.asarray(hgtA[mask]).ravel()
    xB = np.asarray(hgtB[mask]).ravel()

    print(f"Total points: {flagA.size}")
    print(f"Kept (mask from A, valid in both): {mask.sum()}")

    if xA.size == 0:
        raise RuntimeError("No points left after masking. Check MM_Flag values and height validity.")

    # Common bins from combined range
    mn = float(min(xA.min(), xB.min()))
    mx = float(max(xA.max(), xB.max()))
    bin_edges = np.linspace(mn, mx, BINS + 1)

    plt.figure(figsize=(10, 6))
    plt.hist(xA, bins=bin_edges, histtype="step", linewidth=1.8,
             label=f"ORG(n={xA.size})")
    plt.hist(xB, bins=bin_edges, histtype="step", linewidth=1.8,
             label=f"DS (n={xB.size})")

    plt.title(f"{HGT_VAR} comparison \nMask from ORG 5 < {FLAG_VAR} < 8, finite & >0 in BOTH")
    plt.xlabel(f"{HGT_VAR} (m)")
    plt.ylabel("Count")
    plt.legend()
    plt.tight_layout()

    plt.savefig(OUT_PNG, dpi=200)
    print(f"Wrote: {OUT_PNG}")

if __name__ == "__main__":
    main()