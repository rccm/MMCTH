#!/usr/bin/env python3
import numpy as np
import xarray as xr

FILE_A = "/data/gdi/f/gzhao1/mmcth/output/MODMISR_L2_CP_061_A2015176.1650_0625_O082545_MM_V00.nc"
FILE_B = "/data/gdi/f/gzhao1/mmcth/ds_output/debug/MODMISR_L2_CP_061_A2015176.1650_0625_O082545_MM_V00.nc"

FLAG_VAR = "MM_Flag"
HGT_VAR  = "MM_CloudTopHeight"

def mean_mmcth_flag_gt5(nc_path):
    ds = xr.open_dataset(nc_path)
    flag = ds[FLAG_VAR].values
    hgt  = ds[HGT_VAR].values

    mask = (flag > 5) & np.isfinite(hgt) & (hgt > 0)
    if mask.sum() == 0:
        return np.nan, 0

    return float(hgt[mask].mean()), int(mask.sum())

def main():
    meanA, nA = mean_mmcth_flag_gt5(FILE_A)
    meanB, nB = mean_mmcth_flag_gt5(FILE_B)

    print(f"A: mean({HGT_VAR}) for {FLAG_VAR} > 5 (finite & >0), n={nA}: {meanA}")
    print(f"B: mean({HGT_VAR}) for {FLAG_VAR} > 5 (finite & >0), n={nB}: {meanB}")

if __name__ == "__main__":
    main()