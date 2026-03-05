import logging
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg") 
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../src')))
from src.data_readers  import mod_read
in_file = 't1.17261.1500.1000m.hdf'
out_file = 't1.17261.1500.1000m.destripe.hdf'
#in_file = 'test.17261.1500.1000m.hdf'
# in_file ='/data/gdi/satellite/TerraDataArchive/MODIS/MOD021KM/2015/176/MOD021KM.A2015176.1650.061.2017321013930.hdf'


testin_org = mod_read.MODL1Granule(in_file, destripe_flag=False)
testin_org_rad = {f'rad_{band}': testin_org.get_band_data(str(band), 'radiance') for band in [36, 35, 34, 33, 31]}


testin_ds = mod_read.MODL1Granule(in_file, destripe_flag=True)
testin_ds_rad = {f'rad_{band}': testin_ds.get_band_data(str(band), 'radiance') for band in [36, 35, 34, 33, 31]}
c
testout = mod_read.MODL1Granule(out_file,destripe_flag=False)
testout_ds = {f'rad_{band}': testout.get_band_data(str(band), 'radiance') for band in [36, 35, 34, 33, 31]}

band = '36'
testin_org_rad = testin_org_rad[f"rad_{band}"]
testin_ds_rad = testin_ds_rad[f"rad_{band}"]
testout_ds_rad = testout_ds[f"rad_{band}"]



def save_gray(img, out_png, title, mask=None, p_lo=1, p_hi=99):
    """Save a grayscale image with percentile-based contrast limits."""
    arr = np.array(img, dtype=np.float64)

    if mask is None:
        mask = np.isfinite(arr)
    else:
        mask = mask & np.isfinite(arr)

    if not np.any(mask):
        raise ValueError(f"No valid pixels to plot for {out_png}")

    vmin, vmax = np.percentile(arr[mask], [p_lo, p_hi])

    plt.figure(figsize=(6, 5))
    plt.imshow(arr, cmap="gray", vmin=vmin, vmax=vmax, origin="upper")
    plt.colorbar(label="Radiance")
    plt.title(title)
    plt.axis("off")
    plt.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close()

save_gray(testin_org_rad,f"testin_band{band}_org_gray.png", "Band 36 Radiance (ORG) TestIn")
save_gray(testin_ds_rad, f"testin_band{band}_ds_gray.png", "Band 36 Radiance (DS) TestIn")
save_gray(testout_ds_rad, f"testout_band{band}_ds_gray.png", "Band 36 Radiance (DS) TestOut")

diff = testin_ds_rad - testout_ds_rad
save_gray(diff, f"band{band}_diff_gray.png", f"Band {band} Difference ( Testin_DS  - TestOut_DS)")





# diff = arka - dong

# mask = (arka_ds['rad_36'] >= 0) & np.isfinite(diff)


# eps = 1e-12  # avoid divide-by-zero; adjust if radiance can be extremely small
# pct = 100.0 * (arka - dong) / (dong + eps)

# mask = (arka >= 0) & np.isfinite(pct) & (np.abs(dong) > eps)

# # extrema (in percent)
# pct_for_max = np.where(mask, pct, -np.inf)
# pct_for_min = np.where(mask, pct,  np.inf)

# imax = np.unravel_index(np.argmax(pct_for_max), pct.shape)
# imin = np.unravel_index(np.argmin(pct_for_min), pct.shape)

# print("MAX % diff (relative to dong):")
# print("  idx:", imax)
# print("  arka:", arka[imax])
# print("  dong:", dong[imax])
# print("  pct:", pct[imax])

# print("\nMIN % diff (relative to dong):")
# print("  idx:", imin)
# print("  arka:", arka[imin])
# print("  dong:", dong[imin])
# print("  pct:", pct[imin])

# pct_1d = pct[mask].ravel()

# binwidth = 0.1
# mn, mx = pct_1d.min(), pct_1d.max()
# bins = np.arange(mn, mx + binwidth, binwidth)

# weights = np.ones(pct_1d.shape, dtype=np.float64) / pct_1d.size

# # plt.figure()
# # counts, edges, patches = plt.hist(pct_1d, bins=bins, weights=weights, log=True)
# # plt.xlabel("% difference = 100*(Arka - Dongwei)/Dongwei")
# # plt.ylabel("Fraction of pixels (log scale)")
# # plt.title(f"Fraction-per-bin histogram of % diff (binwidth={binwidth})")
# # plt.savefig("test_his_pct_frac_logy.png", dpi=200, bbox_inches="tight")

# # print("Sum of bin fractions (should be ~1):", counts.sum())


# # plt.figure()
# # counts, edges, patches = plt.hist(pct_1d, bins=bins, weights=weights, log=True)
# # plt.xlabel("% difference = 100*(Arka - Dongwei)/Dongwei")
# # plt.ylabel("Fraction of pixels (log scale)")
# # plt.title(f"Fraction-per-bin histogram of % diff (binwidth={binwidth})")
# # plt.savefig("test_his_pct_frac_logy.png", dpi=200, bbox_inches="tight")

# # print("Sum of bin fractions (should be ~1):", counts.sum())

# mask_diff = (arka >= 0) & np.isfinite(diff)
# diff_1d = diff[mask_diff].ravel()

# binwidth = 0.001   # choose units appropriate for radiance
# mn, mx = diff_1d.min(), diff_1d.max()
# bins = np.arange(mn, mx + binwidth, binwidth)



# weights = np.ones(diff_1d.shape, dtype=np.float64) / diff_1d.size

# plt.figure()
# counts, edges, patches = plt.hist(diff_1d, bins=bins, weights=weights, log=True)
# plt.xlabel("Raw difference = Arka - Dongwei (radiance units)")
# plt.ylabel("Fraction of pixels (log scale)")
# plt.title(f"Fraction-per-bin histogram of raw diff (binwidth={binwidth})")
# plt.savefig("test_his_rawdiff_frac_logy.png", dpi=200, bbox_inches="tight")

# print("Sum of bin fractions (should be ~1):", counts.sum())


