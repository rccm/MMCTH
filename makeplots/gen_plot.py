#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Multipanel MM plot with RGB in first panel + optional lat/lon lines (no projection).

Features
--------
- First panel: MODIS RGB from rad_1, rad_4, rad_3 with contrast enhancement.
- Remaining panels: user-chosen scalar variables with individual vmin/vmax.
- Optional latitude/longitude contour lines.
- Optional pixel-index overlay line drawn on all panels.
- Swath is auto-flipped so that north is up.
- Figure size roughly scales with screen size and number of columns.
"""

import numpy as np
import xarray as xr
import matplotlib
matplotlib.use("Agg")   # non-interactive backend
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import BoundaryNorm, ListedColormap
import matplotlib.ticker as mticker


# ---------- RGB helpers ----------

def bytescale(data, cmin=None, cmax=None, high=255, low=0):
    """Scale array to uint8 [low, high]."""
    if data.dtype == np.uint8:
        return data

    if high < low:
        raise ValueError("`high` should be larger than `low`.")

    if cmin is None:
        cmin = np.nanmin(data)
    if cmax is None:
        cmax = np.nanmax(data)

    cscale = cmax - cmin
    if cscale < 0:
        raise ValueError("`cmax` should be larger than `cmin`.")
    elif cscale == 0:
        cscale = 1

    scale = float(high - low) / cscale
    bytedata = (data - cmin) * scale + 0.4999
    bytedata = np.clip(bytedata, 0, high)
    return bytedata.astype(np.uint8) + np.uint8(low)


def get_enhanced_RGB(RGB):
    """Apply enhancement to a 3-band uint8 RGB (0–255)."""

    def scale_image(image):
        along_track, cross_track = image.shape

        x = np.array([0, 30, 60, 120, 190, 255], dtype=np.uint8)
        y = np.array([0, 110, 160, 210, 240, 255], dtype=np.uint8)

        scaled = np.zeros((along_track, cross_track), dtype=np.uint8)

        for i in range(len(x) - 1):
            x1, x2 = x[i], x[i + 1]
            y1, y2 = y[i], y[i + 1]
            m = (y2 - y1) / float(x2 - x1)
            b = y2 - m * x2

            mask = (image >= x1) & (image < x2)
            if np.any(mask):
                temp = np.asarray(m * image + b, dtype=np.uint8)
                scaled[mask] = temp[mask]

        mask = image >= x[-1]
        scaled[mask] = np.uint8(255)

        return scaled

    enhanced = np.zeros_like(RGB, dtype=np.uint8)
    for i in range(3):
        enhanced[:, :, i] = scale_image(bytescale(RGB[:, :, i]))

    return enhanced.astype(np.float32) / 255.0


def _make_rgb(r1, r4, r3):
    """Build enhanced RGB from three bands."""
    rgbdata = np.dstack((r1, r4, r3))
    return get_enhanced_RGB(rgbdata)


# ---------- Main plotting function ----------

def plot_mm_multipanel(
    ds,
    varnames,
    nrows=2,
    ncols=3,
    cmap="turbo",
    cmaps=None,
    norms=None,
    discrete=None,
    add_latlon_lines=False,
    show_lat_labels=False,
    show_lon_labels=False,
    lat_step=5.0,
    lon_step=10.0,
    vlims=None,
    titles=None,
    rgb_title="MODIS RGB",
    output="output.png",
    draw_latlon_label=False,
    RGB_AX0=True,
    flip_cbar_vars=None,
    overlay_xy=None,          # NEW: (x_idx, y_idx) line in pixel coordinates
    overlay_kwargs=None,      # NEW: dict of ax.plot kwargs
):
    cmaps = cmaps or {}
    norms = norms or {}
    discrete = discrete or {}
    flip_cbar_vars = set(flip_cbar_vars or [])

    # --- latitude / longitude and mask ---
    lat = ds["Latitude"].values
    lon = ds["Longitude"].values
    lat_fv = ds["Latitude"].attrs.get("_FillValue", None)
    lon_fv = ds["Longitude"].attrs.get("_FillValue", None)

    mask = ~np.isfinite(lat) | ~np.isfinite(lon)
    if lat_fv is not None:
        mask |= (lat == lat_fv)
    if lon_fv is not None:
        mask |= (lon == lon_fv)

    lat_m = np.ma.masked_where(mask, lat)
    lon_m = np.ma.masked_where(mask, lon)

    # --- ensure north is up: flip vertically if needed ---
    ny, nx = lat_m.shape
    row_mean_lat = np.array([np.nanmean(lat_m[i, :]) for i in range(ny)])
    flip_vert = row_mean_lat[0] > row_mean_lat[-1]
    if flip_vert:
        lat_m = np.flipud(lat_m)
        lon_m = np.flipud(lon_m)
        mask = np.flipud(mask)

    # --- optional overlay line (pixel indices) ---
    overlay_x = overlay_y = None
    if overlay_xy is not None:
        overlay_x = np.asarray(overlay_xy[0], dtype=float)
        overlay_y = np.asarray(overlay_xy[1], dtype=float)

        if overlay_x.shape != overlay_y.shape:
            raise ValueError("overlay_xy must be (x_array, y_array) with the same length")

        # If image is flipped vertically, transform y indices to match displayed image
        if flip_vert:
            overlay_y = (ny - 1) - overlay_y

    def draw_overlay(ax):
        if overlay_x is None or overlay_y is None:
            return
        kw = dict(color="magenta", linewidth=1.8, linestyle="-", marker="o", markersize=2)
        if overlay_kwargs:
            kw.update(overlay_kwargs)
        ax.plot(overlay_x, overlay_y, **kw)

    # --- build RGB, apply same flip & mask ---
    rgb = _make_rgb(ds["rad_1"].values, ds["rad_4"].values, ds["rad_3"].values)
    if flip_vert:
        rgb = np.flipud(rgb)
    rgb[mask] = np.nan

    # --- label formatters ---
    def _fmt_lat(v): return f"{abs(v):.0f}°" + ("N" if v >= 0 else "S")
    def _fmt_lon(v): return f"{abs(v):.0f}°" + ("E" if v >= 0 else "W")

    # --- lat/lon contour levels ---
    if add_latlon_lines or show_lat_labels or show_lon_labels:
        lat_min, lat_max = np.nanmin(lat_m), np.nanmax(lat_m)
        lon_min, lon_max = np.nanmin(lon_m), np.nanmax(lon_m)
        lat_levels = np.arange(np.floor(lat_min / lat_step) * lat_step,
                               np.ceil(lat_max / lat_step) * lat_step + 0.1,
                               lat_step)
        lon_levels = np.arange(np.floor(lon_min / lon_step) * lon_step,
                               np.ceil(lon_max / lon_step) * lon_step + 0.1,
                               lon_step)
    else:
        lat_levels = lon_levels = None

    def add_latlon_contours(ax):
        if not add_latlon_lines:
            return
        if lat_levels is not None:
            ax.contour(lat_m, levels=lat_levels, colors="k", linewidths=0.4, alpha=0.4, origin="lower")
        if lon_levels is not None:
            ax.contour(lon_m, levels=lon_levels, colors="k", linewidths=0.4, alpha=0.4, origin="lower")

    # --- tick positions (where contours cross frame boundary) ---
    def compute_axis_ticks():
        ny_local, nx_local = lat_m.shape
        yticks, ylabels, xticks, xlabels = [], [], [], []

        j_edge = None
        if lat_levels is not None and show_lat_labels:
            for j in range(nx_local):
                if np.isfinite(lat_m[:, j]).sum() > 5:
                    j_edge = j
                    break

        i_edge = None
        if lon_levels is not None and show_lon_labels:
            for i in range(ny_local):
                if np.isfinite(lon_m[i, :]).sum() > 5:
                    i_edge = i
                    break

        def _edge_crossings(a, lev):
            idxs = []
            diff = a - lev
            for k in range(len(a) - 1):
                if not (np.isfinite(diff[k]) and np.isfinite(diff[k + 1])):
                    continue
                if diff[k] == 0:
                    idxs.append(float(k))
                elif diff[k + 1] == 0:
                    idxs.append(float(k + 1))
                elif diff[k] * diff[k + 1] < 0:
                    frac = -diff[k] / (diff[k + 1] - diff[k])
                    idxs.append(k + frac)
            return idxs

        if j_edge is not None:
            col_lat = lat_m[:, j_edge]
            for lev in lat_levels:
                crosses = _edge_crossings(col_lat, lev)
                if crosses:
                    y = float(np.mean(crosses))
                    if not any(abs(y - y0) < 0.5 for y0 in yticks):
                        yticks.append(y)
                        ylabels.append(_fmt_lat(lev))

        if i_edge is not None:
            row_lon = lon_m[i_edge, :]
            for lev in lon_levels:
                crosses = _edge_crossings(row_lon, lev)
                if crosses:
                    x = float(np.mean(crosses))
                    if not any(abs(x - x0) < 0.5 for x0 in xticks):
                        xticks.append(x)
                        xlabels.append(_fmt_lon(lev))

        return yticks, ylabels, xticks, xlabels

    yticks, ylabels, xticks, xlabels = compute_axis_ticks()

    def compute_right_lon_ticks():
        if lon_levels is None or not show_lon_labels:
            return [], []

        ny_local, nx_local = lon_m.shape

        j_right = None
        for j in range(nx_local - 1, -1, -1):
            if np.isfinite(lon_m[:, j]).sum() > 5:
                j_right = j
                break
        if j_right is None:
            return [], []

        col_lon = lon_m[:, j_right]

        def _edge_crossings(a, lev):
            idxs = []
            diff = a - lev
            for k in range(len(a) - 1):
                if not (np.isfinite(diff[k]) and np.isfinite(diff[k + 1])):
                    continue
                if diff[k] == 0:
                    idxs.append(float(k))
                elif diff[k + 1] == 0:
                    idxs.append(float(k + 1))
                elif diff[k] * diff[k + 1] < 0:
                    frac = -diff[k] / (diff[k + 1] - diff[k])
                    idxs.append(k + frac)
            return idxs

        right_ticks, right_labels = [], []
        for lev in lon_levels:
            crosses = _edge_crossings(col_lon, lev)
            if not crosses:
                continue
            y = float(np.mean(crosses))
            if any(abs(y - y0) < 0.5 for y0 in right_ticks):
                continue
            right_ticks.append(y)
            right_labels.append(_fmt_lon(lev))

        return right_ticks, right_labels

    right_lon_yticks, right_lon_ylabels = compute_right_lon_ticks()

    # --- figure size ---
    dpi = plt.rcParams["figure.dpi"]
    screen_w = 2560
    h_data, w_data = lat.shape
    aspect = h_data / w_data
    base_height = screen_w / dpi * 0.4
    fig_height = base_height * max(nrows, 1)
    fig_width = base_height / aspect * (ncols + 1.5)

    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(fig_width, fig_height),
        gridspec_kw={"left": 0.02, "right": 0.98, "wspace": 0.02, "hspace": 0.25},
    )
    axes = np.atleast_1d(axes).ravel()

    def add_cbar(ax, mappable, *, disc_cfg=None, label=None, labelsize=14, ticksize=14,
                 sci=False, sci_ndp=0, flip=False):
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.15, axes_class=plt.Axes)

        if disc_cfg is None:
            cb = fig.colorbar(mappable, cax=cax)
        else:
            cb = fig.colorbar(
                mappable, cax=cax,
                ticks=disc_cfg.get("ticks", None),
                boundaries=disc_cfg.get("boundaries", None),
                spacing="proportional",
            )
            ticklabels = disc_cfg.get("ticklabels", None)
            if ticklabels is not None and disc_cfg.get("ticks", None) is not None:
                cb.set_ticklabels(ticklabels)

        if sci and (disc_cfg is None):
            cb.update_ticks()
            if cb.orientation == "horizontal":
                cb.ax.xaxis.get_offset_text().set_visible(False)
            else:
                cb.ax.yaxis.get_offset_text().set_visible(False)

        cb.ax.tick_params(labelsize=ticksize)
        if label is not None:
            cb.set_label(label, fontsize=labelsize)

        if flip:
            if cb.orientation == "horizontal":
                cb.ax.invert_xaxis()
            else:
                cb.ax.invert_yaxis()

        return cb

    # --- panel 0: RGB ---
    ax_start_idx = 0
    if RGB_AX0:
        ax_start_idx = 1
        ax0 = axes[0]
        ax0.imshow(rgb, origin="lower")
        add_latlon_contours(ax0)
        draw_overlay(ax0)  # NEW

        ax0.set_yticks(yticks if (show_lat_labels and yticks) else [])
        if show_lat_labels and yticks:
            ax0.set_yticklabels(ylabels, fontsize=10)

        ax0.set_xticks(xticks if (show_lon_labels and xticks) else [])
        if show_lon_labels and xticks:
            ax0.set_xticklabels(xlabels, fontsize=10)

        if right_lon_yticks:
            x_text = rgb.shape[1] - 1 + 20
            for y, lab in zip(right_lon_yticks, right_lon_ylabels):
                ax0.text(
                    x_text, y, lab,
                    va="center", ha="left",
                    rotation=90,
                    rotation_mode="anchor",
                    fontsize=10,
                    clip_on=False
                )
        ax0.set_title(rgb_title, fontsize=12, pad=15)

    # --- scalar panels ---
    vlims = vlims or {}
    titles = titles or {}

    def build_discrete_for_var(name, da_values):
        cfg = discrete.get(name, None)
        if cfg is None:
            return None, None, None

        if isinstance(cfg, dict) and "categories" in cfg:
            cats = np.array(cfg["categories"], dtype=float)
            boundaries = np.concatenate(([cats.min() - 0.5], cats + 0.5))
            boundaries = np.unique(boundaries)

            if "colors" in cfg and cfg["colors"] is not None:
                cmap_obj = ListedColormap(cfg["colors"])
            else:
                cmap_obj = plt.get_cmap(cfg.get("cmap", "tab10"), len(cats))

            norm_obj = BoundaryNorm(boundaries, cmap_obj.N, clip=True)

            ticks = cats
            labels = cfg.get("labels", [str(int(c)) if float(c).is_integer() else str(c) for c in cats])
            cbar_cfg = {"boundaries": boundaries, "ticks": ticks, "ticklabels": labels}
            return cmap_obj, norm_obj, cbar_cfg

        boundaries = np.asarray(cfg["boundaries"], dtype=float)
        if boundaries.ndim != 1 or boundaries.size < 2:
            raise ValueError(f"discrete['{name}']['boundaries'] must be a 1D array with >=2 values")

        cmap_obj = cfg.get("cmap", "tab20")
        if isinstance(cmap_obj, str):
            cmap_obj = plt.get_cmap(cmap_obj, boundaries.size - 1)

        norm_obj = BoundaryNorm(boundaries, cmap_obj.N, clip=False)

        if "ticks" in cfg:
            ticks = cfg["ticks"]
        else:
            ticks = 0.5 * (boundaries[:-1] + boundaries[1:])

        cbar_cfg = {"boundaries": boundaries, "ticks": ticks}
        if "ticklabels" in cfg:
            cbar_cfg["ticklabels"] = cfg["ticklabels"]

        return cmap_obj, norm_obj, cbar_cfg

    for ax, name in zip(axes[ax_start_idx:], varnames):
        if name not in ds:
            print(f"Variable not found in file: {name}")
            ax.set_visible(False)
            continue

        da = ds[name].where(np.isfinite(ds[name])).values
        if flip_vert:
            da = np.flipud(da)
        da = np.where(mask, np.nan, da)

        cmap_i = cmaps.get(name, cmap)
        flip_this_bar = (name in flip_cbar_vars)

        disc_cmap, disc_norm, disc_cbar = build_discrete_for_var(name, da)
        if disc_norm is not None:
            cmap_i = disc_cmap
            norm_i = disc_norm
            vmin = vmax = None
        else:
            norm_i = norms.get(name, None)
            vmin, vmax = vlims.get(name, (None, None))

        im = ax.imshow(da, origin="lower", cmap=cmap_i, vmin=vmin, vmax=vmax, norm=norm_i)
        add_latlon_contours(ax)
        draw_overlay(ax)  # NEW

        if draw_latlon_label:
            ax.set_yticks(yticks if (show_lat_labels and yticks) else [])
            if show_lat_labels and yticks:
                ax.set_yticklabels(ylabels, fontsize=10)

            ax.set_xticks(xticks if (show_lon_labels and xticks) else [])
            if show_lon_labels and xticks:
                ax.set_xticklabels(xlabels, fontsize=10)
        else:
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlabel("")
            ax.set_ylabel("")

        ax.set_title(titles.get(name, name), fontsize=14, pad=15)
        add_cbar(ax, im, disc_cfg=disc_cbar, labelsize=14, ticksize=14, sci=False, sci_ndp=0, flip=flip_this_bar)

    # Hide unused axes correctly whether RGB is on or off
    used_axes = ax_start_idx + len(varnames)
    for ax in axes[used_axes:]:
        ax.set_visible(False)

    fig.savefig(output, dpi=300)
    return fig