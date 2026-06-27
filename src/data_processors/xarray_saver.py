#!/usr/bin/env python3
"""
xarray_saver.py — CF-oriented NetCDF writer for MODIS/MISR/ERA5 fused scenes.

Key CF choices
--------------
- Conventions: CF-1.10
- Swath variables remain on index dimensions:
    ('cell_along_swath_1km', 'cell_cross_swath_1km')
  with 2-D auxiliary coordinates Latitude and Longitude and a 'coordinates' attribute.
- ERA5 variables use an explicit pressure vertical coordinate variable z(z).
- Optional scalar scene time coordinate can be included.
- Optional row_time and pixel_time auxiliary coordinates can be included.
- Integer flags use _FillValue outside valid flag_values.
- Categorical flag/code variables do not carry a units attribute.
"""

from __future__ import annotations

import logging
import re
from datetime import datetime, timezone, timedelta
from typing import Optional, Dict, Any, List

import numpy as np
import xarray as xr


class XarraySaver:
    """
    Simple NetCDF writer with a single, easy-to-edit schema for names and attributes.

    Dimensions
    ----------
    MM / MODIS / MISR 2D:
        ('cell_along_swath_1km', 'cell_cross_swath_1km')

    ERA5 2D:
        ('era5_cell_along_swath_10km', 'era5_cell_cross_swath_10km')

    ERA5 3D:
        ('z', 'era5_cell_along_swath_10km', 'era5_cell_cross_swath_10km')
        where z is pressure in hPa
    """

    MM_SCHEMA: Dict[str, Dict[str, Any]] = {
        "cth": {
            "name": "MM_CloudTopHeight",
            "units": "m",
            "dtype": "int16",
            "valid_range": np.array([0, 20000], dtype=np.int16),
            "long_name": "MODIS/MISR fusion cloud top height",
        },
        "ctp": {
            "name": "MM_CloudTopPressure",
            "units": "hPa",
            "dtype": "float32",
            "long_name": "MODIS/MISR fusion cloud top pressure",
        },
        "emissivity": {
            "name": "MM_CloudEffectiveEmissivity",
            "units": "1",
            "dtype": "float32",
            "long_name": "MODIS/MISR fusion cloud effective emissivity",
        },
        "opt": {
            "name": "MM_CloudOpticalDepth",
            "units": "1",
            "dtype": "float32",
            "long_name": "MODIS/MISR fusion cloud optical depth",
        },
        "cflag": {
            "name": "MM_Flag",
            "dtype": "int8",
            "flag_values": np.array([0, 1, 2, 3, 4, 5, 6, 7,8,9], dtype=np.int8),
            "flag_meanings": (
                "no_retrieval"
                "MODIS_only"
                "MISR_only"
                "MISR_selected; MODIS CTH < MISR CTH"
                "MISR_selected; likely one-layer; BT applied"
                "MODIS_selected; likely multi-layer; BT applied"
                "MISR_selected; MM attempted; likely one-layer; CO2-slicing applied"
                "MODIS_selected; MM attempted likely multi-layer; CO2-slicing applied"
                "MM_success_band_pair_36_35"
                "MM_success_band_pair_35_33"
            ),
            "long_name": "cloud classification and quality flag",
        },
        "ml_zscore": {
            "name": "MM_MultilayerZScore",
            "units": "1",
            "dtype": "float32",
            "long_name": (
                "Z score for MM minus MISR cloud top height difference "
                "under the same-layer hypothesis"
            ),
            "comment": (
                "Computed as [(MM_CloudTopHeight - MISR_CloudTopHeight) "
                "- expected_same_layer_difference] divided by combined "
                "same-layer uncertainty. Positive values indicate increasing "
                "evidence that MM and MISR are retrieving different cloud layers."
            ),
        },
        "ml_flag": {
            "name": "MM_MultilayerConfidenceFlag",
            "dtype": "int8",
            "flag_values": np.array([0, 1, 2, 3, 4], dtype=np.int8),
            "flag_meanings": (
                "not_evaluated"
                "weak_or_no_multilayer_evidence"
                "moderate_multilayer_evidence"
                "strong_multilayer_evidence"
                "very_strong_multilayer_evidence"
            ),
            "long_name": (
                "multilayer confidence flag based on MM and MISR cloud top height difference"
            ),
            "comment": (
                "One-sided confidence under the same-layer null hypothesis. "
                "0 = not evaluated; "
                "1 = confidence less than 90 percent; "
                "2 = confidence from 90 to less than 95 percent; "
                "3 = confidence from 95 to less than 99 percent; "
                "4 = confidence greater than or equal to 99 percent."
            ),
        },
        "failure_reason": {
            "name": "MM_FailureReason",
            "dtype": "int8",
            "flag_values": np.array(
                [0, 10, 20, 21, 22, 23, 24, 30, 40, 41, 42, 50, 51, 52, 60, 61, 62, 70, 71],
                dtype=np.int8,
            ),
            "flag_meanings": (
                "ok "
                "not_attempted "
                "invalid_misr_ctp "
                "invalid_surface_pressure "
                "invalid_modis_radiance "
                "invalid_atmospheric_profile "
                "pressure_grid_mapping_failed "
                "no_vertical_search_domain "
                "signal_below_noise "
                "no_valid_ratio "
                "no_ratio_crossing "
                "candidate_not_above_misr "
                "candidate_outside_band_range "
                "no_selectable_band_pair "
                "candidate_not_higher_than_one_layer "
                "reference_one_layer_failed "
                "candidate_has_no_valid_ctp_digit "
                "ctp_ok_emissivity_denominator_bad "
                "ctp_ok_emissivity_ratio_bad"
            ),
            "long_name": "MM two-layer retrieval failure reason",
            "comment": (
                "Diagnostic reason for the MM two-layer CO2-slicing path. "
                "Use MM_Flag to identify the final cloud-top source. Codes 70 "
                "and 71 mean CTP was accepted but emissivity failed."
            ),
        },
        # "mm_attemp": {
        #     "name": "MM_RetrievalStatus",
        #     "dtype": "int8",
        #     "flag_values": np.array([0, 1, 2, 3], dtype=np.int8),
        #     "flag_meanings": (
        #         "not_attempted "
        #         "successful_and_selected "
        #         "attempted_but_not_selected "
        #         "attempted_but_no_valid_solution"
        #     ),
        #     "long_name": "MM two-layer retrieval attempt status",
        #     "comment": (
        #         "Diagnostic status of the MM two-layer retrieval attempt. "
        #         "This field describes whether the MM retrieval was attempted and whether it was selected. "
        #         "The final cloud-top source is given by MM_Flag."
        #     ),
        # }
    }

    MODIS_SCHEMA: Dict[str, Dict[str, Any]] = {
        "cth": {
            "name": "MODIS_CloudTopHeight",
            "units": "m",
            "dtype": "int16",
            "long_name": "MODIS cloud top height",
        },
        "ctp": {
            "name": "MODIS_CloudTopPressure",
            "units": "hPa",
            "dtype": "float32",
            "long_name": "MODIS cloud top pressure",
        },
        "emissivity": {
            "name": "MODIS_CloudEffectiveEmissivity",
            "units": "1",
            "dtype": "float32",
            "long_name": "MODIS cloud effective emissivity",
        },
        "opt": {
            "name": "MODIS_CloudOpticalDepth",
            "units": "1",
            "dtype": "float32",
            "long_name": "MODIS cloud optical depth",
        },
        "ctm": {
            "name": "MODIS_CloudTopMethod",
            "dtype": "int8",
            "long_name": "MODIS cloud top retrieval method",
        },
        "cmulti": {
            "name": "MODIS_CloudMultiLayerFlag",
            "dtype": "int8",
            "long_name": "MODIS Cloud Multi Layer Identification From Shortwave Observations",
        },
        "cqa": {
            "name": "MODIS_CloudMultilayerFlagFromQA",
            "dtype": "int8",
            "long_name": "MODIS Cloud Multi Layer Identification From Quality Assurance Decode",
        },
        "cphase": {
            "name": "MODIS_CloudPhase",
            "dtype": "int8",
            "flag_values": np.array([0, 1, 2, 6], dtype=np.int8),
            "flag_meanings": "clear liquid_water_clouds ice_clouds undetermined_phase_clouds",
            "long_name": "MODIS cloud phase",
        },
    }

    MISR_SCHEMA: Dict[str, Dict[str, Any]] = {
        "misrcth": {
            "name": "MISR_CloudTopHeight",
            "units": "m",
            "dtype": "int16",
            "long_name": "MISR cloud top height without wind correction",
        },
        "misrctp": {
            "name": "MISR_CloudTopPressure",
            "units": "hPa",
            "dtype": "float32",
            "long_name": "MISR cloud top pressure converted from cloud top height",
        },
        "misrcth_qa": {
            "name": "MISR_CloudTopHeight_QI",
            "dtype": "int16",
            "valid_range": np.array([0, 128], dtype=np.int16),
            "long_name": "MISR cloud top height quality indicator",
            "comment": "0 to 100 are quality scores; 128 indicates bad data",
        },
    }

    ERA5_SCHEMA: Dict[str, Dict[str, Any]] = {
        # 3D variables
        "geopotential_org": {
            "name": "ERA5_Geopotential",
            "units": "m2 s-2",
            "dtype": "float32",
            "ndim": 3,
            "long_name": "ERA5 geopotential",
            "standard_name": "geopotential",
        },
        "temperature_org": {
            "name": "ERA5_Temperature",
            "units": "K",
            "dtype": "float32",
            "ndim": 3,
            "long_name": "ERA5 air temperature",
            "standard_name": "air_temperature",
        },
        "mixingratio_org": {
            "name": "ERA5_WaterVaporMixingRatio",
            "units": "1",
            "dtype": "float32",
            "ndim": 3,
            "long_name": "ERA5 water vapor mixing ratio",
        },
        "u": {
            "name": "ERA5_U",
            "units": "m s-1",
            "dtype": "float32",
            "ndim": 3,
            "long_name": "ERA5 eastward wind",
            "standard_name": "eastward_wind",
        },
        "v": {
            "name": "ERA5_V",
            "units": "m s-1",
            "dtype": "float32",
            "ndim": 3,
            "long_name": "ERA5 northward wind",
            "standard_name": "northward_wind",
        },
        "w": {
            "name": "ERA5_W",
            "units": "Pa s-1",
            "ndim": 3,
            "dtype": "float32",
            "long_name": "ERA5 pressure vertical velocity",
            "standard_name": "lagrangian_tendency_of_air_pressure",
        },

        # 2D variables
        "ERA5Humidity": {
            "name": "ERA5_SurfaceWaterVaporMixingRatio",
            "units": "1",
            "dtype": "float32",
            "ndim": 2,
            "long_name": "ERA5 water vapor mixing ratio at surface or 2 m",
        },
        "dew2m": {
            "name": "ERA5_DewPoint2m",
            "units": "K",
            "dtype": "float32",
            "ndim": 2,
            "long_name": "ERA5 2 m dew point temperature",
            "standard_name": "dew_point_temperature",
        },
        "sst": {
            "name": "ERA5_SeaSurfaceTemperature",
            "units": "K",
            "dtype": "float32",
            "ndim": 2,
            "long_name": "ERA5 sea surface temperature",
            "standard_name": "sea_surface_temperature",
        },
        "temp2m": {
            "name": "ERA5_Temperature2m",
            "units": "K",
            "dtype": "float32",
            "ndim": 2,
            "long_name": "ERA5 2 m air temperature",
            "standard_name": "air_temperature",
        },
        "skint": {
            "name": "ERA5_SkinTemperature",
            "units": "K",
            "dtype": "float32",
            "ndim": 2,
            "long_name": "ERA5 skin temperature",
        },
        "msp": {
            "name": "ERA5_MeanSeaLevelPressure",
            "units": "hPa",
            "dtype": "float32",
            "ndim": 2,
            "long_name": "ERA5 mean sea level pressure",
        },
        "swmr": {
            "name": "ERA5_2mWaterVaporMixingRatio",
            "units": "1",
            "dtype": "float32",
            "ndim": 2,
            "long_name": "ERA5 2 m water vapor mixing ratio",
        },
        "ERA5Latitude": {
            "name": "ERA5_Latitude",
            "units": "degrees_north",
            "dtype": "float32",
            "ndim": 2,
            "standard_name": "latitude",
            "long_name": "ERA5 grid latitude",
        },
        "ERA5Longitude": {
            "name": "ERA5_Longitude",
            "units": "degrees_east",
            "dtype": "float32",
            "ndim": 2,
            "standard_name": "longitude",
            "long_name": "ERA5 grid longitude",
        },
        "surface_pressure": {
            "name": "ERA5_SurfacePressure",
            "units": "hPa",
            "dtype": "float32",
            "ndim": 2,
            "long_name": "ERA5 surface pressure",
        },
    }

    def __init__(
        self,
        filename: str,
        *,
        engine: str = "h5netcdf",
        complevel: int = 5,
        logger: Optional[logging.Logger] = None,
    ):
        self.filename = filename
        self.engine = engine
        self.complevel = int(complevel)
        self.log = logger or logging.getLogger(self.__class__.__name__)

        # ERA5 pressure levels in the requested order:
        # 1000 -> 1 hPa
        self.era5_P_levels = np.array(
            [
                1000, 975, 950, 925, 900, 875, 850, 825, 800, 775,
                750, 700, 650, 600, 550, 500, 450, 400, 350, 300,
                250, 225, 200, 175, 150, 125, 100, 70,
                50, 30, 20, 10, 7, 5, 3, 2, 1,
            ],
            dtype=np.int32,
        )

        self._swath_aux_names = {"Latitude", "Longitude", "row_time", "pixel_time", "VZA"}
        self._era5_aux_names = {"ERA5_Latitude", "ERA5_Longitude"}

    def save_mm(
        self,
        *,
        mm_variables: Dict[str, np.ndarray],
        mm_geo: Dict[str, np.ndarray],
        misr_variables: Optional[Dict[str, np.ndarray]] = None,
        modis_variables: Optional[Dict[str, np.ndarray]] = None,
        era5_variables: Optional[Dict[str, np.ndarray]] = None,
        input_files: Optional[List[str]] = None,
        global_attrs: Optional[Dict[str, Any]] = None,
        var_attrs: Optional[Dict[str, Dict[str, Any]]] = None,
        scene_time: Optional[datetime] = None,
    ) -> None:

        for k in ("lat", "lon", "vza"):
            if k not in mm_geo:
                raise KeyError(f"mm_geo missing required key '{k}'")

        lat = np.asarray(mm_geo["lat"])
        lon = np.asarray(mm_geo["lon"])
        vza = np.asarray(mm_geo["vza"])

        if lat.shape != lon.shape:
            raise ValueError("Latitude and Longitude shapes differ")

        ny, nx = lat.shape

        if vza.shape != (ny, nx):
            raise ValueError("VZA shape must match lat/lon")

        ds = xr.Dataset()

        # Swath dimension coordinates
                # Swath dimension coordinates
        ds = ds.assign_coords(
            cell_along_swath_1km=("cell_along_swath_1km", np.arange(ny, dtype=np.int32)),
            cell_cross_swath_1km=("cell_cross_swath_1km", np.arange(nx, dtype=np.int32)),
        )

        ds.coords["cell_along_swath_1km"].attrs.update({
            "long_name": "along-swath index at 1 km",
            "units": "1",
        })
        ds.coords["cell_cross_swath_1km"].attrs.update({
            "long_name": "cross-swath index at 1 km",
            "units": "1",
        })

        enc_hints: Dict[str, Dict[str, Any]] = {}

        # Optional mapping back to original MODIS granule indices
        if "orig_row_index" in mm_geo:
            ori = np.asarray(mm_geo["orig_row_index"])
            if ori.shape != (ny,):
                raise ValueError(f"orig_row_index expected shape {(ny,)}, got {ori.shape}")
            ds = ds.assign_coords(
                orig_row_index=("cell_along_swath_1km", ori.astype(np.int32))
            )
            ds.coords["orig_row_index"].attrs.update({
                "long_name": "row index in original MODIS granule before MISR-overlap subsetting",
                "units": "1",
            })
            enc_hints["orig_row_index"] = {"dtype": "int32", "fill": None}
        if "orig_col_index" in mm_geo:
            oci = np.asarray(mm_geo["orig_col_index"])
            if oci.shape != (nx,):
                raise ValueError(f"orig_col_index expected shape {(nx,)}, got {oci.shape}")
            ds = ds.assign_coords(
                orig_col_index=("cell_cross_swath_1km", oci.astype(np.int32))
            )
            ds.coords["orig_col_index"].attrs.update({
                "long_name": "column index in original MODIS granule before MISR-overlap subsetting",
                "units": "1",
            })
            enc_hints["orig_col_index"] = {"dtype": "int32", "fill": None}

        # Optional row_time
        if "time1d" in mm_geo:
            t1 = np.asarray(mm_geo["time1d"])
            if t1.shape != (ny,):
                raise ValueError(f"time1d expected shape {(ny,)}, got {t1.shape}")

            ds["row_time"] = xr.DataArray(t1, dims=("cell_along_swath_1km",))
            ds["row_time"].attrs.update({
                "standard_name": "time",
                "long_name": "representative observation time for each along-swath row",
                "units": "seconds since 2000-01-01 00:00:00",
                "calendar": "standard",
            })
            enc_hints["row_time"] = self._enc_from_dtype("row_time", "float64")

        # Optional pixel_time
        if "time2d" in mm_geo:
            t2 = np.asarray(mm_geo["time2d"])
            if t2.shape != (ny, nx):
                raise ValueError(f"time2d expected shape {(ny, nx)}, got {t2.shape}")

            ds["pixel_time"] = xr.DataArray(
                t2,
                dims=("cell_along_swath_1km", "cell_cross_swath_1km"),
            )
            ds["pixel_time"].attrs.update({
                "standard_name": "time",
                "long_name": "pixel-level observation time",
                "units": "seconds since 2000-01-01 00:00:00",
                "calendar": "standard",
            })
            enc_hints["pixel_time"] = self._enc_from_dtype("pixel_time", "float64")

        # MM variables
        for key, arr in (mm_variables or {}).items():
            name, attrs, dtype = self._resolve_from_schema(self.MM_SCHEMA, key, arr)
            self._assert_shape_2d(key, arr, (ny, nx))
            ds[name] = xr.DataArray(
                np.asarray(arr),
                dims=("cell_along_swath_1km", "cell_cross_swath_1km"),
            )
            self._apply_attrs(ds, name, attrs, var_attrs)
            enc_hints[name] = self._enc_from_dtype(name, dtype)

        # Swath auxiliary coordinates and geometry
        ds["Latitude"] = xr.DataArray(lat, dims=("cell_along_swath_1km", "cell_cross_swath_1km"))
        ds["Longitude"] = xr.DataArray(lon, dims=("cell_along_swath_1km", "cell_cross_swath_1km"))
        ds["VZA"] = xr.DataArray(vza, dims=("cell_along_swath_1km", "cell_cross_swath_1km"))

        ds["Latitude"].attrs.update({
            "standard_name": "latitude",
            "long_name": "pixel latitude",
            "units": "degrees_north",
        })
        ds["Longitude"].attrs.update({
            "standard_name": "longitude",
            "long_name": "pixel longitude",
            "units": "degrees_east",
        })
        ds["VZA"].attrs.update({
            "long_name": "sensor viewing zenith angle",
            "units": "degrees",
        })

        enc_hints["Latitude"] = self._enc_from_dtype("Latitude", "float32")
        enc_hints["Longitude"] = self._enc_from_dtype("Longitude", "float32")
        enc_hints["VZA"] = self._enc_from_dtype("VZA", "float32")

        # MODIS variables
        if modis_variables:
            for key, arr in modis_variables.items():
                name, attrs, dtype = self._resolve_from_schema(
                    self.MODIS_SCHEMA, key, arr, default_name=key
                )
                self._assert_shape_2d(key, arr, (ny, nx))
                ds[name] = xr.DataArray(
                    np.asarray(arr),
                    dims=("cell_along_swath_1km", "cell_cross_swath_1km"),
                )
                self._apply_attrs(ds, name, attrs, var_attrs)
                enc_hints[name] = self._enc_from_dtype(name, dtype)

        # MISR variables
        if misr_variables:
            for key, arr in misr_variables.items():
                name, attrs, dtype = self._resolve_from_schema(
                    self.MISR_SCHEMA, key, arr, default_name=key
                )
                self._assert_shape_2d(key, arr, (ny, nx))
                ds[name] = xr.DataArray(
                    np.asarray(arr),
                    dims=("cell_along_swath_1km", "cell_cross_swath_1km"),
                )
                self._apply_attrs(ds, name, attrs, var_attrs)
                enc_hints[name] = self._enc_from_dtype(name, dtype)

        # ERA5 variables
        have_era5_latlon = False
        era5_yx_defined = False

        if era5_variables:
            # First pass to determine ERA5 horizontal shape
            for key, arr in era5_variables.items():
                a = np.asarray(arr)
                if a.ndim == 3:
                    _, ey, ex = a.shape
                elif a.ndim == 2:
                    ey, ex = a.shape
                else:
                    raise ValueError(f"{key}: expected 2D or 3D array, got ndim={a.ndim}")

                if not era5_yx_defined:
                    ds = ds.assign_coords(
                        z=("z", self.era5_P_levels),
                        era5_cell_along_swath_10km=("era5_cell_along_swath_10km", np.arange(ey, dtype=np.int32)),
                        era5_cell_cross_swath_10km=("era5_cell_cross_swath_10km", np.arange(ex, dtype=np.int32)),
                    )
                    era5_yx_defined = True
                    break

            # Explicit z coordinate variable with CF attrs
            if "z" in ds.coords:
                ds.coords["z"].attrs.update({
                    "standard_name": "air_pressure",
                    "long_name": "pressure level",
                    "units": "hPa",
                    "positive": "down",
                    "axis": "Z",
                })
                # add encoding hint so z is written explicitly
                enc_hints["z"] = {"dtype": "int32", "fill": None}

            if "era5_cell_along_swath_10km" in ds.coords:
                ds.coords["era5_cell_along_swath_10km"].attrs.update({
                    "long_name": "along-swath index at 10 km",
                    "units": "1",
                })
                enc_hints["era5_cell_along_swath_10km"] = {"dtype": "int32", "fill": None}

            if "era5_cell_cross_swath_10km" in ds.coords:
                ds.coords["era5_cell_cross_swath_10km"].attrs.update({
                    "long_name": "cross-swath index at 10 km",
                    "units": "1",
                })
                enc_hints["era5_cell_cross_swath_10km"] = {"dtype": "int32", "fill": None}

            # Second pass to write ERA5 variables
            for key, arr in era5_variables.items():
                a = np.asarray(arr)
                schema = self.ERA5_SCHEMA.get(
                    key,
                    {"name": key, "dtype": "float32", "ndim": a.ndim},
                )

                name = schema["name"]
                dtype = schema.get("dtype", "float32")
                ndim = schema.get("ndim", a.ndim)
                attrs = {k: v for k, v in schema.items() if k not in ("name", "dtype", "ndim")}

                if ndim == 3 and a.ndim == 3:
                    nz, ey, ex = a.shape
                    if nz != len(self.era5_P_levels):
                        raise ValueError(
                            f"{name}: z must have length {len(self.era5_P_levels)}, got {nz}"
                        )

                    ds[name] = xr.DataArray(
                        a,
                        dims=("z", "era5_cell_along_swath_10km", "era5_cell_cross_swath_10km"),
                    )

                elif ndim == 2 and a.ndim == 2:
                    ds[name] = xr.DataArray(
                        a,
                        dims=("era5_cell_along_swath_10km", "era5_cell_cross_swath_10km"),
                    )

                else:
                    raise ValueError(f"{key}: schema ndim={ndim}, but array ndim={a.ndim}")

                self._apply_attrs(ds, name, attrs, var_attrs)
                enc_hints[name] = self._enc_from_dtype(name, dtype)

                if name in ("ERA5_Latitude", "ERA5_Longitude"):
                    have_era5_latlon = True

        # Optional scalar scene time coordinate
        if scene_time is not None:
            if scene_time.tzinfo is None:
                scene_time = scene_time.replace(tzinfo=timezone.utc)
            else:
                scene_time = scene_time.astimezone(timezone.utc)

            ds = ds.assign_coords(time=np.array(np.datetime64(scene_time.replace(tzinfo=None))))
            ds["time"].attrs.update({
                "standard_name": "time",
                "long_name": "scene time",
            })
            ds["time"].encoding.update({
                "units": "seconds since 1970-01-01 00:00:00",
                "calendar": "standard",
            })

        # Global attributes
        ds.attrs = self._build_global_attrs(global_attrs or {}, lat, lon, input_files)

        # Attach coordinates attributes only to true data variables
        self._attach_swath_coordinates(ds)
        if era5_variables:
            self._attach_era5_coordinates(ds, have_era5_latlon=have_era5_latlon)

        # Ancillary variables
        if "MM_CloudTopHeight" in ds and "MM_Flag" in ds:
            ds["MM_CloudTopHeight"].attrs.setdefault("ancillary_variables", "MM_Flag")

        if "MODIS_CloudTopHeight" in ds and "MODIS_CloudPhase" in ds:
            ds["MODIS_CloudTopHeight"].attrs.setdefault("ancillary_variables", "MODIS_CloudPhase")

        # Ensure _FillValue is only in encoding for data vars
        for v in ds.data_vars:
            ds[v].attrs.pop("_FillValue", None)

        # Encoding
        base = {"zlib": True, "shuffle": True, "complevel": self.complevel}
        encoding: Dict[str, Dict[str, Any]] = {}

        # Data variables
        for name, hint in enc_hints.items():
            entry: Dict[str, Any] = {**base, "dtype": hint["dtype"]}
            if hint["fill"] is not None:
                entry["_FillValue"] = hint["fill"]
            encoding[name] = entry

        # Coordinate variables: no fill value, and keep them explicit
        for cname in ds.coords:
            if cname not in encoding:
                dt = str(ds.coords[cname].dtype)
                encoding[cname] = {"dtype": dt, "_FillValue": None}

        # Range guard for integer data variable storage
        for name, da in ds.data_vars.items():
            if name not in encoding:
                continue

            dt = np.dtype(encoding[name]["dtype"])
            if np.issubdtype(dt, np.integer):
                vals = np.asarray(da.values)

                if np.issubdtype(vals.dtype, np.floating):
                    finite = np.isfinite(vals)
                    if not np.any(finite):
                        continue
                    vmin = np.nanmin(vals)
                    vmax = np.nanmax(vals)
                else:
                    if vals.size == 0:
                        continue
                    vmin = vals.min()
                    vmax = vals.max()

                info = np.iinfo(dt)
                if vmin < info.min or vmax > info.max:
                    self.log.warning(
                        f"{name}: value range [{vmin}, {vmax}] exceeds storage dtype {dt}"
                    )

        ds.to_netcdf(
            self.filename,
            engine=self.engine,
            format="NETCDF4",
            encoding=encoding,
        )
        self.log.info(f"Saved {self.filename}")

    def _resolve_from_schema(
        self,
        schema: Dict[str, Dict[str, Any]],
        key: str,
        arr: np.ndarray,
        default_name: Optional[str] = None,
    ):
        meta = schema.get(key)
        if meta:
            name = meta.get("name", key)
            dtype = meta.get("dtype", self._guess_dtype(name, arr))
            attrs = {k: v for k, v in meta.items() if k not in ("name", "dtype")}
        else:
            name = default_name or key
            dtype = self._guess_dtype(name, arr)
            attrs = {}
        return name, attrs, dtype

    def _apply_attrs(
        self,
        ds: xr.Dataset,
        name: str,
        attrs: Dict[str, Any],
        var_overrides: Optional[Dict[str, Dict[str, Any]]],
    ):
        clean_attrs = {k: v for k, v in attrs.items() if k != "_FillValue"}
        ds[name].attrs.update(clean_attrs)

        if var_overrides and name in var_overrides:
            ds[name].attrs.update({
                k: v
                for k, v in var_overrides[name].items()
                if k not in ("dtype", "_FillValue")
            })

    def _enc_from_dtype(self, name: str, dtype: str) -> Dict[str, Any]:
        dt = np.dtype(str(dtype))

        if np.issubdtype(dt, np.integer):
            lname = name.lower()
            if "flag" in lname or "phase" in lname or "method" in lname or lname.endswith("_qi"):
                fill = np.iinfo(dt).min
            else:
                info = np.iinfo(dt)
                preferred = -9999
                fill = preferred if info.min <= preferred <= info.max else info.min
            return {"dtype": str(dtype), "fill": np.array(fill, dtype=dt)}

        return {"dtype": str(dtype), "fill": np.array(-9999.0, dtype=dt)}

    def _guess_dtype(self, name: str, arr: np.ndarray) -> str:
        lname = name.lower()
        if "flag" in lname or "phase" in lname or "method" in lname:
            return "int8"

        a = np.asarray(arr)
        if np.issubdtype(a.dtype, np.integer):
            return "int16"
        return "float32"

    def _assert_shape_2d(self, key: str, arr: np.ndarray, expected: tuple[int, int]):
        if np.asarray(arr).shape != expected:
            raise ValueError(f"{key}: expected shape {expected}, got {np.asarray(arr).shape}")

    def _build_global_attrs(
        self,
        base: Dict[str, Any],
        lat2d: np.ndarray,
        lon2d: np.ndarray,
        input_files: Optional[List[str]],
    ) -> Dict[str, Any]:
        out = dict(base)

        now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

        out.setdefault("Conventions", "CF-1.10")
        out.setdefault("title", "Fused MISR and MODIS Cloud Properties")
        out.setdefault("institution", "University of Illinois at Urbana-Champaign")
        out.setdefault("source", "MISR TC_Cloud plus MODIS MOD06 and MOD03 plus ERA5 reanalysis")
        out.setdefault("history", f"{now}: created on keeling.earth.illinois.edu")
        out.setdefault("references", "https://cfconventions.org/")
        out["date_created"] = now

        starting_time = self._infer_starting_time_from_filename(self.filename)
        if starting_time is not None:
            out["starting_time"] = starting_time

        with np.errstate(invalid="ignore"):
            lat_min = float(np.nanmin(lat2d))
            lat_max = float(np.nanmax(lat2d))
            lon_min = float(np.nanmin(lon2d))
            lon_max = float(np.nanmax(lon2d))

        out["geospatial_lat_min"] = lat_min
        out["geospatial_lat_max"] = lat_max
        out["geospatial_lon_min"] = lon_min
        out["geospatial_lon_max"] = lon_max

        out["NORTHBOUNDINGCOORDINATE"] = lat_max
        out["SOUTHBOUNDINGCOORDINATE"] = lat_min
        out["WESTBOUNDINGCOORDINATE"] = lon_min
        out["EASTBOUNDINGCOORDINATE"] = lon_max

        if input_files:
            for i, p in enumerate(input_files[:5]):
                out[f"input_file_{i + 1:03d}"] = str(p)

        return out

    def _infer_starting_time_from_filename(self, path: str) -> Optional[str]:
        s = str(path)
        m = re.search(r"_A(\d{4})(\d{3})\.(\d{2})(\d{2})_", s)
        if not m:
            return None

        year = int(m.group(1))
        doy = int(m.group(2))
        hour = int(m.group(3))
        minute = int(m.group(4))

        dt = datetime(year, 1, 1, tzinfo=timezone.utc) + timedelta(days=doy - 1)
        dt = dt.replace(hour=hour, minute=minute, second=0, microsecond=0)

        return dt.strftime("%Y-%m-%dT%H:%M:%SZ")

    def _attach_swath_coordinates(self, ds: xr.Dataset) -> None:
        """
        Attach coordinates attribute to 2D swath data variables.

        Preferred order:
        - Always: Latitude Longitude
        - If pixel_time exists: add pixel_time
        - Else if row_time exists: add row_time
        - Optional scalar scene time may also be added
        """
        parts = ["Latitude", "Longitude"]

        if "pixel_time" in ds:
            parts.append("pixel_time")
        elif "row_time" in ds:
            parts.append("row_time")
        if "time" in ds.coords:
            parts.append("time")
        coord_string = " ".join(parts)

        for v in ds.data_vars:
            if v in self._swath_aux_names:
                continue
            if ds[v].dims == ("cell_along_swath_1km", "cell_cross_swath_1km"):
                ds[v].attrs["coordinates"] = coord_string

    def _attach_era5_coordinates(self, ds: xr.Dataset, *, have_era5_latlon: bool) -> None:
        for v in ds.data_vars:
            if v in self._era5_aux_names:
                continue

            dims = ds[v].dims

            if dims == ("era5_cell_along_swath_10km", "era5_cell_cross_swath_10km"):
                parts = []
                if have_era5_latlon:
                    parts.extend(["ERA5_Latitude", "ERA5_Longitude"])
                if "time" in ds.coords:
                    parts.append("time")
                if parts:
                    ds[v].attrs["coordinates"] = " ".join(parts)

            elif dims == ("z", "era5_cell_along_swath_10km", "era5_cell_cross_swath_10km"):
                parts = ["z"]
                if have_era5_latlon:
                    parts.extend(["ERA5_Latitude", "ERA5_Longitude"])
                if "time" in ds.coords:
                    parts.append("time")
                ds[v].attrs["coordinates"] = " ".join(parts)
