#!/usr/bin/env python3
"""
xarray_saver.py — minimal-CF NetCDF writer for MODIS/MISR/ERA5 fused scenes.

Key CF choices:
- Conventions: CF-1.10
- Swath variables remain on index dims ('cell_along_swath_1km','cell_cross_swath_1km')
  with 2-D auxiliary coordinates Latitude/Longitude and 'coordinates' attributes.
- ERA5 variables use a pressure vertical coordinate 'z' (units=hPa, positive=down).
- Optional scalar 'time' coordinate can be included.
- Integer flags use _FillValue outside valid flag_values and include flag_meanings.
"""

from __future__ import annotations

import logging
from datetime import datetime
from typing import Optional, Dict, Any, List

import numpy as np
import xarray as xr


class XarraySaver:
    """
    Simple NetCDF writer with a single, easy-to-edit schema for names & attributes.

    Dims:
      - MM/MODIS/MISR 2D: ('cell_along_swath_1km', 'cell_cross_swath_1km')
      - ERA5 2D: ('era5_cell_along_swath_10km', 'era5_cell_cross_swath_10km')
      - ERA5 3D: ('z', 'era5_cell_along_swath_10km', 'era5_cell_cross_swath_10km') with z in hPa
    """

    # ---------- ONE-STOP SCHEMA (edit here) ----------
    MM_SCHEMA: Dict[str, Dict[str, Any]] = {
        'cth':   {'name': 'MM_CloudTopHeight',            'units': 'm',   'dtype': 'int16',
                  'valid_range': [0, 20000], 'long_name': 'multi-mission cloud top height'},
        'ctp':   {'name': 'MM_CloudTopPressure',          'units': 'hPa', 'dtype': 'float32',
                  'long_name': 'multi-mission cloud top pressure'},
        'emissivity': {'name': 'MM_CloudEffectiveEmissivity', 'units': '1', 'dtype': 'float32',
                       'long_name': 'multi-mission cloud effective emissivity'},
        'opt':   {'name': 'MMCloudOpticalDepth',         'units': '1',   'dtype': 'float32',
                  'long_name': 'multi-mission cloud optical depth'},

        'cflag': {'name': 'MM_Flag', 'units': '1', 'dtype': 'int8',
                  'flag_values': np.array([0, 1, 2, 3, 4, 5, 6], dtype='int8'),
                  'flag_meanings': 'clear low_cloud mid_cloud high_cloud multilayer ambiguous invalid',
                  'long_name': 'cloud classification/quality flag'},
        'dflag': {'name': 'MM_DiagFlag', 'units': '1', 'dtype': 'int8',
                  'flag_values': np.array([0, 1, 2, 3], dtype='int8'),
                  'flag_meanings': 'ok radiance_qc_failed geometry_qc_failed other_diagnostic',
                  'long_name': 'diagnostic flag'}
    }

    MODIS_SCHEMA: Dict[str, Dict[str, Any]] = {
        'cth':   {'name': 'MODIS_CloudTopHeight',        'units': 'm',   'dtype': 'int16',
                  'long_name': 'MODIS cloud top height'},
        'ctp':   {'name': 'MODIS_CloudTopPressure',      'units': 'hPa', 'dtype': 'float32',
                  'long_name': 'MODIS cloud top pressure'},
        'emissivity': {'name': 'MODIS_CloudEmissivity',  'units': '1',   'dtype': 'float32',
                       'long_name': 'MODIS cloud emissivity'},
        'opt':   {'name': 'MODIS_CloudOpticalDepth',     'units': '1',   'dtype': 'float32',
                  'long_name': 'MODIS cloud optical depth'},
        'cphase':{'name': 'MODIS_CloudPhase',            'units': '1',   'dtype': 'int8',
                  'flag_values': np.array([0, 1, 2, 3], dtype='int8'),
                  'flag_meanings': 'water ice mixed unknown',
                  'long_name': 'MODIS cloud phase flag'}
    }

    MISR_SCHEMA: Dict[str, Dict[str, Any]] = {
        'misrcth': {'name': 'MISR_CloudTopHeight', 'units': 'm', 'dtype': 'int16',
                    'long_name': 'MISR cloud top height without wind correction'},
        'misrcth_qa': {'name': 'MISR_CloudTopHeight_QI', 'units': 'm', 'dtype': 'int8',
                    'long_name': 'MISR cloud top height Qualty Indicator',
                    'flag_meaning': '0: Worst Quality, 100: Best Quatliy, 128: Bad Data',
                    },
        
    }

    ERA5_SCHEMA: Dict[str, Dict[str, Any]] = {
        # 3D inputs (z,y,x)
        'geopotential_org': {'name': 'ERA5_Geopotential', 'units': '"m**2 s**-2', 'dtype': 'float32', 'ndim': 3,
                             'long_name': 'ERA5 Geopotential'},
        'temperature_org':  {'name': 'ERA5_Temperature',  'units': 'K',      'dtype': 'float32', 'ndim': 3,
                             'long_name': 'ERA5 temperature'},
        'mixingratio_org':  {'name': 'ERA5_Humidity',     'units': 'kg kg**-1',      'dtype': 'float32', 'ndim': 3,
                             'long_name': 'ERA5 water vapor mixing ratio'},
        'u': {'name': 'ERA5_U','units': 'm s**-1',  'dtype': 'float32', 'ndim': 3,
                             'long_name': 'ERA5 U component of wind)', "short_name": 'eastward_wind'},
        'v': {'name': 'ERA5_V','units': 'm s**-1',  'dtype': 'float32', 'ndim': 3,
                             'long_name': 'ERA5 V component of wind',"short_name": 'northward_wind' },
        'w': {'name': 'ERA5_W','units': 'm s**-1',  'dtype': 'float32', 'ndim': 3,
                             'long_name': 'ERA5 Vertical velocity',"short_name": 'w'},



        # 2D inputs (y,x)
        'ERA5Humidity':     {'name': 'ERA5_Humidity',         'units': '1', 'dtype': 'float32', 'ndim': 2,
                             'long_name': 'ERA5 water vapor mixing ratio at surface/2m'},
        'dew2m':            {'name': 'ERA5_DewPoint2m',       'units': 'K', 'dtype': 'float32', 'ndim': 2,
                             'long_name': 'ERA5 2m dew point temperature'},
        'sst':              {'name': 'ERA5_SeaSurfaceTemperature','units':'K','dtype':'float32','ndim':2,
                             'long_name': 'ERA5 sea surface temperature'},
        'temp2m':           {'name': 'ERA5_Temperature2m',    'units': 'K', 'dtype': 'float32', 'ndim': 2,
                             'long_name': 'ERA5 2m air temperature'},
        'skint':            {'name': 'ERA5_SkinTemperature',  'units': 'K', 'dtype': 'float32', 'ndim': 2,
                             'long_name': 'ERA5 skin temperature'},
        'msp':              {'name': 'ERA5_MeanSeaLevelPressure','units':'hPa','dtype':'float32','ndim':2,
                             'long_name': 'ERA5 mean sea level pressure'},
        'swmr':             {'name': 'ERA5_2mWaterMixingRatio','units':'1', 'dtype':'float32','ndim':2,
                             'long_name': 'ERA5 2m water vapor mixing ratio'},
        'ERA5Latitude':     {'name': 'ERA5_Latitude', 'units': 'degrees_north', 'dtype': 'float32', 'ndim': 2,
                             'standard_name': 'latitude',  'long_name': 'ERA5 grid latitude'},
        'ERA5Longitude':    {'name': 'ERA5_Longitude','units': 'degrees_east', 'dtype': 'float32', 'ndim': 2,
                             'standard_name': 'longitude', 'long_name': 'ERA5 grid longitude'},
        'surface_pressure': {'name': 'ERA5_SurfacePressure','units': 'hPa',  'dtype': 'float32', 'ndim': 2,
                             'long_name': 'ERA5 surface pressure (replicated to levels for convenience)'},
    }

    def __init__(self, filename: str, *, engine: str = "h5netcdf", complevel: int = 5,
                 logger: Optional[logging.Logger] = None):
        self.filename  = filename
        self.engine    = engine
        self.complevel = int(complevel)
        self.log       = logger or logging.getLogger(self.__class__.__name__)
        # Pressure levels (hPa), increasing downward → positive="down"
        self.era5_P_levels = np.array(
            [1, 2, 3, 5, 7, 10, 20, 30, 50, 70,
             100, 125, 150, 175, 200, 225, 250, 300, 350, 400,
             450, 500, 550, 600, 650, 700, 750, 775, 800, 825,
             850, 875, 900, 925, 950, 975, 1000], dtype=np.int32
        )

    # -------------- public API --------------
    def save_mm(
        self,
        *,  # MUST USE KEYWORD-ONLY PARAMETERS
        mm_variables: Dict[str, np.ndarray],
        mm_geo: Dict[str, np.ndarray],  # requires keys: 'lat', 'lon', 'vza'
        misr_variables: Optional[Dict[str, np.ndarray]] = None,
        modis_variables: Optional[Dict[str, np.ndarray]] = None,
        era5_variables: Optional[Dict[str, np.ndarray]] = None,
        input_files: Optional[List[str]] = None,
        global_attrs: Optional[Dict[str, Any]] = None,
        var_attrs: Optional[Dict[str, Dict[str, Any]]] = None,  # optional per-var overrides
        scene_time: Optional[datetime] = None
    ) -> None:

        # ---- sanity / shapes ----
        for k in ("lat", "lon", "vza"):
            if k not in mm_geo:
                raise KeyError(f"mm_geo missing '{k}'")
        lat, lon, vza = mm_geo["lat"], mm_geo["lon"], mm_geo["vza"]
        if lat.shape != lon.shape:
            raise ValueError("Latitude and Longitude shapes differ.")
        ny, nx = lat.shape
        if vza.shape != (ny, nx):
            raise ValueError("VZA shape must match lat/lon.")

        ds = xr.Dataset()

        # Swath index coords (dimension variables)
        ds = ds.assign_coords(
            cell_along_swath_1km=("cell_along_swath_1km", np.arange(ny, dtype=np.int32)),
            cell_cross_swath_1km=("cell_cross_swath_1km", np.arange(nx, dtype=np.int32)),
        )
        ds.coords["cell_along_swath_1km"].attrs.update({"long_name": "along-swath index at 1 km", "units": "1"})
        ds.coords["cell_cross_swath_1km"].attrs.update({"long_name": "cross-swath index at 1 km", "units": "1"})
     
        enc_hints: Dict[str, Dict[str, Any]] = {}  # var name -> {'dtype':..., 'fill':...}

        # --- MM (generic loop, uses schema) ---
        for key, arr in (mm_variables or {}).items():
            name, attrs, dtype = self._resolve_from_schema(self.MM_SCHEMA, key, arr)
            self._assert_shape_2d(key, arr, (ny, nx))
            ds[name] = xr.DataArray(arr, dims=("cell_along_swath_1km", "cell_cross_swath_1km"))
            self._apply_attrs(ds, name, attrs, var_attrs)
            enc_hints[name] = self._enc_from_dtype(name, dtype)

        # --- GEO: 2-D lat/lon auxiliary coordinates + VZA ---
        ds["Latitude"]  = xr.DataArray(lat, dims=("cell_along_swath_1km", "cell_cross_swath_1km"))
        ds["Longitude"] = xr.DataArray(lon, dims=("cell_along_swath_1km", "cell_cross_swath_1km"))
        ds["VZA"]       = xr.DataArray(vza, dims=("cell_along_swath_1km", "cell_cross_swath_1km"))

        ds["Latitude"].attrs.update({"standard_name": "latitude",  "units": "degrees_north"})
        ds["Longitude"].attrs.update({"standard_name": "longitude", "units": "degrees_east"})
        ds["VZA"].attrs.setdefault("long_name", "sensor_zenith_angle")
        ds["VZA"].attrs.setdefault("units", "degree")

        enc_hints["Latitude"]  = self._enc_from_dtype("Latitude",  "float32")
        enc_hints["Longitude"] = self._enc_from_dtype("Longitude", "float32")
        enc_hints["VZA"]       = self._enc_from_dtype("VZA",       "float32")

        # --- MODIS (optional) ---
        if modis_variables:
            for key, arr in modis_variables.items():
                name, attrs, dtype = self._resolve_from_schema(self.MODIS_SCHEMA, key, arr, default_name=key)
                self._assert_shape_2d(key, arr, (ny, nx))
                ds[name] = xr.DataArray(arr, dims=("cell_along_swath_1km", "cell_cross_swath_1km"))
                self._apply_attrs(ds, name, attrs, var_attrs)
                enc_hints[name] = self._enc_from_dtype(name, dtype)

        # --- MISR (optional) ---
        if misr_variables:
            for key, arr in misr_variables.items():
                name, attrs, dtype = self._resolve_from_schema(self.MISR_SCHEMA, key, arr, default_name=key)
                self._assert_shape_2d(key, arr, (ny, nx))
                ds[name] = xr.DataArray(arr, dims=("cell_along_swath_1km", "cell_cross_swath_1km"))
                self._apply_attrs(ds, name, attrs, var_attrs)
                enc_hints[name] = self._enc_from_dtype(name, dtype)

        # --- ERA5 (2D/3D, label z) ---
        have_era5_latlon = False
        if era5_variables:
            for key, arr in era5_variables.items():
                a = np.asarray(arr)
                schema = self.ERA5_SCHEMA.get(key, {'name': key, 'units': '1', 'dtype': 'float32', 'ndim': a.ndim})
                name, dtype, ndim = schema['name'], schema.get('dtype', 'float32'), schema.get('ndim', a.ndim)
                attrs = {k: v for k, v in schema.items() if k not in ('name', 'dtype', 'ndim')}

                if ndim == 3 and a.ndim == 3:
                    nz, ey, ex = a.shape
                    if nz != len(self.era5_P_levels):
                        raise ValueError(f"{name}: z must have length {len(self.era5_P_levels)}; got {nz}")
                    ds[name] = xr.DataArray(
                        a,
                        dims=("z", "era5_cell_along_swath_10km", "era5_cell_cross_swath_10km"),
                        coords={
                            "z": self.era5_P_levels,
                            "era5_cell_along_swath_10km": np.arange(ey, dtype=np.int32),
                            "era5_cell_cross_swath_10km": np.arange(ex, dtype=np.int32),
                        },
                    )
                elif ndim == 2 and a.ndim == 2:
                    ey, ex = a.shape
                    ds[name] = xr.DataArray(
                        a,
                        dims=("era5_cell_along_swath_10km", "era5_cell_cross_swath_10km"),
                        coords={
                            "era5_cell_along_swath_10km": np.arange(ey, dtype=np.int32),
                            "era5_cell_cross_swath_10km": np.arange(ex, dtype=np.int32),
                        },
                    )
                else:
                    # ignore unexpected dims
                    continue

                self._apply_attrs(ds, name, attrs, var_attrs)
                enc_hints[name] = self._enc_from_dtype(name, dtype)

                if name in ("ERA5Latitude", "ERA5Longitude"):
                    have_era5_latlon = True

            # Add attrs to pressure coord if present
            if "z" in ds.coords:
                ds["z"].attrs.update({
                    "standard_name": "air_pressure",
                    "long_name": "pressure_level",
                    "units": "hPa",
                    "positive": "down",
                    "axis": "Z"
                })
            ds.coords["era5_cell_along_swath_10km"].attrs.update({"long_name": "along-swath index at 10 km", "units": "1"})
            ds.coords["era5_cell_cross_swath_10km"].attrs.update({"long_name": "cross-swath index at 10 km", "units": "1"})


        # --- optional scalar time coordinate ---
        if scene_time is not None:
            ds = ds.assign_coords(time=np.array(np.datetime64(scene_time)))
            ds["time"].attrs.update({"standard_name": "time", "long_name": "scene_time"})
            # Let xarray handle CF time encoding; explicit encoding is optional:
            ds["time"].encoding.update({
                "units": "seconds since 1970-01-01 00:00:00",
                "calendar": "standard"
            })

        # --- global attributes (CF/ACDD-minimum) ---
        ds.attrs = self._build_global_attrs(global_attrs or {}, lat, lon, input_files)

        # --- CRS (optional but helpful for GIS/CF tooling) ---
        ds["crs"] = xr.DataArray(0)
        ds["crs"].attrs.update({
            "grid_mapping_name": "latitude_longitude",
            "long_name": "WGS84 geographic coordinate reference system",
            "semi_major_axis": 6378137.0,
            "inverse_flattening": 298.257223563
        })

        # --- attach coordinates/grid_mapping to variables ---
        self._attach_swath_coordinates(ds)  # adds 'coordinates' to (ny,nx) vars
        if era5_variables:
            self._attach_era5_coordinates(ds, have_era5_latlon=have_era5_latlon)
        for v in ds.data_vars:
            ds[v].attrs.setdefault("grid_mapping", "crs")

        # Link flags as ancillary_variables if present
        if "MMCloudTopHeight" in ds and "MMFlag" in ds:
            av = "MMFlag"
            if "MMDiagFlag" in ds:
                av += " MMDiagFlag"
            ds["MMCloudTopHeight"].attrs.setdefault("ancillary_variables", av)
        if "MODISCloudTopHeight" in ds and "MODISCloudPhase" in ds:
            ds["MODISCloudTopHeight"].attrs.setdefault("ancillary_variables", "MODISCloudPhase")

        # --- ensure no stray _FillValue in attrs (we set via encoding only) ---
        for v in ds.data_vars:
            ds[v].attrs.pop("_FillValue", None)

        # --- encoding & write ---
        base = {'zlib': True, 'shuffle': True, 'complevel': self.complevel}
        encoding = {name: {**base, 'dtype': h['dtype'], '_FillValue': h['fill']}
                    for name, h in enc_hints.items()}

        # int-range guard
        for name, da in ds.data_vars.items():
            if name not in encoding:
                continue
            dt = np.dtype(encoding[name]['dtype'])
            if np.issubdtype(dt, np.integer):
                with np.errstate(invalid="ignore"):
                    vmin = np.nanmin(da.values)  # may warn if all-NaN; ignore
                    vmax = np.nanmax(da.values)
                info = np.iinfo(dt)
                if np.isfinite(vmin) and np.isfinite(vmax) and (vmin < info.min or vmax > info.max):
                    self.log.warning(f"{name}: value range [{vmin}, {vmax}] exceeds {dt}")

        ds.to_netcdf(self.filename, engine=self.engine, format="NETCDF4", encoding=encoding)
        self.log.info(f"Saved {self.filename}")

    # -------------- small helpers --------------
    def _resolve_from_schema(self, schema: Dict[str, Dict[str, Any]], key: str, arr: np.ndarray,
                             default_name: Optional[str] = None):
        """Return (output_name, attrs, dtype) using schema entry if present, else sensible defaults."""
        meta = schema.get(key)
        if meta:
            name = meta.get('name', key)
            dtype = meta.get('dtype', self._guess_dtype(name, arr))
            attrs = {k: v for k, v in meta.items() if k not in ('name', 'dtype')}
        else:
            name = default_name or key
            dtype = self._guess_dtype(name, arr)
            attrs = {}
        return name, attrs, dtype

    def _apply_attrs(self, ds: xr.Dataset, name: str, attrs: Dict[str, Any],
                     var_overrides: Optional[Dict[str, Dict[str, Any]]]):
        # schema attrs (except _FillValue)
        ds[name].attrs.update({k: v for k, v in attrs.items() if k != "_FillValue"})
        # optional overrides (except dtype/_FillValue)
        if var_overrides and name in var_overrides:
            ds[name].attrs.update({k: v for k, v in var_overrides[name].items()
                                   if k not in ('_FillValue', 'dtype')})

    def _enc_from_dtype(self, name: str, dtype: str) -> Dict[str, Any]:
        """
        Decide on storage dtype and a safe _FillValue that fits that dtype.
        - flags -> dtype minimum (e.g., -128 for int8) so it never collides with valid flag_values
        - other integers -> -9999 if it fits; otherwise use dtype minimum
        - floats -> -9999.0
        """
        dt = np.dtype(str(dtype))
        if np.issubdtype(dt, np.integer):
            if "flag" in name.lower():
                fill = np.iinfo(dt).min  # avoid 0 which is a valid flag
            else:
                info = np.iinfo(dt)
                preferred = -9999
                fill = preferred if info.min <= preferred <= info.max else info.min
            return {"dtype": str(dtype), "fill": np.array(fill, dtype=dt)}
        # float-like
        return {"dtype": str(dtype), "fill": np.array(-9999.0, dtype=dt)}

    def _guess_dtype(self, name: str, arr: np.ndarray) -> str:
        if "flag" in name.lower():
            return "int8"
        return "int16" if np.issubdtype(np.asarray(arr).dtype, np.integer) else "float32"

    def _assert_shape_2d(self, key: str, arr: np.ndarray, expected: tuple[int, int]):
        if np.asarray(arr).shape != expected:
            raise ValueError(f"{key}: expected shape {expected}, got {np.asarray(arr).shape}")

    def _build_global_attrs(self, base: Dict[str, Any], lat2d: np.ndarray, lon2d: np.ndarray,
                            input_files: Optional[List[str]]) -> Dict[str, Any]:
        out = dict(base)
        # CF/ACDD-minimum discovery attrs
        out.setdefault('Conventions', 'CF-1.10')
        out.setdefault('title', 'Fused MISR and MODIS Cloud Properties')
        out.setdefault('institution', 'University of Illinois at Urbana-Champaign')
        out.setdefault('source', 'MISR AM1 + MODIS MOD06/MOD03 + ERA5 reanalysis')
        out.setdefault('history', f"{datetime.utcnow().strftime('%Y-%m-%dT%H:%M:%SZ')}: created by XarraySaver")
        out.setdefault('references', 'https://cfconventions.org/')
        out['date_created'] = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")

        # Nice-to-have ACDD geospatial extents (not required but harmless/handy)
        with np.errstate(invalid="ignore"):
            out['geospatial_lat_min'] = float(np.nanmin(lat2d))
            out['geospatial_lat_max'] = float(np.nanmax(lat2d))
            out['geospatial_lon_min'] = float(np.nanmin(lon2d))
            out['geospatial_lon_max'] = float(np.nanmax(lon2d))

        # Retain your legacy bounds (OK to include alongside ACDD)
        out['NORTHBOUNDINGCOORDINATE'] = out['geospatial_lat_min']
        out['SOUTHBOUNDINGCOORDINATE'] = out['geospatial_lat_max']
        out['WESTBOUNDINGCOORDINATE']  = out['geospatial_lon_min']
        out['EASTBOUNDINGCOORDINATE']  = out['geospatial_lon_max']

        if input_files:
            for i, p in enumerate(input_files[0:5]):
                out[f"input_file_{i+1:03d}"] = str(p)
        return out

    def _attach_swath_coordinates(self, ds: xr.Dataset) -> None:
        """Attach 'coordinates' attribute to all 2D swath variables."""
        coo = "Latitude Longitude"
        if "time" in ds.coords:
            coo += " time"
        for v in ds.data_vars:
            if ds[v].dims == ("cell_along_swath_1km", "cell_cross_swath_1km"):
                ds[v].attrs.setdefault("coordinates", coo)

    def _attach_era5_coordinates(self, ds: xr.Dataset, *, have_era5_latlon: bool) -> None:
        """Attach 'coordinates' attribute to ERA5 variables (2D and 3D)."""
        for v in ds.data_vars:
            dims = ds[v].dims
            if dims == ("era5_cell_along_swath_10km", "era5_cell_cross_swath_10km"):
                parts = []
                if have_era5_latlon:
                    parts = ["ERA5Latitude", "ERA5Longitude"]
                if "time" in ds.coords:
                    parts.append("time")
                if parts:
                    ds[v].attrs.setdefault("coordinates", " ".join(parts))
            elif dims == ("z", "era5_cell_along_swath_10km", "era5_cell_cross_swath_10km"):
                parts = ["z"]
                if have_era5_latlon:
                    parts += ["ERA5Latitude", "ERA5Longitude"]
                if "time" in ds.coords:
                    parts.append("time")
                ds[v].attrs.setdefault("coordinates", " ".join(parts))