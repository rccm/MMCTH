import logging
from datetime import datetime, timezone

import numpy as np
import xarray as xr

from src.data_processors.process_main import MainProcessor
from src.data_processors.xarray_saver import XarraySaver
from src.data_readers.mod_read import MOD06Granule


def test_decode_snow_ice_mask_preserves_source_polarity_and_missing_values():
    raw_u8 = np.array([[0b11011111, 0b11111111, 0]], dtype=np.uint8)
    raw_i8 = raw_u8.view(np.int8)

    decoded = MOD06Granule._decode_snow_ice_first_byte(raw_i8, fill_value=0)

    np.testing.assert_array_equal(decoded, np.array([[0, 1, 2]], dtype=np.int8))


def test_decode_snow_ice_mask_requires_a_determined_cloud_mask():
    # The last byte has bit 5 set, but bit 0 is clear: snow/ice is unknown.
    raw = np.array([[0b00000001, 0b00100001, 0b00100000]], dtype=np.uint8)

    decoded = MOD06Granule._decode_snow_ice_first_byte(raw, fill_value=np.uint8(0))
    cloud = MOD06Granule._decode_cloud_mask_first_byte(raw, fill_value=np.uint8(0))

    np.testing.assert_array_equal(decoded, np.array([[0, 1, 2]], dtype=np.int8))
    np.testing.assert_array_equal(cloud, np.array([[0, 0, 4]], dtype=np.int8))


def test_decode_snow_ice_mask_uses_unknown_for_explicit_source_missing():
    raw = np.array([[0b00100001]], dtype=np.uint8)

    decoded_from_mask = MOD06Granule._decode_snow_ice_first_byte(
        raw,
        source_missing=np.array([[True]]),
    )
    decoded_from_masked_array = MOD06Granule._decode_snow_ice_first_byte(
        np.ma.array(raw, mask=[[True]]),
    )

    np.testing.assert_array_equal(decoded_from_mask, np.array([[2]], dtype=np.int8))
    np.testing.assert_array_equal(
        decoded_from_masked_array,
        np.array([[2]], dtype=np.int8),
    )


def test_decode_snow_ice_mask_handles_masked_and_nonfinite_input():
    raw = np.ma.array(
        [[0b11011111, 0b11111111, np.nan]],
        mask=[[False, True, False]],
    )

    decoded = MOD06Granule._decode_snow_ice_first_byte(raw)

    np.testing.assert_array_equal(decoded, np.array([[0, 2, 2]], dtype=np.int8))


def test_snow_ice_output_validation_does_not_modify_other_modis_fields():
    processor = MainProcessor.__new__(MainProcessor)
    processor.log = logging.getLogger("test_snow_ice_mask")
    shape = (2, 3)
    ctp = np.arange(6, dtype=np.float32).reshape(shape)
    modis_vars = {
        "ctp": ctp.copy(),
        "snow_ice_mask": np.array([[0, 1, 9], [2, 1, 0]], dtype=np.int8),
    }

    result = processor._modis_variables_with_snow_ice_mask(modis_vars, shape)

    np.testing.assert_array_equal(result["ctp"], ctp)
    np.testing.assert_array_equal(
        result["snow_ice_mask"],
        np.array([[0, 1, 2], [2, 1, 0]], dtype=np.int8),
    )
    np.testing.assert_array_equal(modis_vars["snow_ice_mask"], np.array([[0, 1, 9], [2, 1, 0]]))

    missing = processor._modis_variables_with_snow_ice_mask(None, shape)
    np.testing.assert_array_equal(missing["snow_ice_mask"], np.full(shape, 2, dtype=np.int8))


def test_saver_adds_only_snow_ice_mask_to_existing_data_variables(tmp_path):
    shape = (2, 3)
    geo = {
        "lat": np.arange(6, dtype=np.float32).reshape(shape),
        "lon": np.arange(6, dtype=np.float32).reshape(shape) + 10.0,
        "vza": np.zeros(shape, dtype=np.float32),
        "elevation": np.arange(6, dtype=np.float32).reshape(shape) * 10.0,
        "land_sea_mask": np.array([[0, 1, 2], [3, 4, 5]], dtype=np.int16),
    }
    mm = {"cth": np.full(shape, 1000, dtype=np.int16)}
    existing_modis = {
        "ctp": np.full(shape, 500.0, dtype=np.float32),
        "cloud_mask_1km": np.array([[0, 1, 2], [3, 4, 5]], dtype=np.int8),
        "cmulti": np.zeros(shape, dtype=np.int8),
    }
    era5 = {"surface_pressure": np.full((1, 1), 1000.0, dtype=np.float32)}
    era5_time = datetime(2020, 1, 2, 3, 4, tzinfo=timezone.utc)
    snow_ice = np.array([[0, 1, 2], [1, 0, 1]], dtype=np.int8)
    before_path = tmp_path / "before.nc"
    after_path = tmp_path / "after.nc"

    XarraySaver(str(before_path)).save_mm(
        mm_variables=mm,
        mm_geo=geo,
        modis_variables=existing_modis,
        era5_variables=era5,
        era5_time=era5_time,
    )
    XarraySaver(str(after_path)).save_mm(
        mm_variables=mm,
        mm_geo=geo,
        modis_variables={**existing_modis, "snow_ice_mask": snow_ice},
        era5_variables=era5,
        era5_time=era5_time,
    )

    with xr.open_dataset(before_path, mask_and_scale=False) as before, xr.open_dataset(
        after_path, mask_and_scale=False
    ) as after:
        assert set(after.data_vars) == set(before.data_vars) | {"SnowIceMask"}
        for name in before.data_vars:
            np.testing.assert_array_equal(after[name].values, before[name].values)
            assert set(after[name].attrs) == set(before[name].attrs)
            for attr_name, before_value in before[name].attrs.items():
                np.testing.assert_equal(after[name].attrs[attr_name], before_value)

        snow = after["SnowIceMask"]
        np.testing.assert_array_equal(snow.values, snow_ice)
        np.testing.assert_array_equal(snow.attrs["flag_values"], np.array([0, 1, 2], dtype=np.int8))
        assert snow.attrs["flag_meanings"] == (
            "snow_or_ice_background no_snow_or_ice_background no_result_or_missing"
        )

        assert {
            "SnowIceMask",
            "SurfaceElevation",
            "LandSeaMask",
            "ERA5_Time",
            "MODIS_CloudMultilayerFlag",
        }.issubset(after.variables)
        assert {
            "Snow_Ice_Mask",
            "elevation",
            "land_sea_mask",
            "ERA5_time",
            "MODIS_CloudMultiLayerFlag",
        }.isdisjoint(after.variables)
        assert "ERA5_Time" in after["ERA5_SurfacePressure"].coords
        cloud_mask = after["MODIS_CloudMask"]
        assert cloud_mask.attrs["_FillValue"] == np.int8(5)
        assert cloud_mask.attrs["cloud_mask_decoder_version"] == (
            "modis_cloud_mask_byte0_bit0_gate_native_bits1_2_v2"
        )
