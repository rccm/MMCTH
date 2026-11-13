"""MM correction provenance, field validity, and production NetCDF regression checks."""

import copy
import logging
from types import SimpleNamespace

import numpy as np
import pytest
import xarray as xr

from src.data_processors import process_main
from src.data_processors.process_main import MainProcessor


CORRECTED_KEYS = {
    "cth": "cth_mmcorrected",
    "ctp": "ctp_mmcorrected",
    "emissivity": "emissivity_mmcorrected",
}
CORRECTED_NAMES = {
    "MODIS_CloudTopHeight_MMCorrected",
    "MODIS_CloudTopPressure_MMCorrected",
    "MODIS_CloudEffectiveEmissivity_MMCorrected",
}


def inputs(shape=(1, 1)):
    original = {
        "cth": np.full(shape, 6000.0, dtype=np.float32),
        "ctp": np.full(shape, 450.0, dtype=np.float32),
        "emissivity": np.full(shape, 0.65, dtype=np.float32),
    }
    mm = {
        "cth": np.full(shape, 10000.0, dtype=np.float32),
        "ctp": np.full(shape, 250.0, dtype=np.float32),
        "emissivity": np.full(shape, 0.25, dtype=np.float32),
        "cflag": np.full(shape, 8, dtype=np.int8),
        "failure_reason": np.zeros(shape, dtype=np.int8),
    }
    return original, mm


def merge(original, mm, shape=(1, 1)):
    processor = MainProcessor.__new__(MainProcessor)
    processor._add_modis_mm_corrected_variables(original, mm, shape)
    return original


@pytest.mark.parametrize("flag", [8, 9])
def test_accepted_correction_replaces_each_field_without_aliasing(flag):
    original, mm = inputs()
    mm["cflag"][:] = flag
    originals_before = copy.deepcopy(original)
    mm_before = copy.deepcopy(mm)

    merge(original, mm)

    for key, corrected in CORRECTED_KEYS.items():
        np.testing.assert_array_equal(original[corrected], mm[key])
        np.testing.assert_array_equal(original[key], originals_before[key])
        np.testing.assert_array_equal(mm[key], mm_before[key])
        assert not np.shares_memory(original[corrected], original[key])
        assert not np.shares_memory(original[corrected], mm[key])
        original[corrected][:] = -123.0
        np.testing.assert_array_equal(original[key], originals_before[key])
        np.testing.assert_array_equal(mm[key], mm_before[key])


@pytest.mark.parametrize("reason", [70, 71])
def test_emissivity_failure_keeps_accepted_cloud_top_pair(reason):
    original, mm = inputs()
    mm["failure_reason"][:] = reason
    merge(original, mm)

    np.testing.assert_array_equal(original["cth_mmcorrected"], [[10000.0]])
    np.testing.assert_array_equal(original["ctp_mmcorrected"], [[250.0]])
    np.testing.assert_array_equal(original["emissivity_mmcorrected"], original["emissivity"])


@pytest.mark.parametrize("flag", range(8))
def test_finite_modis_or_misr_fallbacks_are_not_corrections(flag):
    original, mm = inputs()
    mm["cflag"][:] = flag
    # A finite fallback and even an OK diagnostic do not establish MM selection.
    merge(original, mm)
    for key, corrected in CORRECTED_KEYS.items():
        np.testing.assert_array_equal(original[corrected], original[key])
        assert not np.shares_memory(original[corrected], original[key])


@pytest.mark.parametrize("reason", [10, 20, 40, 60, 61, 62])
def test_success_flag_with_inconsistent_failure_reason_keeps_modis(reason):
    original, mm = inputs()
    mm["failure_reason"][:] = reason
    merge(original, mm)
    for key, corrected in CORRECTED_KEYS.items():
        np.testing.assert_array_equal(original[corrected], original[key])


@pytest.mark.parametrize(
    "key,value",
    [
        ("cth", np.nan), ("cth", np.inf), ("cth", -9999.0),
        ("cth", -500.1), ("cth", 20000.1),
        ("ctp", np.nan), ("ctp", -np.inf), ("ctp", -9999.0),
        ("ctp", 0.0), ("ctp", 1100.1),
    ],
)
def test_bad_corrected_cloud_top_value_preserves_both_original_fields(key, value):
    original, mm = inputs()
    mm[key][:] = value
    merge(original, mm)
    np.testing.assert_array_equal(original["cth_mmcorrected"], original["cth"])
    np.testing.assert_array_equal(original["ctp_mmcorrected"], original["ctp"])
    # Emissivity has a separate numerical-validity outcome.
    np.testing.assert_array_equal(original["emissivity_mmcorrected"], mm["emissivity"])


@pytest.mark.parametrize("height,pressure", [(-500.0, 1100.0), (20000.0, 0.01)])
def test_corrected_cloud_top_range_endpoints_are_valid(height, pressure):
    original, mm = inputs()
    mm["cth"][:] = height
    mm["ctp"][:] = pressure
    merge(original, mm)
    np.testing.assert_array_equal(original["cth_mmcorrected"], mm["cth"])
    np.testing.assert_array_equal(original["ctp_mmcorrected"], mm["ctp"])


@pytest.mark.parametrize(
    "value,accepted",
    [(np.nan, False), (np.inf, False), (-9999.0, False), (0.0, False),
     (0.009, False), (np.nextafter(np.float32(0.01), np.float32(0.0)), False),
     (np.float32(0.01), True), (1.0, True), (1.001, False)],
)
def test_emissivity_uses_its_own_validity_range(value, accepted):
    original, mm = inputs()
    mm["emissivity"][:] = value
    merge(original, mm)
    expected = mm["emissivity"] if accepted else original["emissivity"]
    np.testing.assert_array_equal(original["emissivity_mmcorrected"], expected)
    np.testing.assert_array_equal(original["cth_mmcorrected"], mm["cth"])
    np.testing.assert_array_equal(original["ctp_mmcorrected"], mm["ctp"])


@pytest.mark.parametrize(
    "key,value",
    [("cth", np.nan), ("cth", np.inf), ("cth", -9999.0),
     ("cth", -1.0), ("cth", 20001.0), ("ctp", np.nan),
     ("ctp", -9999.0), ("ctp", 0.0), ("ctp", 1101.0)],
)
def test_missing_or_invalid_original_cloud_top_does_not_get_backfilled(key, value):
    original, mm = inputs()
    original[key][:] = value
    merge(original, mm)
    for cloud_top_key in ("cth", "ctp"):
        np.testing.assert_array_equal(
            original[CORRECTED_KEYS[cloud_top_key]], original[cloud_top_key]
        )
    np.testing.assert_array_equal(original["emissivity_mmcorrected"], mm["emissivity"])


@pytest.mark.parametrize("value", [np.nan, np.inf, -9999.0, -0.01, 1.01])
def test_missing_or_invalid_original_emissivity_does_not_get_backfilled(value):
    original, mm = inputs()
    original["emissivity"][:] = value
    merge(original, mm)
    np.testing.assert_array_equal(original["emissivity_mmcorrected"], original["emissivity"])
    np.testing.assert_array_equal(original["cth_mmcorrected"], mm["cth"])
    np.testing.assert_array_equal(original["ctp_mmcorrected"], mm["ctp"])


@pytest.mark.parametrize("value", [0.0, 1.0])
def test_original_emissivity_endpoints_can_be_corrected(value):
    original, mm = inputs()
    original["emissivity"][:] = value
    merge(original, mm)
    np.testing.assert_array_equal(original["emissivity_mmcorrected"], mm["emissivity"])


@pytest.mark.parametrize("missing", ["cflag", "failure_reason"])
def test_missing_retrieval_provenance_keeps_original_values(missing):
    original, mm = inputs()
    del mm[missing]
    merge(original, mm)
    for key, corrected in CORRECTED_KEYS.items():
        np.testing.assert_array_equal(original[corrected], original[key])


@pytest.mark.parametrize("missing", ["cth", "ctp", "emissivity"])
def test_absent_corrected_fields_keep_the_corresponding_original_values(missing):
    original, mm = inputs()
    del mm[missing]
    merge(original, mm)
    np.testing.assert_array_equal(original[CORRECTED_KEYS[missing]], original[missing])
    if missing in ("cth", "ctp"):
        other = "ctp" if missing == "cth" else "cth"
        np.testing.assert_array_equal(original[CORRECTED_KEYS[other]], original[other])
        np.testing.assert_array_equal(original["emissivity_mmcorrected"], mm["emissivity"])


def test_absent_original_fields_stay_missing_on_the_output_grid():
    shape = (2, 3)
    _, mm = inputs(shape)
    original = {}
    merge(original, mm, shape)
    for corrected in CORRECTED_KEYS.values():
        assert original[corrected].shape == shape
        assert np.isnan(original[corrected]).all()


@pytest.mark.parametrize("missing", ["cth", "ctp", "emissivity"])
def test_partly_missing_original_fields_do_not_get_backfilled(missing):
    original, mm = inputs()
    del original[missing]
    merge(original, mm)
    assert np.isnan(original[CORRECTED_KEYS[missing]]).all()
    if missing in ("cth", "ctp"):
        other = "ctp" if missing == "cth" else "cth"
        np.testing.assert_array_equal(original[CORRECTED_KEYS[other]], original[other])


@pytest.mark.parametrize("masked", ["cth", "ctp", "emissivity", "cflag", "failure_reason"])
def test_masked_corrected_values_and_provenance_are_missing(masked):
    original, mm = inputs()
    mm[masked] = np.ma.array(mm[masked], mask=np.ones((1, 1), dtype=bool))
    merge(original, mm)
    if masked in ("cth", "ctp", "cflag", "failure_reason"):
        for key in ("cth", "ctp"):
            np.testing.assert_array_equal(original[CORRECTED_KEYS[key]], original[key])
    if masked in ("emissivity", "cflag", "failure_reason"):
        np.testing.assert_array_equal(original["emissivity_mmcorrected"], original["emissivity"])


@pytest.mark.parametrize("masked", ["cth", "ctp", "emissivity"])
def test_masked_original_values_remain_missing(masked):
    original, mm = inputs()
    original[masked] = np.ma.array(original[masked], mask=np.ones((1, 1), dtype=bool))
    merge(original, mm)
    assert np.isnan(original[CORRECTED_KEYS[masked]]).all()
    if masked in ("cth", "ctp"):
        other = "ctp" if masked == "cth" else "cth"
        np.testing.assert_array_equal(original[CORRECTED_KEYS[other]], original[other])


def production_inputs():
    shape = (20, 20)
    modis, mm = inputs(shape)
    # Cover success, field-specific failure, fallbacks, and missing original data.
    mm["cflag"][0, :4] = [8, 9, 3, 7]
    mm["failure_reason"][0, :4] = [0, 70, 10, 40]
    modis["cth"][1, 0] = np.nan
    modis["ctp"][1, 1] = np.nan
    modis["emissivity"][1, 2] = np.nan
    mm["emissivity"][0, 1] = -9999.0
    modis.update({
        "opt": np.full(shape, 4.0, dtype=np.float32),
        "ctm": np.ones(shape, dtype=np.int8),
        "cphase": np.full(shape, 2, dtype=np.int8),
        "cloud_mask_1km": np.zeros(shape, dtype=np.int8),
        "snow_ice_mask": np.ones(shape, dtype=np.int8),
    })
    mm["opt"] = np.full(shape, 0.5, dtype=np.float32)
    geo = {
        "lat": np.linspace(20, 21, 400, dtype=np.float32).reshape(shape),
        "lon": np.linspace(-100, -99, 400, dtype=np.float32).reshape(shape),
        "vza": np.zeros(shape, dtype=np.float32),
        "sza": np.full(shape, 45, dtype=np.float32),
        "elevation": np.zeros(shape, dtype=np.float32),
        "land_sea_mask": np.ones(shape, dtype=np.int16),
        "time1d": np.arange(20, dtype=np.float64),
        "time2d": np.arange(400, dtype=np.float64).reshape(shape),
        "orig_row_index": np.arange(20, dtype=np.int32),
        "orig_col_index": np.arange(20, dtype=np.int32),
        "orig_nrows": 20,
        "orig_ncols": 20,
    }
    misr = {
        "misrcth": np.full(shape, 2000.0, dtype=np.float32),
        "misrctp": np.full(shape, 800.0, dtype=np.float32),
        "misrcth_qa": np.full(shape, 90, dtype=np.int16),
    }
    era5 = {"surface_pressure": np.full(shape, 100000.0, dtype=np.float32)}
    return modis, mm, geo, misr, era5


def production_processor():
    processor = MainProcessor.__new__(MainProcessor)
    processor.log = logging.getLogger(__name__)
    processor.mm_processor = SimpleNamespace(orbit="O123456")
    processor.orbit = "O123456"
    processor.has_mod06 = True
    processor.has_tc = True
    processor.input_files = ["MOD06.hdf", "MOD03.hdf", "MISR.hdf"]
    processor.met_date = [2020, 1, 2, 0]
    processor.timestamp = "0304"
    return processor


def test_save_pixels_adds_only_three_variables_and_preserves_existing_product(tmp_path, monkeypatch):
    modis, mm, geo, misr, era5 = production_inputs()
    processor = production_processor()
    monkeypatch.setattr(process_main, "create_json_sidecar", lambda *args, **kwargs: None)
    before_path = tmp_path / "before.nc"
    after_path = tmp_path / "after.nc"
    with monkeypatch.context() as baseline_patch:
        baseline_patch.setattr(
            MainProcessor, "_add_modis_mm_corrected_variables", lambda *args, **kwargs: None
        )
        processor.save_pixels(mm, modis, geo, misr, era5, str(before_path))
    processor.save_pixels(mm, modis, geo, misr, era5, str(after_path))

    with xr.open_dataset(before_path, decode_cf=False) as before, xr.open_dataset(
        after_path, decode_cf=False
    ) as after:
        assert set(after.variables) == set(before.variables) | CORRECTED_NAMES
        for name in before.variables:
            xr.testing.assert_identical(after[name], before[name])
            assert after[name].dtype == before[name].dtype
        for key in before.attrs:
            if key not in {"history", "date_created"}:
                np.testing.assert_equal(after.attrs[key], before.attrs[key])
        assert set(after.attrs) == set(before.attrs)

        expected = copy.deepcopy(modis)
        merge(expected, mm, geo["lat"].shape)
        for key, name, dtype, units in (
            ("cth_mmcorrected", "MODIS_CloudTopHeight_MMCorrected", "int16", "m"),
            ("ctp_mmcorrected", "MODIS_CloudTopPressure_MMCorrected", "float32", "hPa"),
            ("emissivity_mmcorrected", "MODIS_CloudEffectiveEmissivity_MMCorrected", "float32", "1"),
        ):
            variable = after[name]
            assert variable.dims == ("cell_along_swath_1km", "cell_cross_swath_1km")
            assert variable.dtype == np.dtype(dtype)
            assert variable.attrs["units"] == units
            assert variable.attrs["_FillValue"] == -9999
            assert variable.attrs["coordinates"] == "Latitude Longitude pixel_time"
            expected_values = np.where(np.isnan(expected[key]), -9999, expected[key]).astype(dtype)
            np.testing.assert_array_equal(variable.values, expected_values)


def test_save_pixels_without_modis_writes_missing_corrected_fields(tmp_path, monkeypatch):
    _, mm, geo, misr, era5 = production_inputs()
    processor = production_processor()
    processor.has_mod06 = False
    monkeypatch.setattr(process_main, "create_json_sidecar", lambda *args, **kwargs: None)
    path = tmp_path / "misr_only.nc"
    processor.save_pixels(mm, None, geo, misr, era5, str(path))
    with xr.open_dataset(path, decode_cf=False) as dataset:
        for name in CORRECTED_NAMES:
            assert dataset[name].shape == geo["lat"].shape
            assert (dataset[name].values == dataset[name].attrs["_FillValue"]).all()
