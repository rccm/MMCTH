import numpy as np
import xarray as xr

from src.data_processors.process_main import MainProcessor
from src.data_processors.xarray_saver import XarraySaver


def test_tiny_overlap_keeps_matching_zero_sized_block_dimensions():
    processor = MainProcessor.__new__(MainProcessor)
    lat = np.zeros((5, 170), dtype=np.float32)
    lon = np.zeros_like(lat)

    averaged = processor.block_average_dict({"surface_pressure": lat})
    block_lat, block_lon = processor.calculate_block_lat_lon(lat, lon)

    assert averaged["surface_pressure"].shape == (0, 17)
    assert block_lat.shape == averaged["surface_pressure"].shape
    assert block_lon.shape == averaged["surface_pressure"].shape


def test_tiny_overlap_keeps_matching_shape_when_cross_track_is_too_small():
    processor = MainProcessor.__new__(MainProcessor)
    lat = np.zeros((200, 5), dtype=np.float32)
    lon = np.zeros_like(lat)

    averaged = processor.block_average_dict({"surface_pressure": lat})
    block_lat, block_lon = processor.calculate_block_lat_lon(lat, lon)

    assert averaged["surface_pressure"].shape == (20, 0)
    assert block_lat.shape == averaged["surface_pressure"].shape
    assert block_lon.shape == averaged["surface_pressure"].shape


def test_tiny_overlap_zero_sized_era5_grid_can_be_saved(tmp_path):
    swath_shape = (5, 170)
    era5_shape = (0, 17)
    output = tmp_path / "tiny_overlap.nc"

    XarraySaver(str(output)).save_mm(
        mm_variables={"cth": np.zeros(swath_shape, dtype=np.int16)},
        mm_geo={
            "lat": np.zeros(swath_shape, dtype=np.float32),
            "lon": np.zeros(swath_shape, dtype=np.float32),
            "vza": np.zeros(swath_shape, dtype=np.float32),
        },
        era5_variables={
            "surface_pressure": np.empty(era5_shape, dtype=np.float32),
            "ERA5Latitude": np.empty(era5_shape, dtype=np.float32),
            "ERA5Longitude": np.empty(era5_shape, dtype=np.float32),
        },
    )

    with xr.open_dataset(output) as dataset:
        assert dataset.sizes["era5_cell_along_swath_10km"] == 0
        assert dataset.sizes["era5_cell_cross_swath_10km"] == 17
