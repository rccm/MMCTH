import numpy as np

from src.data_processors.xarray_saver import XarraySaver
from src.data_readers.mod_read import MOD06Granule


def test_cloud_mask_decoder_gates_confidence_with_determination_bit():
    # Low three bits are [confidence bit 2, confidence bit 1, determined].
    raw = np.arange(8, dtype=np.uint8)[None, :]

    decoded = MOD06Granule._decode_cloud_mask_first_byte(raw)

    np.testing.assert_array_equal(
        decoded,
        np.array([[4, 0, 4, 1, 4, 2, 4, 3]], dtype=np.int8),
    )


def test_cloud_mask_decoder_preserves_signed_bytes_without_using_packed_fill_equality():
    raw_u8 = np.array(
        [[0b11100001, 0b11100011, 0b11100101, 0b11100111, 0xFF]],
        dtype=np.uint8,
    )
    raw_i8 = raw_u8.view(np.int8)

    decoded_signed = MOD06Granule._decode_cloud_mask_first_byte(
        raw_i8,
        fill_value=np.int8(-1),
    )
    decoded_unsigned = MOD06Granule._decode_cloud_mask_first_byte(
        raw_u8,
        fill_value=np.uint8(0),
    )

    np.testing.assert_array_equal(
        decoded_signed,
        np.array([[0, 1, 2, 3, 3]], dtype=np.int8),
    )
    np.testing.assert_array_equal(decoded_signed, decoded_unsigned)


def test_cloud_mask_decoder_does_not_treat_packed_fill_value_as_source_missing():
    raw = np.array([[0, 2, 4, 6]], dtype=np.uint8)

    decoded = MOD06Granule._decode_cloud_mask_first_byte(
        raw,
        fill_value=np.uint8(0),
    )

    np.testing.assert_array_equal(
        decoded,
        np.array([[4, 4, 4, 4]], dtype=np.int8),
    )


def test_cloud_mask_decoder_uses_class_five_only_for_explicit_source_missing():
    raw = np.array([[0b00000001]], dtype=np.uint8)

    decoded_from_mask = MOD06Granule._decode_cloud_mask_first_byte(
        raw,
        source_missing=np.array([[True]]),
    )
    decoded_from_masked_array = MOD06Granule._decode_cloud_mask_first_byte(
        np.ma.array(raw, mask=[[True]]),
    )

    np.testing.assert_array_equal(decoded_from_mask, np.array([[5]], dtype=np.int8))
    np.testing.assert_array_equal(
        decoded_from_masked_array,
        np.array([[5]], dtype=np.int8),
    )


def test_cloud_mask_metadata_matches_decoded_bits_and_values():
    attrs = XarraySaver.MODIS_SCHEMA["cloud_mask_1km"]

    np.testing.assert_array_equal(
        attrs["flag_values"],
        np.array([0, 1, 2, 3, 4], dtype=np.int8),
    )
    np.testing.assert_array_equal(
        attrs["valid_values"],
        np.array([0, 1, 2, 3, 4], dtype=np.int8),
    )
    assert attrs["flag_meanings"] == (
        "confident_cloudy probably_cloudy probably_clear confident_clear "
        "not_determined"
    )
    assert "bit 0" in attrs["description"]
    assert "bits 1-2" in attrs["description"]
    assert attrs["cloud_mask_decoder_version"] == (
        "modis_cloud_mask_byte0_bit0_gate_native_bits1_2_v2"
    )

    encoding = XarraySaver("unused.nc")._enc_from_dtype("MODIS_CloudMask", "int8")
    assert encoding["fill"] == np.int8(5)
