"""Regression tests for even/odd single-stack centering conventions."""

from __future__ import annotations

import unittest

import numpy as np

from geometry import (
    symmetrize_map as boss_symmetrize_map,
    two_halo_template as boss_two_halo_template,
    validate_radially_symmetrized_map as boss_validate_radially_symmetrized_map,
)
from sim_utils import (
    make_two_halo_template,
    radial_symmetrize_map,
    single_stack_offsets,
    symmetrized_radial_interpolator,
    validate_radially_symmetrized_map,
)


class RadialCenteringTests(unittest.TestCase):
    def test_even_grid_origin_is_between_four_central_pixels(self):
        axis, _, _ = single_stack_offsets(box_size_hmpc=100.0, grid_size=100)
        self.assertEqual(axis[49], -0.5)
        self.assertEqual(axis[50], 0.5)

    def test_even_grid_symmetrization_uses_physical_center(self):
        rng = np.random.default_rng(17)
        original = rng.normal(size=(100, 100))
        result = radial_symmetrize_map(original)

        np.testing.assert_array_equal(result, result[::-1, :])
        np.testing.assert_array_equal(result, result[:, ::-1])
        np.testing.assert_allclose(result, boss_symmetrize_map(original), rtol=0, atol=0)
        validate_radially_symmetrized_map(result)

    def test_odd_grid_center_remains_central_pixel(self):
        original = np.zeros((101, 101), dtype=float)
        original[50, 50] = 1.0
        result = radial_symmetrize_map(original)

        self.assertEqual(np.unravel_index(np.argmax(result), result.shape), (50, 50))
        np.testing.assert_array_equal(result, result[::-1, :])
        np.testing.assert_array_equal(result, result[:, ::-1])

    def test_old_even_grid_product_is_rejected(self):
        yy, xx = np.indices((100, 100))
        old_product = np.exp(-np.hypot(xx - 50.0, yy - 50.0) / 5.0)

        with self.assertRaisesRegex(ValueError, "old N//2 centering"):
            validate_radially_symmetrized_map(old_product)
        with self.assertRaisesRegex(ValueError, "old N//2 centering"):
            symmetrized_radial_interpolator(old_product)
        axis = np.arange(-50.0, 51.0)
        gx, gy = np.meshgrid(axis, axis)
        with self.assertRaisesRegex(ValueError, "N//2 centering"):
            boss_two_halo_template(old_product, gx, gy, 10.0)

    def test_validators_agree_across_the_two_trees(self):
        """lib/geometry.py and sim_utils.py carry duplicate validators.

        They cannot be merged: the simulation tree runs with
        ``PYTHONPATH=analysis/sim`` and cannot import ``lib``.  So assert they
        behave identically instead -- this is the guard that keeps the two
        copies from drifting, which is how the centering conventions diverged in
        the first place.
        """
        yy, xx = np.indices((100, 100))
        cases = {
            "good": radial_symmetrize_map(np.exp(-np.hypot(xx - 49.5, yy - 49.5) / 5.0)),
            "mis_centered": np.exp(-np.hypot(xx - 50.0, yy - 50.0) / 5.0),
            "raw": np.random.default_rng(3).normal(size=(100, 100)),
            "non_square": np.zeros((100, 101)),
        }
        for name, arr in cases.items():
            with self.subTest(case=name):
                errs = []
                for fn in (validate_radially_symmetrized_map,
                           boss_validate_radially_symmetrized_map):
                    try:
                        fn(arr)
                        errs.append(None)
                    except ValueError as exc:
                        errs.append(str(exc))
                self.assertEqual(errs[0], errs[1],
                                 f"validators disagree on the {name} case")

    def test_validate_false_allows_a_deliberately_raw_map(self):
        """--no-symmetrize-single hands in a raw map on purpose.

        The symmetry check exists to catch archived mis-centered products, so it
        must not reject a caller that opted out.  Guarding this because adding
        the check silently broke that flag.
        """
        raw = np.random.default_rng(5).normal(size=(100, 100))
        axis = np.arange(-50.0, 51.0)
        gx, gy = np.meshgrid(axis, axis)

        with self.assertRaises(ValueError):
            boss_two_halo_template(raw, gx, gy, 10.0)
        for template in (boss_two_halo_template(raw, gx, gy, 10.0, validate=False),
                         make_two_halo_template(raw, gx, gy, 10.0, validate=False)):
            self.assertEqual(template.shape, (101, 101))
            self.assertTrue(np.isfinite(template).all())

    def test_two_halo_template_is_reflection_symmetric(self):
        yy, xx = np.indices((100, 100))
        center = 49.5
        single = radial_symmetrize_map(
            np.exp(-np.hypot(xx - center, yy - center) / 5.0)
        )
        axis = np.arange(-50.0, 51.0)
        gx, gy = np.meshgrid(axis, axis)
        template = make_two_halo_template(single, gx, gy, 10.0)

        np.testing.assert_allclose(template, template[::-1, :], rtol=0, atol=1e-15)
        np.testing.assert_allclose(template, template[:, ::-1], rtol=0, atol=1e-15)


if __name__ == "__main__":
    unittest.main()
