"""Regression tests for the shared stacking-grid convention.

Two conventions used to be spelled out separately in sim_utils, stack_single_jk
and combine_jackknife.  They agreed at 100 pixels and disagreed at 101: the
cell-centred formula gives a spacing of 100/101 = 0.990 on an odd grid, so a
101-pixel single stack came out 1% smaller than the 101-pixel pair stack it was
subtracted from.  These tests pin the single definition everything now shares.
"""

from __future__ import annotations

import unittest

import numpy as np

from geometry import (assert_unit_pixel_scale, stack_axis,
                      symmetrized_radial_interpolator)

BOX = 100.0


class StackAxisTests(unittest.TestCase):
    def test_spacing_is_unity_at_both_parities(self):
        for n in (100, 101):
            with self.subTest(grid=n):
                self.assertAlmostEqual(stack_axis(BOX, n)[1], 1.0, places=12)

    def test_parity_decides_whether_a_pixel_sits_on_zero(self):
        even, _ = stack_axis(BOX, 100)
        odd, _ = stack_axis(BOX, 101)
        self.assertFalse(np.any(even == 0.0))
        self.assertTrue(np.any(odd == 0.0))

    def test_endpoints(self):
        even, _ = stack_axis(BOX, 100)
        odd, _ = stack_axis(BOX, 101)
        self.assertAlmostEqual(even[0], -49.5)
        self.assertAlmostEqual(even[-1], +49.5)
        self.assertAlmostEqual(odd[0], -50.0)
        self.assertAlmostEqual(odd[-1], +50.0)

    def test_axis_is_symmetric(self):
        for n in (100, 101):
            axis, _ = stack_axis(BOX, n)
            np.testing.assert_allclose(axis, -axis[::-1], atol=1e-12)

    def test_single_and_pair_axes_agree_at_101(self):
        """The bug this module exists to prevent."""
        from sim_utils import pair_stack_offsets, single_stack_offsets
        for n in (100, 101):
            with self.subTest(grid=n):
                s, _, _ = single_stack_offsets(BOX, n)
                p, _, _ = pair_stack_offsets(BOX, n)
                np.testing.assert_allclose(s, p, atol=1e-12)

    def test_every_consumer_shares_the_convention(self):
        from combine_jackknife import grid_radius_hmpc
        from deconvolve_pair_profile import axis_for
        from sim_utils import single_stack_offsets
        for n in (100, 101):
            with self.subTest(grid=n):
                axis, _ = stack_axis(BOX, n)
                np.testing.assert_allclose(axis_for(n), axis, atol=1e-12)
                np.testing.assert_allclose(single_stack_offsets(BOX, n)[0], axis,
                                           atol=1e-12)
                gx, gy = np.meshgrid(axis, axis)
                np.testing.assert_allclose(
                    grid_radius_hmpc(n).reshape(n, n), np.hypot(gx, gy),
                    atol=1e-12)

    def test_unit_pixel_scale_guard(self):
        self.assertAlmostEqual(assert_unit_pixel_scale(100), 1.0)
        self.assertAlmostEqual(assert_unit_pixel_scale(101), 1.0)
        with self.assertRaises(ValueError):
            assert_unit_pixel_scale(64)          # 100/64 = 1.5625 h^-1 Mpc

    def test_radial_interpolator_rejects_non_unit_pixels(self):
        """It measures radius in pixels and is evaluated in h^-1 Mpc."""
        with self.assertRaises(ValueError):
            symmetrized_radial_interpolator(np.zeros((64, 64)), validate=False)


class StackSingleJkGridTests(unittest.TestCase):
    """The stacker must follow the runtime grid, not the imported default."""

    def tearDown(self):
        import stack_single_jk as sjk
        sjk.set_grid(100)

    def test_set_grid_updates_offsets_and_width(self):
        import stack_single_jk as sjk
        for n in (100, 101):
            with self.subTest(grid=n):
                sjk.set_grid(n)
                self.assertEqual(sjk.GRID, n)
                self.assertEqual(sjk.OFF_X.size, n * n)
                axis, _ = stack_axis(BOX, n)
                np.testing.assert_allclose(sjk.OFFSETS, axis, atol=1e-12)

    def test_accumulator_width_follows_the_grid(self):
        import stack_single_jk as sjk
        for n in (100, 101):
            with self.subTest(grid=n):
                sjk.set_grid(n)
                state = sjk.init_reducer(n_regions=3, n_chunks=1)
                self.assertEqual(state["sum_wk"].shape, (3, n * n))
                self.assertEqual(state["sum_w"].shape, (3, n * n))

    def test_mean_map_reshapes_to_the_grid(self):
        import stack_single_jk as sjk
        for n in (100, 101):
            with self.subTest(grid=n):
                sjk.set_grid(n)
                out = sjk.mean_map(np.ones(n * n), np.ones(n * n))
                self.assertEqual(out.shape, (n, n))

    def test_set_grid_rejects_non_unit_pixels(self):
        import stack_single_jk as sjk
        with self.assertRaises(ValueError):
            sjk.set_grid(64)


if __name__ == "__main__":
    unittest.main()
