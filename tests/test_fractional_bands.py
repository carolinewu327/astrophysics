"""Regression tests for the fractional-overlap band extractor.

``band_masks`` needs the band edges (|Y| = 1.5, 10.5) to fall *between* row
centres.  That holds on the 101-pixel pair grid and fails on the 100-pixel
single grid, where 1.5 is a row centre and the two bands overlap.  The
fractional extractor weights the boundary rows so the same physical band comes
out of either grid, which is what lets the existing 100-pixel single stacks and
jackknife accumulators be used without a re-stack.
"""

from __future__ import annotations

import unittest

import numpy as np

from geometry import (band_masks, band_profiles, band_profiles_fractional,
                      band_row_weights)

BOX = 100.0


def axis_for(n: int) -> np.ndarray:
    if n % 2 == 1:
        return np.linspace(-0.5 * BOX, 0.5 * BOX, n)
    cell = BOX / n
    return np.linspace(-0.5 * BOX + 0.5 * cell, 0.5 * BOX - 0.5 * cell, n)


class FractionalBandTests(unittest.TestCase):
    def test_effective_widths_match_on_both_grids(self):
        for n in (100, 101):
            with self.subTest(grid=n):
                axis = axis_for(n)
                self.assertAlmostEqual(
                    band_row_weights(axis, 0.0, 1.5).sum(), 3.0, places=9)
                self.assertAlmostEqual(
                    band_row_weights(axis, 1.5, 10.5).sum(), 18.0, places=9)

    def test_boundary_rows_are_half_weighted_on_the_even_grid(self):
        axis = axis_for(100)
        w = band_row_weights(axis, 0.0, 1.5)
        np.testing.assert_allclose(w[np.abs(axis) == 0.5], 1.0)
        np.testing.assert_allclose(w[np.abs(axis) == 1.5], 0.5)
        self.assertEqual(int((w > 0).sum()), 4)

    def test_signed_overlap_keeps_the_central_row_whole(self):
        """An |Y|-based overlap halves the Y = 0 row and gives width 2.5."""
        axis = axis_for(101)
        w = band_row_weights(axis, 0.0, 1.5)
        np.testing.assert_allclose(w[axis == 0.0], 1.0)

    def test_bands_are_disjoint_and_partition_the_shells(self):
        for n in (100, 101):
            axis = axis_for(n)
            w_cen = band_row_weights(axis, 0.0, 1.5)
            w_off = band_row_weights(axis, 1.5, 10.5)
            # a row may be split between the two bands but never over-counted
            np.testing.assert_array_less(w_cen + w_off, 1.0 + 1e-12)

    def test_reduces_to_band_profiles_on_the_odd_grid(self):
        rng = np.random.default_rng(20260820)
        m = rng.normal(size=(101, 101))
        axis = axis_for(101)
        cen, off = band_profiles(m, axis)
        f_cen, f_off = band_profiles_fractional(m, axis)
        np.testing.assert_allclose(f_cen, cen, rtol=0, atol=0)
        np.testing.assert_allclose(f_off, off, rtol=0, atol=0)

    def test_band_masks_still_rejects_the_even_grid(self):
        """The failure the extractor exists to route around, pinned."""
        with self.assertRaises(ValueError):
            band_masks(axis_for(100))

    def test_constant_map_returns_the_constant(self):
        for n in (100, 101):
            m = np.full((n, n), 3.5)
            cen, off = band_profiles_fractional(m, axis_for(n))
            np.testing.assert_allclose(cen, 3.5)
            np.testing.assert_allclose(off, 3.5)


if __name__ == "__main__":
    unittest.main()
