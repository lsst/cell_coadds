# This file is part of cell_coadds.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import unittest

import numpy as np

import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.geom as geom
import lsst.utils.tests
from lsst.cell_coadds import CoaddApCorrMapStacker


class TestCoaddApCorrMapStacker(lsst.utils.tests.TestCase):
    """Test the CoaddApCorrMap stacker."""

    @classmethod
    def setUpClass(cls):  # noqa: D102
        cls.bbox = geom.Box2I(
            minimum=geom.Point2I(-10, -10),
            maximum=geom.Point2I(10, 10),
        )
        cls.evaluation_point = geom.Point2D(0, 0)
        cls.visit_count = 5
        cls.ap_corr_map_list = [afwImage.ApCorrMap() for _ in range(cls.visit_count)]
        # Set them up to the spatially constant BoundedFields.
        for idx, ap_corr_map in enumerate(cls.ap_corr_map_list, start=2):
            ap_corr_map.set(
                "base_PsfFlux_instFlux", afwMath.ChebyshevBoundedField(cls.bbox, np.array([[1 / idx]]))
            )
            ap_corr_map.set(
                "base_PsfFlux_instFluxErr", afwMath.ChebyshevBoundedField(cls.bbox, np.array([[0.05 / idx]]))
            )
            ap_corr_map.set(
                "base_GaussianFlux_instFlux", afwMath.ChebyshevBoundedField(cls.bbox, np.array([[1.5 / idx]]))
            )
            ap_corr_map.set(
                "base_GaussianFlux_instFluxErr",
                afwMath.ChebyshevBoundedField(cls.bbox, np.array([[0.6 / idx]])),
            )

    def test_add_direct_ap_corr(self):
        """Test that the values are what we expect when adding directly."""
        stacker = CoaddApCorrMapStacker(
            evaluation_point=self.evaluation_point,
            do_coadd_inverse_ap_corr=False,
        )
        for idx in range(self.visit_count):
            stacker.add(
                self.ap_corr_map_list[idx],
                weight=1.0,
            )

        final_ap_corr_map = stacker.final_ap_corr_map

        self.assertFloatsAlmostEqual(
            final_ap_corr_map["base_PsfFlux_instFlux"],
            0.29,  # mean of 1/2, 1/3, 1/4, 1/5, 1/6
            atol=1e-8,
        )
        self.assertFloatsAlmostEqual(
            final_ap_corr_map["base_PsfFlux_instFluxErr"],
            0.05 * 0.14019827229875395,
            atol=1e-8,
        )
        self.assertFloatsAlmostEqual(
            final_ap_corr_map["base_GaussianFlux_instFlux"],
            1.5 * 0.29,
            atol=1e-8,
        )
        self.assertFloatsAlmostEqual(
            final_ap_corr_map["base_GaussianFlux_instFluxErr"],
            0.6 * 0.14019827229875395,
            atol=1e-8,
        )

    def test_add_inverse_ap_corr(self):
        """Test that the values are what we expect when adding inverse."""
        stacker = CoaddApCorrMapStacker(
            evaluation_point=self.evaluation_point,
            do_coadd_inverse_ap_corr=True,
        )
        for idx in range(self.visit_count):
            stacker.add(
                self.ap_corr_map_list[idx],
                weight=1.0,
            )

        final_ap_corr_map = stacker.final_ap_corr_map

        self.assertFloatsAlmostEqual(
            final_ap_corr_map["base_PsfFlux_instFlux"],
            0.25,
            atol=1e-8,
        )

        self.assertFloatsAlmostEqual(
            final_ap_corr_map["base_GaussianFlux_instFlux"],
            1.5 * 0.25,
            atol=1e-8,
        )

        self.assertFloatsAlmostEqual(
            final_ap_corr_map["base_PsfFlux_instFluxErr"],
            0.05 * (0.25**2) * np.sqrt(90) / 5,
            atol=1e-8,
        )

        self.assertFloatsAlmostEqual(
            final_ap_corr_map["base_GaussianFlux_instFluxErr"],
            0.6 / (1.5**2) * ((1.5 * 0.25) ** 2) * np.sqrt(90) / 5,
            atol=1e-8,
        )


class TestCoaddApCorrMapMemoryTestCase(lsst.utils.tests.MemoryTestCase):
    """Test the CoaddApCorrMapStacker for memory leaks."""


def setup_module(module):  # noqa: D103
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
