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

    @lsst.utils.tests.methodParameters(
        do_coadd_inverse_ap_corr=[True, False],
    )
    def test_weight_normalization(self, do_coadd_inverse_ap_corr):
        """Test that the coadd aperture corrections stay the same when the
        relative weights are the same.
        """
        stacker = CoaddApCorrMapStacker(
            evaluation_point=self.evaluation_point,
            do_coadd_inverse_ap_corr=do_coadd_inverse_ap_corr,
        )

        normalized_weights = np.random.rand(self.visit_count)
        normalized_weights /= np.sum(normalized_weights)

        for idx in range(self.visit_count):
            stacker.add(
                self.ap_corr_map_list[idx],
                weight=normalized_weights[idx],
            )

        reference_ap_corr_map = stacker.final_ap_corr_map

        for scale in (
            1.0,
            0.5,
            2.0,
        ):
            stacker = CoaddApCorrMapStacker(
                evaluation_point=self.evaluation_point,
                do_coadd_inverse_ap_corr=do_coadd_inverse_ap_corr,
            )
            weights = scale * normalized_weights
            for idx in range(self.visit_count):
                stacker.add(
                    self.ap_corr_map_list[idx],
                    weight=weights[idx],
                )
            final_ap_corr_map = stacker.final_ap_corr_map

            for field_name in reference_ap_corr_map:
                with self.subTest(scale=scale, field_name=field_name):
                    self.assertFloatsAlmostEqual(
                        final_ap_corr_map[field_name],
                        reference_ap_corr_map[field_name],
                        atol=1e-9,
                    )

    @lsst.utils.tests.methodParameters(
        do_coadd_inverse_ap_corr=[True, False],
    )
    def test_constant_weight(self, do_coadd_inverse_ap_corr):
        """Test that the coadd aperture corrections stay the same when the
        weights are all equal.
        """
        stacker = CoaddApCorrMapStacker(
            evaluation_point=self.evaluation_point,
            do_coadd_inverse_ap_corr=do_coadd_inverse_ap_corr,
        )

        for idx in range(self.visit_count):
            stacker.add(
                self.ap_corr_map_list[idx],
                weight=1.0,
            )

        reference_ap_corr_map = stacker.final_ap_corr_map

        stacker = CoaddApCorrMapStacker(
            evaluation_point=self.evaluation_point,
            do_coadd_inverse_ap_corr=do_coadd_inverse_ap_corr,
        )
        for idx in range(self.visit_count):
            stacker.add(
                self.ap_corr_map_list[idx],
                weight=2.0,
            )
        final_ap_corr_map = stacker.final_ap_corr_map

        for algorithm_name in stacker.ap_corr_names:
            field_name = f"{algorithm_name}_instFlux"
            with self.subTest(field_name=field_name):
                self.assertFloatsAlmostEqual(
                    final_ap_corr_map[field_name],
                    reference_ap_corr_map[field_name],
                    atol=1e-8,
                )

            field_name = f"{algorithm_name}_instFluxErr"
            # Errors are down by sqrt(visit_count).
            conversion_factor = 1.0  # np.sqrt(self.visit_count)
            with self.subTest(field_name=field_name):
                self.assertFloatsAlmostEqual(
                    final_ap_corr_map[field_name] / conversion_factor,
                    reference_ap_corr_map[field_name],
                    atol=1e-8,
                )

    @lsst.utils.tests.methodParameters(
        do_coadd_inverse_ap_corr=[True, False],
    )
    def test_constant_apcorr(self, do_coadd_inverse_ap_corr):
        """Test that the coadd aperture corrections are constant when the
        input aperture corrections are constant.
        """
        stacker = CoaddApCorrMapStacker(
            evaluation_point=self.evaluation_point,
            do_coadd_inverse_ap_corr=do_coadd_inverse_ap_corr,
        )

        weights = np.random.rand(self.visit_count)
        initial_ap_corr_map = self.ap_corr_map_list[0]
        for idx in range(self.visit_count):
            stacker.add(
                initial_ap_corr_map,  # Add the same over and over again.
                weight=weights[idx],
            )

        final_ap_corr_map = stacker.final_ap_corr_map

        for algorithm_name in stacker.ap_corr_names:
            field_name = f"{algorithm_name}_instFlux"
            with self.subTest(field_name=field_name):
                self.assertFloatsAlmostEqual(
                    final_ap_corr_map[field_name],
                    initial_ap_corr_map[field_name].evaluate(geom.Point2D(self.evaluation_point)),
                    atol=1e-8,
                )

            field_name = f"{algorithm_name}_instFluxErr"
            # Errors are down by conversion_factor.
            conversion_factor = np.sqrt(np.sum(weights**2)) / np.sum(weights)
            with self.subTest(field_name=field_name):
                self.assertFloatsAlmostEqual(
                    final_ap_corr_map[field_name] / conversion_factor,
                    initial_ap_corr_map[field_name].evaluate(geom.Point2D(self.evaluation_point)),
                    atol=1e-8,
                )

        # Evaluate aperture corrections at several points.
        evaluate_at_x = np.linspace(self.bbox.getMinX(), self.bbox.getMaxX(), 5)
        evaluate_at_y = np.linspace(self.bbox.getMinY(), self.bbox.getMaxY(), 5)
        for algorithm_name in stacker.ap_corr_names:
            field_name = f"{algorithm_name}_instFlux"
            reference_value = initial_ap_corr_map[field_name].evaluate(geom.Point2D(self.evaluation_point))
            aperture_correction_values = initial_ap_corr_map[field_name].evaluate(
                evaluate_at_x,
                evaluate_at_y,
            )
            with self.subTest(field_name=field_name):
                np.testing.assert_array_equal(aperture_correction_values, reference_value)


class TestCoaddApCorrMapMemoryTestCase(lsst.utils.tests.MemoryTestCase):
    """Test the CoaddApCorrMapStacker for memory leaks."""


def setup_module(module):  # noqa: D103
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
