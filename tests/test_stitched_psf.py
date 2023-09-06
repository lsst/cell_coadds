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

import pickle
import unittest

import lsst.geom as geom
import lsst.meas.base.tests
import lsst.utils.tests
import numpy as np
from lsst.afw.detection import GaussianPsf
from lsst.afw.image import ImageD
from lsst.cell_coadds import GridContainer, StitchedPsf, UniformGrid
from lsst.meas.algorithms import SingleGaussianPsf
from lsst.skymap import Index2D


class StitchedPsfTestCase(lsst.utils.tests.TestCase):
    """Test the methods of StitchedPsf class."""

    @classmethod
    def setUpClass(cls) -> None:  # noqa: D102
        super().setUpClass()

        shape = Index2D(x=3, y=2)

        cls.psf_size = 15
        cls.psf_sigmas = {
            Index2D(x=0, y=0): 1.2,
            Index2D(x=0, y=1): 0.7,
            Index2D(x=1, y=0): 0.9,
            Index2D(x=1, y=1): 1.1,
            Index2D(x=2, y=0): 1.3,
            Index2D(x=2, y=1): 0.8,
        }

        gc = GridContainer[ImageD](shape)
        for key, sigma in cls.psf_sigmas.items():
            # It does not matter where we compute the kernel image.
            # Pick any point.
            gc[key] = GaussianPsf(cls.psf_size, cls.psf_size, sigma).computeKernelImage(geom.Point2D(1, 1))

        grid = UniformGrid(cell_size=geom.Extent2I(cls.psf_size, cls.psf_size), shape=shape)

        cls.psf = StitchedPsf(gc, grid)

        cls.test_positions = (
            (geom.Point2D(5, 4), Index2D(x=0, y=0)),  # inner point in lower left
            (geom.Point2D(6, 24), Index2D(x=0, y=1)),  # inner point in upper left
            (geom.Point2D(23.2, 7.8), Index2D(x=1, y=0)),  # inner point in lower middle
            (geom.Point2D(21, 22), Index2D(x=1, y=1)),  # inner point in upper middle
            (geom.Point2D(37, 9.4), Index2D(x=2, y=0)),  # inner point in lower right
            (geom.Point2D(42, 24), Index2D(x=2, y=1)),  # inner point in upper right
            # Some points that lie on the border
            (geom.Point2D(31, 24), Index2D(x=2, y=1)),  # inner point in upper right
            (geom.Point2D(44, 0), Index2D(x=2, y=0)),  # inner point in lower right
            (geom.Point2D(17, 16), Index2D(x=1, y=1)),  # inner point in upper middle
            (geom.Point2D(15, 8), Index2D(x=1, y=0)),  # inner point in lower middle
            (geom.Point2D(0, 29), Index2D(x=0, y=1)),  # inner point in upper left
            (geom.Point2D(0, 0), Index2D(x=0, y=0)),  # inner point in lower left
        )

    def test_resized(self):
        """Test that the resized method works as it should."""
        original_bbox = self.psf.computeBBox(self.psf.getAveragePosition())
        # Resizing to original size should leave everything unchanged.
        psf = self.psf.resized(self.psf_size, self.psf_size)
        self.assertEqual(psf.getAveragePosition(), self.psf.getAveragePosition())
        self.assertEqual(
            psf.computeBBox(psf.getAveragePosition()), self.psf.computeBBox(self.psf.getAveragePosition())
        )
        self.assertImagesEqual(
            psf.computeKernelImage(psf.getAveragePosition()),
            self.psf.computeKernelImage(self.psf.getAveragePosition()),
        )
        self.assertEqual(psf, self.psf)

        # Resize to a rectangular postage stamp
        psf = self.psf.resized(25, 21)
        self.assertEqual(self.psf.computeBBox(self.psf.getAveragePosition()), original_bbox)
        self.assertEqual(psf.computeBBox(psf.getAveragePosition()).getDimensions(), geom.Extent2I(25, 21))
        self.assertEqual(
            psf.computeBBox(psf.getAveragePosition()).getCenter(),
            self.psf.computeBBox(self.psf.getAveragePosition()).getCenter(),
        )
        self.assertNotEqual(psf, self.psf)

        # Test that resizing to even dimensions throws an error
        with self.assertRaises(ValueError):
            psf = self.psf.resized(26, 22)  # even, even
        with self.assertRaises(ValueError):
            psf = self.psf.resized(25, 22)  # odd, even
        with self.assertRaises(ValueError):
            psf = self.psf.resized(26, 21)  # even, odd

    def test_clone(self):
        """Test that the clone method works."""
        psf = self.psf
        cloned = psf.clone()
        # Check that the cloned version is an exact replica of the original.
        # self.assertEqual(psf, cloned)  # Should there be an __eq__ method?
        self.assertEqual(psf.getAveragePosition(), cloned.getAveragePosition())
        self.assertEqual(
            psf.computeBBox(psf.getAveragePosition()), cloned.computeBBox(cloned.getAveragePosition())
        )
        self.assertImagesEqual(
            psf.computeKernelImage(psf.getAveragePosition()),
            cloned.computeKernelImage(cloned.getAveragePosition()),
        )
        self.assertEqual(psf, cloned)
        psf = psf.resized(41, 41)
        # Now that one of them has been resized, they should not be equal.
        self.assertNotEqual(
            psf.computeBBox(psf.getAveragePosition()), cloned.computeBBox(cloned.getAveragePosition())
        )
        with self.assertRaises(TypeError):
            self.assertImagesEqual(
                psf.computeKernelImage(psf.getAveragePosition()),
                cloned.computeKernelImage(cloned.getAveragePosition()),
            )
        self.assertNotEqual(psf, cloned)

    def test_computeKernelImage(self):
        """Test the computeKernelImage method for a StitchedPsf object."""
        stitched_psf = self.psf
        psf_bbox = geom.Box2I(
            geom.Point2I(-(self.psf_size // 2), -(self.psf_size // 2)),
            geom.Extent2I(self.psf_size, self.psf_size),
        )

        for position, cell_index in self.test_positions:
            image1 = stitched_psf.computeKernelImage(position)
            image2 = SingleGaussianPsf(
                self.psf_size, self.psf_size, self.psf_sigmas[cell_index]
            ).computeKernelImage(position)
            # Small differences may exist due to differences in evaluating
            # GaussianPsf vs. SingleGaussianPsf
            self.assertImagesAlmostEqual(image1, image2, atol=2e-16)
            self.assertEqual(image1.getBBox(), psf_bbox)

    def test_computeBBox(self):
        """Test the computeBBox method for a StitchedPsf object."""
        psf = self.psf
        psf_bbox = geom.Box2I(
            geom.Point2I(-(self.psf_size // 2), -(self.psf_size // 2)),
            geom.Extent2I(self.psf_size, self.psf_size),
        )

        for position, _ in self.test_positions:
            bbox = psf.computeBBox(position)
            self.assertEqual(bbox, psf_bbox)

    def test_computeShape(self):
        """Test the results from the computeShape method on a StitchedPsf
        object matches the true input.
        """
        stitched_psf = self.psf
        for position, cell_index in self.test_positions:
            psf_shape = stitched_psf.computeShape(position)  # check we can compute shape
            self.assertIsNot(psf_shape.getIxx(), np.nan)
            self.assertIsNot(psf_shape.getIyy(), np.nan)
            self.assertIsNot(psf_shape.getIxy(), np.nan)

            # Moments measured from pixellated images are significantly
            # underestimated for small PSFs.
            if self.psf_sigmas[cell_index] >= 1.0:
                self.assertAlmostEqual(psf_shape.getIxx(), self.psf_sigmas[cell_index] ** 2, delta=1e-3)
                self.assertAlmostEqual(psf_shape.getIyy(), self.psf_sigmas[cell_index] ** 2, delta=1e-3)
                self.assertAlmostEqual(psf_shape.getIxy(), 0.0)

    def test_computeApertureFlux(self):
        """Test that the results from the computeApertureFlux method on a
        StitchedPsf object returns the analytical results for a Gaussian PSF.
        """
        stitched_psf = self.psf
        for position, cell_index in self.test_positions:
            flux1sigma = stitched_psf.computeApertureFlux(self.psf_sigmas[cell_index], position=position)
            self.assertAlmostEqual(flux1sigma, 0.39, delta=5e-2)

            flux3sigma = stitched_psf.computeApertureFlux(
                3.0 * self.psf_sigmas[cell_index], position=position
            )
            self.assertAlmostEqual(flux3sigma, 0.97, delta=2e-2)

    def test_computeImage(self):
        """Test the computeImage method for a StitchedPsf object produces
        the same result as that on GaussianPsf for Gaussian PSFs.
        """
        stitched_psf = self.psf
        psf_extent = geom.Extent2I(self.psf_size, self.psf_size)

        for position, cell_index in self.test_positions:
            image1 = stitched_psf.computeImage(position)
            image2 = GaussianPsf(self.psf_size, self.psf_size, self.psf_sigmas[cell_index]).computeImage(
                position
            )
            self.assertImagesEqual(image1, image2)
            self.assertEqual(image1.getBBox().getDimensions(), psf_extent)

    def test_computeImage_computeKernelImage(self):
        """Test that computeImage called at integer points gives the same
        result as calling computeKernelImage.
        """
        stitched_psf = self.psf
        for position, cell_index in self.test_positions:
            pos = geom.Point2D(geom.Point2I(position))  # round to integer
            image1 = stitched_psf.computeKernelImage(pos)
            image2 = stitched_psf.computeImage(pos)
            self.assertImagesEqual(image1, image2)

    def test_pickle(self):
        """Test that StitchedPsf objects can be pickled and unpickled."""
        self.assertTrue(self.psf.isPersistable())
        stream = pickle.dumps(self.psf)
        psf = pickle.loads(stream)
        self.assertEqual(psf.getAveragePosition(), self.psf.getAveragePosition())
        self.assertEqual(
            psf.computeBBox(psf.getAveragePosition()), self.psf.computeBBox(self.psf.getAveragePosition())
        )
        self.assertImagesEqual(
            psf.computeKernelImage(psf.getAveragePosition()),
            self.psf.computeKernelImage(self.psf.getAveragePosition()),
        )
        self.assertEqual(psf, self.psf)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    """Test for memory/resource leaks."""


def setup_module(module):  # noqa: D103
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
