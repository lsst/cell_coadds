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
from typing import Iterable, Mapping

import lsst.cell_coadds.test_utils as test_utils
import lsst.geom as geom
import lsst.meas.base.tests
import lsst.utils.tests
import numpy as np
from lsst.afw.geom import Quadrupole
from lsst.afw.image import ImageF
from lsst.cell_coadds import (
    CellCoaddFitsReader,
    CellIdentifiers,
    CoaddUnits,
    CommonComponents,
    ExplodedCoadd,
    MultipleCellCoadd,
    OwnedImagePlanes,
    PatchIdentifiers,
    SingleCellCoadd,
    StitchedCoadd,
    UniformGrid,
)
from lsst.meas.algorithms import SingleGaussianPsf
from lsst.skymap import Index2D


class BaseMultipleCellCoaddTestCase(lsst.utils.tests.TestCase):
    """A base class that provides a common set of methods."""

    psf_size: int
    psf_sigmas: Mapping[Index2D, float]
    border_size: int
    inner_size: int
    outer_size: int
    test_positions: Iterable[tuple[geom.Point2D, Index2D]]
    exposures: Mapping[Index2D, lsst.afw.image.ExposureF]
    multiple_cell_coadd: MultipleCellCoadd

    @classmethod
    def setUpClass(cls) -> None:
        """Setup a multiple cell coadd with 2x2 cells."""
        np.random.seed(42)
        data_id = test_utils.generate_data_id()
        common = CommonComponents(
            units=CoaddUnits.nJy,
            wcs=test_utils.generate_wcs(),
            band=data_id["band"].name,
            identifiers=PatchIdentifiers.from_data_id(data_id),
        )

        nx, ny = 3, 2
        cls.psf_sigmas = {
            Index2D(x=0, y=0): 1.2,
            Index2D(x=0, y=1): 0.7,
            Index2D(x=1, y=0): 0.9,
            Index2D(x=1, y=1): 1.1,
            Index2D(x=2, y=0): 1.3,
            Index2D(x=2, y=1): 0.8,
        }

        cls.border_size = 5
        # In practice, we expect this to be squares.
        # To check for any possible x, y swap, we use different values.
        cls.psf_size_x, cls.psf_size_y = 21, 23
        cls.inner_size_x, cls.inner_size_y = 17, 15
        cls.outer_size_x = cls.inner_size_x + 2 * cls.border_size
        cls.outer_size_y = cls.inner_size_y + 2 * cls.border_size

        patch_outer_bbox = geom.Box2I(
            geom.Point2I(0, 0), geom.Extent2I(nx * cls.inner_size_x, ny * cls.inner_size_y)
        )
        patch_outer_bbox.grow(cls.border_size)

        # Add one star and one galaxy per quadrant.
        # The mapping of positions to cell indices assume inner_size = (17, 15)
        # and border_size = 5. If that is changed, these values need an update.
        sources = (
            # flux, centroid, shape
            (1000.0, geom.Point2D(6.3, 7.2), None),
            (2500.0, geom.Point2D(16.8, 18.3), None),
            (1500.0, geom.Point2D(21.2, 5.1), None),
            (3200.0, geom.Point2D(16.1, 23.9), None),
            (1800.0, geom.Point2D(44.7, 8.9), None),
            (2100.0, geom.Point2D(34.1, 19.2), None),
            (900.0, geom.Point2D(9.1, 13.9), Quadrupole(2.5, 1.5, 0.8)),
            (1250.0, geom.Point2D(19.3, 11.2), Quadrupole(1.5, 2.5, 0.75)),
            (2100.0, geom.Point2D(5.1, 21.2), Quadrupole(1.7, 1.9, 0.05)),
            (2800.0, geom.Point2D(24.1, 19.2), Quadrupole(1.9, 1.7, 0.1)),
            (2350.0, geom.Point2D(40.3, 13.9), Quadrupole(1.8, 1.8, -0.4)),
            (4999.0, geom.Point2D(45.8, 22.0), Quadrupole(1.6, 1.2, 0.2)),
        )

        # The mapping of positions to cell indices assume inner_size = (17, 15)
        # and border_size = 5. If that is changed, these values need an update.
        cls.test_positions = (
            (geom.Point2D(5, 4), Index2D(x=0, y=0)),  # inner point in lower left
            (geom.Point2D(6, 24), Index2D(x=0, y=1)),  # inner point in upper left
            (geom.Point2D(25.2, 7.8), Index2D(x=1, y=0)),  # inner point in lower middle
            (geom.Point2D(23, 22), Index2D(x=1, y=1)),  # inner point in upper middle
            (geom.Point2D(39, 9.4), Index2D(x=2, y=0)),  # inner point in lower right
            (geom.Point2D(44, 24), Index2D(x=2, y=1)),  # inner point in upper right
            # Some points that lie on the border
            (geom.Point2D(33, 24), Index2D(x=1, y=1)),  # inner point in upper right
            (geom.Point2D(46, 0), Index2D(x=2, y=0)),  # inner point in lower right
            (geom.Point2D(19, 16), Index2D(x=1, y=1)),  # inner point in upper middle
            (geom.Point2D(17, 8), Index2D(x=1, y=0)),  # inner point in lower middle
            (geom.Point2D(0, 29), Index2D(x=0, y=1)),  # inner point in upper left
            (geom.Point2D(0, 0), Index2D(x=0, y=0)),  # inner point in lower left
        )

        schema = lsst.meas.base.tests.TestDataset.makeMinimalSchema()

        single_cell_coadds = []
        cls.exposures = dict.fromkeys(cls.psf_sigmas.keys())

        for x in range(nx):
            for y in range(ny):
                identifiers = CellIdentifiers(
                    cell=Index2D(x=x, y=y),
                    skymap=common.identifiers.skymap,
                    tract=common.identifiers.tract,
                    patch=common.identifiers.patch,
                    band=common.identifiers.band,
                )

                outer_bbox = geom.Box2I(
                    geom.Point2I(x * cls.inner_size_x, y * cls.inner_size_y),
                    geom.Extent2I(cls.inner_size_x, cls.inner_size_y),
                )
                outer_bbox.grow(cls.border_size)

                dataset = lsst.meas.base.tests.TestDataset(
                    patch_outer_bbox, psfSigma=cls.psf_sigmas[identifiers.cell]
                )

                for inst_flux, position, shape in sources:
                    dataset.addSource(inst_flux, position, shape)

                # Create a spatially varying variance plane.
                variance = ImageF(
                    # np.random.uniform returns an array with x-y flipped.
                    np.random.uniform(
                        0.8,
                        1.2,
                        (
                            ny * cls.inner_size_y + 2 * cls.border_size,
                            nx * cls.inner_size_x + 2 * cls.border_size,
                        ),
                    ).astype(np.float32),
                    xy0=outer_bbox.getMin(),
                )
                exposure, _ = dataset.realize(variance.getArray() ** 0.5, schema, randomSeed=123456789)
                cls.exposures[identifiers.cell] = exposure
                exposure = exposure[outer_bbox]
                image_plane = OwnedImagePlanes(
                    image=exposure.image, variance=exposure.variance, mask=exposure.mask
                )

                single_cell_coadds.append(
                    SingleCellCoadd(
                        outer=image_plane,
                        psf=SingleGaussianPsf(
                            cls.psf_size_x, cls.psf_size_y, cls.psf_sigmas[Index2D(x=x, y=y)]
                        ).computeKernelImage(outer_bbox.getCenter()),
                        inner_bbox=geom.Box2I(
                            geom.Point2I(x * cls.inner_size_x, y * cls.inner_size_y),
                            geom.Extent2I(cls.inner_size_x, cls.inner_size_y),
                        ),
                        inputs={
                            None,  # type: ignore [arg-type]
                        },
                        common=common,
                        identifiers=identifiers,
                    )
                )

        grid_bbox = geom.Box2I(
            geom.Point2I(0, 0), geom.Extent2I(nx * cls.inner_size_x, ny * cls.inner_size_y)
        )
        grid = UniformGrid.from_bbox_shape(grid_bbox, Index2D(x=nx, y=ny))

        cls.multiple_cell_coadd = MultipleCellCoadd(
            single_cell_coadds,
            grid=grid,
            outer_cell_size=geom.Extent2I(cls.outer_size_x, cls.outer_size_y),
            inner_bbox=None,
            common=common,
            psf_image_size=geom.Extent2I(cls.psf_size_x, cls.psf_size_y),
        )

    @classmethod
    def tearDownClass(cls) -> None:
        del cls.multiple_cell_coadd
        del cls.exposures
        super().tearDownClass()

    def assertMultipleCellCoaddsEqual(self, mcc1: MultipleCellCoadd, mcc2: MultipleCellCoadd) -> None:
        """Check the equality of two instances of `MultipleCellCoadd`.

        Parameters
        ----------
        mcc1 : `MultipleCellCoadd`
            The MultipleCellCoadd created by reading a FITS file.
        mcc2 : `MultipleCellCoadd`
            The reference MultipleCellCoadd for comparison.
        """
        self.assertEqual(mcc1.band, mcc2.band)
        self.assertEqual(mcc1.identifiers, mcc2.identifiers)
        self.assertEqual(mcc1.inner_bbox, mcc2.inner_bbox)
        self.assertEqual(mcc1.outer_bbox, mcc2.outer_bbox)
        self.assertEqual(mcc1.outer_cell_size, mcc2.outer_cell_size)
        self.assertEqual(mcc1.mask_fraction_names, mcc2.mask_fraction_names)
        self.assertEqual(mcc1.n_noise_realizations, mcc2.n_noise_realizations)
        self.assertEqual(mcc1.psf_image_size, mcc2.psf_image_size)
        self.assertEqual(mcc1.units, mcc2.units)
        self.assertEqual(mcc1.wcs.getFitsMetadata().toString(), mcc2.wcs.getFitsMetadata().toString())

        # Check that the individual cells are identical.
        self.assertEqual(mcc1.cells.keys(), mcc2.cells.keys())
        for idx in mcc1.cells.keys():
            self.assertImagesEqual(mcc1.cells[idx].outer.image, mcc2.cells[idx].outer.image)
            self.assertMasksEqual(mcc1.cells[idx].outer.mask, mcc2.cells[idx].outer.mask)
            self.assertImagesEqual(mcc1.cells[idx].outer.variance, mcc2.cells[idx].outer.variance)
            self.assertImagesEqual(mcc1.cells[idx].psf_image, mcc2.cells[idx].psf_image)

            self.assertEqual(mcc1.cells[idx].band, mcc1.band)
            self.assertEqual(mcc1.cells[idx].common, mcc1.common)
            self.assertEqual(mcc1.cells[idx].units, mcc2.units)
            self.assertEqual(mcc1.cells[idx].wcs, mcc1.wcs)
            # Identifiers differ because of the ``cell`` component.
            # Check the other attributes within the identifiers.
            for attr in ("skymap", "tract", "patch", "band"):
                self.assertEqual(getattr(mcc1.cells[idx].identifiers, attr), getattr(mcc1.identifiers, attr))


class MultipleCellCoaddTestCase(BaseMultipleCellCoaddTestCase):
    """Test the construction and interfaces of MultipleCellCoadd."""

    def test_fits(self):
        """Test that we can write a coadd to a FITS file and read it."""
        with lsst.utils.tests.getTempFilePath(".fits") as filename:
            self.multiple_cell_coadd.write_fits(filename)
            mcc1 = MultipleCellCoadd.read_fits(filename)  # Test the readFits method.

            # Test the reader class.
            reader = CellCoaddFitsReader(filename)
            mcc2 = reader.readAsMultipleCellCoadd()

            wcs = reader.readWcs()

        self.assertMultipleCellCoaddsEqual(mcc1, self.multiple_cell_coadd)
        self.assertMultipleCellCoaddsEqual(mcc2, self.multiple_cell_coadd)
        # By transititve property of equality, mcc1 == mcc2.

        self.assertEqual(self.multiple_cell_coadd.band, self.multiple_cell_coadd.common.band)
        self.assertEqual(
            wcs.getFitsMetadata().toString(), self.multiple_cell_coadd.wcs.getFitsMetadata().toString()
        )


class ExplodedCoaddTestCase(BaseMultipleCellCoaddTestCase):
    """Test the construction and methods of an ExplodedCoadd instance."""

    exploded_coadd: ExplodedCoadd

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.exploded_coadd = cls.multiple_cell_coadd.explode()

    @classmethod
    def tearDownClass(cls) -> None:
        del cls.exploded_coadd
        super().tearDownClass()

    def test_exploded_psf_image(self):
        """Show that psf_image sizes are absurd."""
        self.assertEqual(
            self.exploded_coadd.psf_image.getBBox().getDimensions(),
            geom.Extent2I(3 * self.psf_size_x, 2 * self.psf_size_y),
        )
        for pad_psfs_with in (-999, -4, 0, 4, 8, 21, 40, 100):
            exploded_coadd = self.multiple_cell_coadd.explode(pad_psfs_with=pad_psfs_with)
            self.assertEqual(
                exploded_coadd.psf_image.getBBox().getDimensions(),
                geom.Extent2I(3 * self.outer_size_x, 2 * self.outer_size_y),
            )


class StitchedCoaddTestCase(BaseMultipleCellCoaddTestCase):
    """Test the construction and methods of a StitchedCoadd instance."""

    stitched_coadd: StitchedCoadd

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.stitched_coadd = cls.multiple_cell_coadd.stitch()

    @classmethod
    def tearDownClass(cls) -> None:
        del cls.stitched_coadd
        super().tearDownClass()

    def test_computeBBox(self):
        """Test the computeBBox method for a StitchedPsf object."""
        stitched_psf = self.stitched_coadd.psf

        psf_bbox = geom.Box2I(
            geom.Point2I(-(self.psf_size_x // 2), -(self.psf_size_y // 2)),
            geom.Extent2I(self.psf_size_x, self.psf_size_y),
        )

        for position, _ in self.test_positions:
            bbox = stitched_psf.computeBBox(position)
            self.assertEqual(bbox, psf_bbox)

    def test_computeShape(self):
        """Test the computeShape method for a StitchedPsf object."""
        stitched_psf = self.stitched_coadd.psf
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

    def test_computeKernelImage(self):
        """Test the computeKernelImage method for a StitchedPsf object."""
        stitched_psf = self.stitched_coadd.psf
        psf_bbox = geom.Box2I(
            geom.Point2I(-(self.psf_size_x // 2), -(self.psf_size_y // 2)),
            geom.Extent2I(self.psf_size_x, self.psf_size_y),
        )

        for position, cell_index in self.test_positions:
            image1 = stitched_psf.computeKernelImage(position)
            image2 = SingleGaussianPsf(
                self.psf_size_x, self.psf_size_y, self.psf_sigmas[cell_index]
            ).computeKernelImage(position)
            self.assertImagesEqual(image1, image2)
            self.assertEqual(image1.getBBox(), psf_bbox)

    def test_computeImage(self):
        """Test the computeImage method for a StitchedPsf object."""
        stitched_psf = self.stitched_coadd.psf
        psf_extent = geom.Extent2I(self.psf_size_x, self.psf_size_y)

        for position, cell_index in self.test_positions:
            image1 = stitched_psf.computeImage(position)
            image2 = SingleGaussianPsf(
                self.psf_size_x, self.psf_size_y, self.psf_sigmas[cell_index]
            ).computeImage(position)
            self.assertImagesEqual(image1, image2)
            self.assertEqual(image1.getBBox().getDimensions(), psf_extent)

    def test_computeImage_computeKernelImage(self):
        """Test that computeImage called at integer points gives the same
        result as calling computeKernelImage.
        """
        stitched_psf = self.stitched_coadd.psf
        for position, cell_index in self.test_positions:
            pos = geom.Point2D(geom.Point2I(position))  # round to integer
            image1 = stitched_psf.computeKernelImage(pos)
            image2 = stitched_psf.computeImage(pos)
            self.assertImagesEqual(image1, image2)

    def test_computeApetureFlux(self):
        """Test the computeApertureFlux method for a StitchedPsf object."""
        stitched_psf = self.stitched_coadd.psf
        for position, cell_index in self.test_positions:
            flux1sigma = stitched_psf.computeApertureFlux(self.psf_sigmas[cell_index], position=position)
            self.assertAlmostEqual(flux1sigma, 0.39, delta=5e-2)

            flux3sigma = stitched_psf.computeApertureFlux(
                3.0 * self.psf_sigmas[cell_index], position=position
            )
            self.assertAlmostEqual(flux3sigma, 0.97, delta=2e-2)

    def test_asExposure(self):
        """Test the asExposure method for a StitchedCoadd object."""
        exposure = self.stitched_coadd.asExposure()

        for y in range(2):
            for x in range(3):
                bbox = geom.Box2I(
                    geom.Point2I(x * self.inner_size_x, y * self.inner_size_y),
                    geom.Extent2I(self.inner_size_x, self.inner_size_y),
                )
                index = Index2D(x=x, y=y)
                self.assertImagesEqual(exposure.image[bbox], self.exposures[index].image[bbox])
                self.assertImagesEqual(exposure.variance[bbox], self.exposures[index].variance[bbox])
                self.assertImagesEqual(exposure.mask[bbox], self.exposures[index].mask[bbox])


class TestMemory(lsst.utils.tests.MemoryTestCase):
    pass


def setup_module(module):
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
