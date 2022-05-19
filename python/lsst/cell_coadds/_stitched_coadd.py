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

from __future__ import annotations

__all__ = ("StitchedCoadd",)

from typing import TYPE_CHECKING, Optional

from lsst.afw.image import ExposureF, FilterLabel, ImageF, Mask, PhotoCalib
from lsst.geom import Box2I

from ._cell_coadds import StitchedPsf, UniformGrid
from ._common_components import CoaddUnits, CommonComponents, CommonComponentsProperties
from ._image_planes import ImagePlanes, ImagePlaneTag

from . import typing_helpers

if TYPE_CHECKING:
    from ._multiple_cell_coadd import MultipleCellCoadd


class StitchInnerCoaddCells:
    def __init__(self, coadd: MultipleCellCoadd, bbox: Box2I):
        self._coadd = coadd
        self._bbox = bbox

    def __call__(self, tag: ImagePlaneTag) -> typing_helpers.ImageLike:
        result = tag.image_type(self._bbox)
        for cell in self._coadd.cells:
            common_bbox = cell.inner.bbox.clippedTo(self._bbox)
            if not common_bbox.isEmpty():
                result[common_bbox] = tag.get(cell.inner)[common_bbox]
        return result


class StitchedCoadd(ImagePlanes, CommonComponentsProperties):
    """Coadd that stitches together images from adjacent cells, considering
    their inner regions only.
    """

    def __init__(
        self,
        image: ImageF,
        mask: Mask,
        variance: ImageF,
        grid: UniformGrid,
        psf: StitchedPsf,
        common: CommonComponents,
    ):
        super().__init__(image, mask, variance)
        self._grid = grid
        self._psf = psf
        self._common = common

    @classmethod
    def build(cls, coadd: MultipleCellCoadd, bbox: Optional[Box2I] = None) -> StitchedCoadd:
        if bbox is None:
            bbox = coadd.inner_bbox
        psf = StitchedPsf(
            coadd.cells.rebuild_transformed(lambda cell: cell.psf_image.convertD()).finish(),
            coadd.grid,
        )
        return cls.from_callback(
            StitchInnerCoaddCells(bbox, coadd),
            coadd.mask_fraction_names,
            coadd.n_noise_realizations,
            grid=coadd.grid,
            psf=psf,
            common=coadd.common,
        )

    @property
    def grid(self) -> UniformGrid:
        """Object that defines the piecewise grid (of inner cell regions) that
        this object stitches together.

        This may include cells outside the region covered by these image
        planes.
        """
        return self._grid

    @property
    def psf(self) -> StitchedPsf:
        """The piecewise PSF of this image."""
        return self._psf

    @property
    def common(self) -> CommonComponents:
        # Docstring inherited.
        return self._common

    def asExposure(self) -> ExposureF:
        """Return an `lsst.afw.image.Exposure` view of this piecewise image."""
        result = ExposureF(self.asMaskedImage())
        # Exposure components derived from "common" components are all simple.
        result.setWcs(self._common.wcs)
        result.setFilterLabel(FilterLabel(band=self.band))
        if self.units is CoaddUnits.nJy:
            result.setPhotoCalib(PhotoCalib(1.0))

        # We can't do result.setId here, because:
        #
        # - we don't know whether this should be a packed tract+patch+band ID
        #   or just a tract+patch ID;
        #
        # - we don't know how to pack the information we have anyway.
        #
        # Maybe DM-31924 will provide a solution to at least the latter.
        result.setId(self._common.identifiers.patch)

        # We could add CoaddInputs here, but without WCS, PSF, etc in them;
        # it's not clear that's good enough or even useful, given that the cell
        # provide a much more clean view of what the inputs are at any given
        # point.

        # PSF is the first of many components that need piecewise
        # implementations.  More to do here for at least aperture corrections
        # and transmission curves.
        result.setPsf(self.psf)

        return result
