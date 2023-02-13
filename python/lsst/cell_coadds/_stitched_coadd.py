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

from typing import TYPE_CHECKING, AbstractSet, Iterator

from lsst.afw.image import ExposureF, FilterLabel, PhotoCalib
from lsst.geom import Box2I

from ._cell_coadds import StitchedPsf
from ._common_components import CoaddUnits, CommonComponents, CommonComponentsProperties
from ._image_planes import ImagePlanes
from ._stitched_image_planes import StitchedImagePlanes
from ._uniform_grid import UniformGrid

if TYPE_CHECKING:
    from ._multiple_cell_coadd import MultipleCellCoadd


class StitchedCoadd(StitchedImagePlanes, CommonComponentsProperties):
    """A lazy-evaluation coadd that stitches together images from adjacent
    cells.

    Parameters
    ----------
    cell_coadd : `MultipleCellCoadd`
        Cell-based coadd to stitch together.
    bbox : `Box2I`, optional
        The region over which a contiguous coadd is desired.  Defaults to
        ``cell_coadd.inner_bbox``.

    Notes
    -----
    This class simply inserts subimages from each cell into the full image,
    doing so when an attribute is first accessed to avoid stitching together
    planes that may never be accessed.
    """

    def __init__(self, cell_coadd: MultipleCellCoadd, *, bbox: Box2I | None = None):
        super().__init__()
        if bbox is None:
            bbox = cell_coadd.inner_bbox
        elif not cell_coadd.inner_bbox.contains(bbox):
            raise ValueError(
                f"Cell coadd inner bounding box {cell_coadd.inner_bbox} does not "
                f"contain stitch target area {bbox}."
            )
        self._bbox = bbox
        self._cell_coadd = cell_coadd
        self._psf: StitchedPsf | None = None
        self._common = cell_coadd.common

    @property
    def bbox(self) -> Box2I:
        # Docstring inherited.
        return self._bbox

    @property
    def grid(self) -> UniformGrid:
        """Object that defines the piecewise grid (of inner cell regions) that
        this object stitches together.

        This may include cells outside the region covered by these image
        planes.
        """
        return self._cell_coadd.grid

    def _iter_cell_planes(self) -> Iterator[ImagePlanes]:
        # Docstring inherited.
        for cell in self._cell_coadd.cells:
            yield cell.inner

    @property
    def n_noise_realizations(self) -> int:
        # Docstring inherited.
        return self._cell_coadd.n_noise_realizations

    @property
    def mask_fraction_names(self) -> AbstractSet[str]:
        # Docstring inherited.
        return self._cell_coadd.mask_fraction_names

    @property
    def psf(self) -> StitchedPsf:
        """The piecewise PSF of this image."""
        if self._psf is None:
            self._psf = StitchedPsf(
                self._cell_coadd.cells.rebuild_transformed(lambda cell: cell.psf_image.convertD()),
                self._cell_coadd.grid,
            )
        return self._psf

    @property
    def common(self) -> CommonComponents:
        # Docstring inherited.
        return self._cell_coadd.common

    def asExposure(self) -> ExposureF:
        """Return an `lsst.afw.image.Exposure` view of this piecewise image."""
        result = ExposureF(self.asMaskedImage())
        # Exposure components derived from "common" components are all simple.
        result.setWcs(self._cell_coadd.wcs)
        result.setFilter(FilterLabel(band=self.band))
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
        # result.setId(self._cell_coadd.identifiers.patch)

        # We could add CoaddInputs here, but without WCS, PSF, etc in them;
        # it's not clear that's good enough or even useful, given that the cell
        # provide a much more clean view of what the inputs are at any given
        # point.

        # PSF is the first of many components that need piecewise
        # implementations.  More to do here for at least aperture corrections
        # and transmission curves.
        result.setPsf(self.psf)

        return result
