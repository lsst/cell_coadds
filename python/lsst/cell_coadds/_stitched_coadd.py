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

from collections.abc import Iterator, Set
from functools import partial
from typing import TYPE_CHECKING

from lsst.afw.image import ExposureF, FilterLabel, PhotoCalib
from lsst.geom import Box2I, Point2I

from ._common_components import CoaddUnits, CommonComponents, CommonComponentsProperties
from ._grid_container import GridContainer
from ._image_planes import ImagePlanes, ViewImagePlanes
from ._stitched_aperture_correction import StitchedApertureCorrection
from ._stitched_image_planes import StitchedImagePlanes
from ._stitched_psf import StitchedPsf
from ._uniform_grid import UniformGrid
from .typing_helpers import StitchedCoaddApCorrMap

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
        ``cell_coadd.outer_bbox``.

    Notes
    -----
    This class simply inserts subimages from each cell into the full image,
    doing so when an attribute is first accessed to avoid stitching together
    planes that may never be accessed.

    A `StitchedCoadd` cannot be serialized in FITS format directly.  Instead,
    the recommended way is to serialize the `MultipleCellCoadd` instance that
    was used to construct the object and reconstruct the `StitchedCoadd` by
    calling the `stitch` method on it. A less recommended way is to call the
    `asExposure` method to get an `lsst.afw.image.Exposure` object and persist
    that to the disk.
    """

    def __init__(self, cell_coadd: MultipleCellCoadd, *, bbox: Box2I | None = None):
        super().__init__()
        if bbox is None:
            bbox = cell_coadd.outer_bbox
        elif not cell_coadd.outer_bbox.contains(bbox):
            raise ValueError(
                f"Cell coadd inner bounding box {cell_coadd.outer_bbox} does not "
                f"contain stitch target area {bbox}."
            )
        self._bbox = bbox
        self._cell_coadd = cell_coadd
        self._psf: StitchedPsf | None = None
        # self._ap_corr_map: StitchedCoaddApCorrMap | None = None
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
        x_max, y_max = self._cell_coadd.cells.last.identifiers.cell

        for cell in self._cell_coadd.cells.values():
            bbox = cell.inner.bbox
            if cell.identifiers.cell.x == 0:
                # This is a special case for the first column of cells.
                bbox.include(Point2I(cell.outer.bbox.beginX, cell.outer.bbox.centerY))
            elif cell.identifiers.cell.x == x_max:
                # This is a special case for the last column of cells.
                bbox.include(Point2I(cell.outer.bbox.endX, cell.outer.bbox.centerY))

            if cell.identifiers.cell.y == 0:
                # This is a special case for the last row of cells.
                bbox.include(Point2I(cell.outer.bbox.centerX, cell.outer.bbox.beginY))
            elif cell.identifiers.cell.y == y_max:
                # This is a special case for the first row of cells.
                bbox.include(Point2I(cell.outer.bbox.centerX, cell.outer.bbox.endY))

            bbox.clip(cell.outer.bbox)
            make_view = partial(cell.make_view, bbox=bbox)
            yield ViewImagePlanes(cell.outer, bbox=bbox, make_view=make_view)

    @property
    def n_noise_realizations(self) -> int:
        # Docstring inherited.
        return self._cell_coadd.n_noise_realizations

    @property
    def mask_fraction_names(self) -> Set[str]:
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
            result.metadata["BUNIT"] = "nJy"

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

    # @property
    # def ap_corr_map(self) -> StitchedCoaddApCorrMap:
        """Stitch the aperture correction maps from the cell coadd.

        This converts the aperture correction maps from each cell into a single
        `ApCorrMap` that quacks like `lsst.afw.image.ApCorrMap`. The resulting
        object has the fields to correct as keys and a
        `StitchedApertureCorrection` for each field.

        Notes
        -----
        These cannot be attached to an `~lsst.afw.image.Exposure` object.
        """
        # if self._ap_corr_map is None:
            # ap_corr_map: dict[str, StitchedApertureCorrection] = {}
            # field_names = self._cell_coadd.cells.first.aperture_corrected_algorithms
            # for field_name in field_names:
                # gc = GridContainer[float](shape=self.grid.shape)
                # for scc in self._cell_coadd.cells.values():
                    # gc[scc.identifiers.cell] = scc.aperture_correction_map[field_name]
                # ap_corr_map[field_name] = StitchedApertureCorrection(self.grid, gc)

            # self._ap_corr_map = ap_corr_map

        # return self._ap_corr_map
