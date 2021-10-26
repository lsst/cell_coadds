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

__all__ = ("SingleCellCoadd",)


from lsst.afw.image import ImageF
from lsst.geom import Box2I

from ._image_planes import ImagePlanes, OwnedImagePlanes, ViewImagePlanes


class SingleCellCoadd:
    """A single coadd cell, built only from input images that completely
    overlap that cell.

    Parameters
    ----------
    outer: `OwnedImagePlanes`
        The actual coadded images.
    psf : `ImageF`
        The coadded PSF image.
    inner_bbox: `Box2I`
        The bounding box of the inner region of this cell; must be disjoint
        with but adjacent to all other cell inner regions.

    Notes
    -----
    At present we assume a single PSF image per cell is sufficient to capture
    spatial variability, which seems adequate given the results we have so far
    and the cell sizes we intend to use.
    """

    def __init__(self, outer: OwnedImagePlanes, psf: ImageF, inner_bbox: Box2I):
        assert outer.bbox.contains(
            inner_bbox
        ), f"Cell inner bbox {inner_bbox} is not contained by outer bbox {outer.bbox}."
        self._outer = outer
        self._psf = psf
        self._inner = ViewImagePlanes(outer, inner_bbox)

    @property
    def inner(self) -> ImagePlanes:
        """Image planes within the inner region of this cell that is disjoint
        with all other cell inner regions.
        """
        return self._inner

    @property
    def outer(self) -> ImagePlanes:
        """Image planes within the full outer region of this cell."""
        return self._outer

    @property
    def psf_image(self) -> ImageF:
        """The coadded PSF image."""
        return self._psf
