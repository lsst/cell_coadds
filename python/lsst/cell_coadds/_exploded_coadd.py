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

__all__ = ("ExplodedCoadd",)

from typing import TYPE_CHECKING, Optional

from lsst.geom import Box2I
from lsst.afw.image import ImageF, Mask

from ._cell_coadds import UniformGrid
from ._common_components import CommonComponents, CommonComponentsProperties
from ._image_planes import ImagePlaneTag, ImagePlanes
from . import typing_helpers

if TYPE_CHECKING:
    from ._multiple_cell_coadd import MultipleCellCoadd


class ExplodeOuterCoaddCells:
    def __init__(self, coadd: MultipleCellCoadd, grid: UniformGrid, bbox: Box2I):
        self._coadd = coadd
        self._grid = grid
        self._bbox = bbox

    def __call__(self, tag: ImagePlaneTag) -> typing_helpers.ImageLike:
        result = tag.image_type(self._bbox)
        for cell in self._coadd.cells:
            cell_bbox = self._grid.bbox_of(cell.identifiers.cell.index)
            common_bbox = cell_bbox.clippedTo(self._bbox)
            if not common_bbox.isEmpty():
                cell_image: typing_helpers.ImageLike = tag.get(cell.outer).clone(False)  # shallow copy
                cell_image.setXY0(cell_bbox.getMin())
                result[common_bbox] = cell_image[common_bbox]
        return result


def stitch_cell_psf_images(
    coadd: MultipleCellCoadd, exploded_grid: UniformGrid, pad_psfs_with: Optional[float] = None
) -> ImageF:
    if pad_psfs_with is None:
        psf_grid = UniformGrid(coadd.psf_image_size, coadd.grid.shape)
    elif (coadd.psf_image_size > coadd.outer_cell_size).any():
        raise ValueError(
            f"PSF image dimensions {coadd.psf_image_size} are larger than "
            f"outer cell dimensions {coadd.outer_cell_size}; cannot pad."
        )
    else:
        psf_grid = exploded_grid
    stitched_psf_image = ImageF(psf_grid.bbox)
    if pad_psfs_with is not None:
        stitched_psf_image.set(pad_psfs_with)
    for cell in coadd.cells:
        target_subimage = stitched_psf_image[psf_grid.bbox_of(cell.identifiers.cell.index)]
        target_dimensions = target_subimage.getDimensions()
        source_dimensions = cell.psf_image.getDimensions()
        offset = (target_dimensions - source_dimensions) // 2
        # Use numpy views instead of Image methods because the images
        # are not in the same coordinate system to begin with, so xy0
        # doesn't help us.
        target_subimage.array[
            offset.y : offset.y + source_dimensions.y, offset.x : offset.x + source_dimensions.x
        ] = cell.psf_image.array
    return stitched_psf_image, psf_grid


class ExplodedCoadd(ImagePlanes, CommonComponentsProperties):
    """A coadd that stitches together the outer regions of the cells in a
    `MultipleCellCoadd` (including multiple values for most pixels).

    Parameters
    ----------
    cell_coadd : `MultipleCellCoadd`
        Cell-based coadd to stitch together.
    pad_psfs_with : `float` or `None`, optional
        A floating-point value to pad PSF images with so each PSF-image cell
        has the same dimensions as the image (outer) cell it corresponds to.
        If `None`, PSF images will not be padded and the full PSF image will
        generally be smaller than the exploded image it corresponds to.
    """

    def __init__(
        self,
        image: ImageF,
        mask: Mask,
        variance: ImageF,
        grid: UniformGrid,
        psf_grid: UniformGrid,
        psf_image: ImageF,
        common: CommonComponents,
    ):
        super().__init__(image, mask, variance)
        self._grid = grid
        self._psf_grid = psf_grid
        self._psf_image = psf_image
        self._common = common

    @classmethod
    def build(cls, coadd: MultipleCellCoadd, pad_psfs_with: Optional[float] = None) -> ExplodedCoadd:
        grid = UniformGrid(coadd.outer_cell_size, coadd.grid.shape)
        psf_image, psf_grid = stitch_cell_psf_images(coadd, grid, pad_psfs_with)
        return cls.from_callback(
            ExplodeOuterCoaddCells(coadd, grid, grid.bbox),
            coadd.mask_fraction_names,
            coadd.n_noise_realizations,
            grid=grid,
            psf_grid=psf_grid,
            psf_image=psf_image,
            common=coadd.common,
        )

    @property
    def grid(self) -> UniformGrid:
        """Object that defines the piecewise grid (of outer cell regions) that
        this object stitches together.

        This grid always starts at `(0, 0)`; because there is no way to align
        this "exploded" grid with the inner cell grid over more than one cell,
        no attempt is made to align it with the inner cell grid's overall
        offset.
        """
        return self._grid

    @property
    def psf_grid(self) -> UniformGrid:
        """Object that describes the grid on which PSF model images are
        stitched together.
        """
        return self._psf_grid

    @property
    def psf_image(self) -> ImageF:
        """A stitched-together image of the PSF models for each cell."""
        return self._psf_image
