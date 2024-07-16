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

from typing import TYPE_CHECKING, Self

import lsst.shoefits as shf
import numpy as np

from ._image_planes import ImagePlanes
from ._to_upstream import PixelShape
from ._uniform_grid import UniformGrid

if TYPE_CHECKING:
    from ._multiple_cell_coadd import MultipleCellCoadd


class ExplodedCoadd(ImagePlanes):
    """A coadd that stitches together the outer regions of the cells in a
    `MultipleCellCoadd` (including multiple values for most pixels).

    Notes
    -----
    `ExplodedCoadd` is intended primarily for serialization and visualization,
    and with this in mind it stores its PSF as a single stitched image as well.

    An `ExplodedCoadd` and a `MultipleCellCoadd` may share data, but generally
    only if the `MultipleCellCoadd` is constructed from the `ExplodedCoadd`,
    not the other way around.
    """

    grid: UniformGrid
    """Object that defines the piecewise grid (of outer cell regions) that
    this object stitches together.

    This grid always starts at ``(0, 0)``; because there is no way to align
    this "exploded" grid with the inner cell grid over more than one cell, no
    attempt is made to align it with the inner cell grid's overall offset.
    """

    psf_cell_size: PixelShape
    """Shape of one cell's PSF model image in `psf_image`."""

    psf_image: shf.Image
    """A stitched-together image of the PSF models for each cell."""

    @classmethod
    def from_cell_coadd(cls, cell_coadd: MultipleCellCoadd) -> Self:
        """Build an `ExplodedCoadd` from a `MultipleCellCoadd`.

        Parameters
        ----------
        cell_coadd : `MultipleCellCoadd`
            Cell-based coadd to stitch together.

        Returns
        -------
        exploded : `ExplodedCoadd`
            Stitched outer-cell coadd.
        """
        exploded_grid = UniformGrid.from_cell_size_shape(
            cell_size=cell_coadd.outer_cell_size, shape=cell_coadd.grid.shape, padding=0
        )
        psf_image_bbox = shf.Box.factory[
            : cell_coadd.psf_image_size.y * exploded_grid.shape.y,
            : cell_coadd.psf_image_size.x * exploded_grid.shape.x,
        ]
        psf_image = shf.Image(0.0, bbox=psf_image_bbox, dtype=np.float32)
        result = cls.from_bbox(
            bbox=exploded_grid.bbox,
            mask_schema=cell_coadd.mask_schema,
            include_mask_fractions=cell_coadd.has_mask_fractions,
            n_noise_realizations=cell_coadd.n_noise_realizations,
            psf_native_size=cell_coadd.psf_image_size,
            psf_image=psf_image,
        )
        for index, cell in cell_coadd.cells.items():
            exploded_cell_bbox = exploded_grid.bbox_of(index)
            exploded_cell_view = ImagePlanes.view_of(result, exploded_cell_bbox)
            exploded_cell_view.copy_from(cell.outer)
            iy = index.y * cell_coadd.psf_image_size.y
            ix = index.x * cell_coadd.psf_image_size.x
            psf_image.array[iy : iy + cell_coadd.psf_image_size.y, ix : ix + cell_coadd.psf_image_size.x] = (
                cell.psf_image.array
            )
        return result
