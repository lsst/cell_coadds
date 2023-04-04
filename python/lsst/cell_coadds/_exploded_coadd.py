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

from typing import TYPE_CHECKING, AbstractSet, Iterator

from lsst.afw.image import ImageF
from lsst.geom import Box2I

from ._image_planes import ImagePlanes, ViewImagePlanes
from ._stitched_image_planes import StitchedImagePlanes
from ._uniform_grid import UniformGrid
from .typing_helpers import ImageLike

if TYPE_CHECKING:
    from ._multiple_cell_coadd import MultipleCellCoadd


class ExplodedCoadd(StitchedImagePlanes):
    """A lazy-evaluation coadd that stitches together the outer regions of the
    cells in a `MultipleCellCoadd` (including multiple values for most pixels).

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

    def __init__(self, cell_coadd: MultipleCellCoadd, *, pad_psfs_with: float | None = None):
        super().__init__()
        self._grid = UniformGrid(cell_coadd.outer_cell_size, cell_coadd.grid.shape)
        if pad_psfs_with is None:
            self._psf_grid = UniformGrid(cell_coadd.psf_image_size, cell_coadd.grid.shape)
        elif (cell_coadd.psf_image_size.x > cell_coadd.outer_cell_size.x) or (
            cell_coadd.psf_image_size.y > cell_coadd.outer_cell_size.y
        ):
            raise ValueError(
                f"PSF image dimensions {cell_coadd.psf_image_size} are larger than "
                f"outer cell dimensions {cell_coadd.outer_cell_size}; cannot pad."
            )
        else:
            self._psf_grid = self._grid
        self._cell_coadd = cell_coadd
        self._pad_psfs_with = pad_psfs_with
        self._psf_image: ImageF | None = None

    @property
    def bbox(self) -> Box2I:
        # Docstring inherited.
        return self._grid.bbox

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
    def n_noise_realizations(self) -> int:
        # Docstring inherited.
        return self._cell_coadd.n_noise_realizations

    @property
    def mask_fraction_names(self) -> AbstractSet[str]:
        # Docstring inherited.
        return self._cell_coadd.mask_fraction_names

    def _iter_cell_planes(self) -> Iterator[ImagePlanes]:
        # Docstring inherited.
        for cell in self._cell_coadd.cells.values():
            new_bbox = self._grid.bbox_of(cell.identifiers.cell.index)

            def _make_view(original: ImageLike) -> ImageLike:
                result = original[:, :]  # copy bbox, share pixel data.
                result.setXY0(new_bbox.getMin())
                return result

            yield ViewImagePlanes(cell.outer, _make_view, bbox=new_bbox)

    @property
    def psf_image(self) -> ImageF:
        """A stitched-together image of the PSF models for each cell."""
        if self._psf_image is None:
            stitched_psf_image = ImageF(self.psf_grid.bbox)
            if self._pad_psfs_with is not None:
                stitched_psf_image.set(self._pad_psfs_with)
            for cell in self._cell_coadd.cells.values():
                target_subimage = stitched_psf_image[self.psf_grid.bbox_of(cell.identifiers.cell)]
                target_dimensions = target_subimage.getDimensions()
                source_dimensions = cell.psf_image.getDimensions()
                offset = (target_dimensions - source_dimensions) // 2
                # Use numpy views instead of Image methods because the images
                # are not in the same coordinate system to begin with, so xy0
                # doesn't help us.
                target_subimage.array[
                    offset.y : offset.y + source_dimensions.y, offset.x : offset.x + source_dimensions.x
                ] = cell.psf_image.array
            self._psf_image = stitched_psf_image
        return self._psf_image
