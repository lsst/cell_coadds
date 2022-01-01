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

__all__ = ("MultipleCellCoadd",)

from typing import AbstractSet, Iterable, Optional, Set

from lsst.geom import Box2I, Extent2I

from ._cell_coadds import GridContainer, GridContainerBuilder, UniformGrid
from ._common_components import CommonComponents, CommonComponentsProperties
from ._exploded_coadd import ExplodedCoadd
from ._single_cell_coadd import SingleCellCoadd
from ._stitched_coadd import StitchedCoadd


class MultipleCellCoadd(CommonComponentsProperties):
    """A data structure for coadds built from many overlapping cells.

    Notes
    -----
    `MultipleCellCoadd` is designed to be used both by measurement algorithms
    that are able to take advantage of cell boundaries and overlap regions
    (which can use the ``.cells`` attribute to access `SingleCellCoadd` objects
    directly) and measurement algorithms that just want one image and don't
    care (or don't care much) about discontinuities (which can use `stitch` to
    obtain such an image).

    Indexing with `Box2I` yields a `MultipleCellCoadd` view containing just the
    cells that overlap that region.
    """

    def __init__(
        self,
        cells: Iterable[SingleCellCoadd],
        grid: UniformGrid,
        outer_cell_size: Extent2I,
        psf_image_size: Extent2I,
        *,
        common: CommonComponents,
        inner_bbox: Optional[Box2I] = None,
    ):
        self._grid = grid
        self._outer_cell_size = outer_cell_size
        self._psf_image_size = psf_image_size
        self._common = common
        cells_builder: GridContainerBuilder[SingleCellCoadd] = GridContainerBuilder(self._grid.shape)
        self._mask_fraction_names: Set[str] = set()
        for cell in cells:
            index = cell.identifiers.cell.index
            cells_builder[index] = cell
            if cell.inner.bbox != self._grid.bbox_of(index):
                raise ValueError(
                    f"Cell at index {index} has inner bbox {cell.inner.bbox}, "
                    f"but grid expects {self._grid.bbox_of(index)}."
                )
            if cell.outer.bbox.getDimensions() != self._outer_cell_size:
                raise ValueError(
                    f"Cell at index {index} has outer dimensions {cell.outer.bbox.getDimensions()}, "
                    f"but coadd expects {self._outer_cell_size}."
                )
            if cell.psf_image.getDimensions() != self._psf_image_size:
                raise ValueError(
                    f"Cell at index {index} has PSF image with dimensions {cell.psf_image.getDimensions}, "
                    f"but coadd expects {self._psf_image_size}."
                )
            self._mask_fraction_names.update(cell.outer.mask_fractions.keys())
        self._cells = cells_builder.finish()
        n_noise_realizations = {len(cell.outer.noise_realizations) for cell in self._cells}
        self._n_noise_realizations = n_noise_realizations.pop()
        if n_noise_realizations:
            n_noise_realizations.add(self._n_noise_realizations)
            raise ValueError(
                f"Inconsistent number of noise realizations ({n_noise_realizations}) betwen cells."
            )
        max_inner_bbox = Box2I(self._cells.first.inner.bbox.getMin(), self._cells.last.inner.bbox.getMax())
        if inner_bbox is None:
            inner_bbox = max_inner_bbox
        elif not max_inner_bbox.contains(inner_bbox):
            raise ValueError(
                f"Requested inner bounding box {inner_bbox} is not fully covered by these "
                f"cells (bbox is {max_inner_bbox}."
            )
        self._inner_bbox = inner_bbox

    @property
    def cells(self) -> GridContainer[SingleCellCoadd]:
        """The grid of single-cell coadds, indexed by (y, x)."""
        return self._cells

    @property
    def n_noise_realizations(self) -> int:
        """The number of noise realizations cells are guaranteed to have."""
        return self._n_noise_realizations

    @property
    def mask_fraction_names(self) -> AbstractSet[str]:
        """The names of all mask planes whose fractions were propagated in any
        cell.

        Cells that do not have a mask fraction for a particular name may be
        assumed to have the fraction for that mask plane uniformly zero.
        """
        return self._mask_fraction_names

    @property
    def grid(self) -> UniformGrid:
        """Object that defines the inner geometry for all cells."""
        return self._grid

    @property
    def outer_cell_size(self) -> Extent2I:
        """Dimensions of the outer region of each cell."""
        return self._outer_cell_size

    @property
    def psf_image_size(self) -> Extent2I:
        """Dimensions of PSF model images."""
        return self._psf_image_size

    @property
    def outer_bbox(self) -> Box2I:
        """The rectangular region fully covered by all cell outer bounding
        boxes."""
        return Box2I(self.cells.first.outer.bbox.getMin(), self.cells.last.outer.bbox.getMax())

    @property
    def inner_bbox(self) -> Box2I:
        """The rectangular region fully covered by all cell inner bounding
        boxes."""
        return self._inner_bbox

    @property
    def common(self) -> CommonComponents:
        # Docstring inherited.
        return self._common

    def stitch(self, bbox: Optional[Box2I] = None) -> StitchedCoadd:
        """Return a contiguous (but in general discontinuous) coadd by
        stitching together inner cells.

        Parameters
        ----------
        bbox : `Box2I`, optional
            Region for the returned coadd; default is ``self.inner_bbox``.

        Returns
        -------
        stitched : `StitchedCellCoadd`
            Contiguous coadd covering the given area.  Each image plane is
            actually constructed when first accessed, not when this method
            is called.
        """
        # In the future, stitching algorithms that apply ramps to smooth
        # discontinuities may also be provided; we'd implement that by having
        # this return different types (from a common ABC), perhaps dispatched
        # by an enum.
        return StitchedCoadd(self, bbox=bbox)

    def explode(self, pad_psfs_with: Optional[float] = None) -> ExplodedCoadd:
        """Return a coadd whose image planes stitch together the outer regions
        of each cell, duplicating pixels in the overlap regions.

        Parameters
        ----------
        pad_psfs_with : `float` or `None`, optional
            A floating-point value to pad PSF images with so each PSF-image
            cell has the same dimensions as the image (outer) cell it
            corresponds to.  If `None`, PSF images will not be padded and the
            full PSF image will generally be smaller than the exploded image it
            corresponds to.

        Returns
        -------
        exploded : `ExplodedCoadd`
            Exploded version of the coadd.
        """
        return ExplodedCoadd(self, pad_psfs_with=pad_psfs_with)
