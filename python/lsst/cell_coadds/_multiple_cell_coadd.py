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

from typing import AbstractSet, Any, Iterable, Optional, Set, Tuple

import numpy as np
from lsst.geom import Box2I

from ._common_components import CommonComponents, CommonComponentsProperties
from ._single_cell_coadd import SingleCellCoadd
from ._stitched import StitchedCellCoadd

# Need stubs for compiled modules.
from ._cell_coadds import UniformGrid  # type: ignore


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
        *,
        common: CommonComponents,
        inner_bbox: Optional[Box2I] = None,
    ):
        self._grid = grid
        self._common = common
        self._cells = np.empty(self._grid.shape, dtype=object)
        self._cells.flags.writeable = False  # don't allow cells to be reassigned.
        self._n_noise_realizations = None
        self._mask_fraction_names: Set[str] = set()
        for cell in cells:
            index = self._grid.index(cell.inner.bbox.getBegin())
            self._cells[index] = cell
            if cell.inner.bbox != self._grid.bbox_of(index):
                raise ValueError(
                    f"Cell at index row={index[0]}, col={index[1]} has "
                    f"inner bbox {cell.inner.bbox}, but grid expects {self._grid.bbox_of(index)}."
                )
            if self._n_noise_realizations is None:
                self._n_noise_realizations = len(cell.outer.noise_realizations)
            else:
                if self._n_noise_realizations != len(cell.outer.noise_realizations):
                    raise ValueError("Inconsistent number of noise realizations between cells.")
            self._mask_fraction_names.update(cell.outer.mask_fractions.keys())
        max_inner_bbox = Box2I(self._cells[0, 0].inner.bbox.getMin(), self._cells[-1, -1].inner.bbox.getMax())
        if inner_bbox is None:
            inner_bbox = max_inner_bbox
        elif not max_inner_bbox.contains(inner_bbox):
            raise ValueError(
                f"Requested inner bounding box {inner_bbox} is not fully covered by these "
                f"cells (bbox is {max_inner_bbox}."
            )
        self._inner_bbox = inner_bbox

    # TODO: remove "type: ignore" below after we adopt np 1.21 and turn on its
    # mypy plugin.
    @property
    def cells(self) -> np.ndarray[Tuple[Any, Any], SingleCellCoadd]:  # type: ignore
        """The grid of single-cell coadds, indexed by (y, x)."""
        return self._cells

    @property
    def n_noise_realizations(self) -> int:
        """The number of noise realizations cells are guaranteed to have."""
        if self._n_noise_realizations is None:
            return 0
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

    def __getitem__(self, bbox: Box2I) -> MultipleCellCoadd:
        min_index = self._grid.index(bbox.getMin())
        max_index = self._grid.index(bbox.getMax())
        result = MultipleCellCoadd.__new__()  # type: ignore
        result._grid = self._grid
        result._common = self._common
        result._cells = self._cells[min_index[0] : max_index[0] + 1, min_index[1] : max_index[1] + 1].copy()
        result._inner_bbox = bbox
        result._n_noise_realizations = self._n_noise_realizations
        result._mask_fraction_names = self._mask_fraction_names
        return result

    @property
    def outer_bbox(self) -> Box2I:
        """The rectangular region fully covered by all cell outer bounding
        boxes."""
        return Box2I(self.cells[0, 0].outer.bbox.getMin(), self.cells[-1, -1].outer.bbox.getMax())

    @property
    def inner_bbox(self) -> Box2I:
        """The rectangular region fully covered by all cell inner bounding
        boxes."""
        return self._inner_bbox

    @property
    def common(self) -> CommonComponents:
        # Docstring inherited.
        return self._common

    def stitch(self, bbox: Optional[Box2I] = None) -> StitchedCellCoadd:
        """Return a contiguous (but in general discontinuous) coadd by
        stitching together inner cells.

        Parameters
        ----------
        bbox : `Box2I`, optional
            Region for the returned coadd; default is ``self.inner_bbox``.

        Returns
        -------
        stitched : `StichedCellCoadd`
            Contiguous coadd covering the given area.  Each image plane is
            actually constructed when first accessed, not when this method
            is called.
        """
        # In the future, stitching algorithms that apply ramps to smooth
        # discontinuities may also be provided; we'd implement that by having
        # this return different types (from a common ABC), perhaps dispatched
        # by an enum.
        return StitchedCellCoadd(self, bbox=bbox)
