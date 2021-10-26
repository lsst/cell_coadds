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

from typing import Any, Iterable, Optional, Tuple

import numpy as np
from lsst.geom import Box2I

from ._single_cell_coadd import SingleCellCoadd
from ._stitched import StitchedCellCoadd


class MultipleCellCoadd:
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

    def __init__(self, cells: Iterable[SingleCellCoadd], *, inner_bbox: Optional[Box2I] = None):
        cells = list(cells)
        x_starts_set = {cell.outer.bbox.getBeginX() for cell in cells}
        y_starts_set = {cell.outer.bbox.getBeginY() for cell in cells}
        if len(x_starts_set) * len(y_starts_set) != len(cells):
            raise ValueError("Irregular or incomplete cell grid.")
        self._x_starts = np.array(sorted(x_starts_set), dtype=np.int64)
        self._y_starts = np.array(sorted(y_starts_set), dtype=np.int64)
        self._cells = np.empty((len(self._y_starts), len(self._x_starts)), dtype=object)
        self._cells.flags.writeable = False  # don't allow cells to be reassigned.
        self._n_noise_realizations = None
        for cell in cells:
            x_index = self._x_starts[cell.outer.bbox.getBeginX()]
            y_index = self._y_starts[cell.outer.bbox.getBeginY()]
            self._cells[y_index, x_index] = cell
            if self._n_noise_realizations is None:
                self._n_noise_realizations = len(cell.outer.noise_realizations)
            else:
                if self._n_noise_realizations != len(cell.outer.noise_realizations):
                    raise ValueError("Inconsistent number of noise realizations between cells.")
        # TODO: check for gaps
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

    def __getitem__(self, bbox: Box2I) -> MultipleCellCoadd:
        # TODO: test for off-by-one errors here, especially for ends
        x_index_begin = np.searchsorted(self._x_starts, bbox.getBeginX())
        x_index_end = np.searchsorted(self._x_starts, bbox.getEndX())
        y_index_begin = np.searchsorted(self._y_starts, bbox.getBeginY())
        y_index_end = np.searchsorted(self._y_starts, bbox.getEndY())
        result = MultipleCellCoadd.__new__()  # type: ignore
        result._x_starts = self._x_starts[x_index_begin:x_index_end].copy()
        result._y_starts = self._y_starts[y_index_begin:y_index_end].copy()
        result.cells = self.cells[y_index_begin:y_index_end, x_index_begin:x_index_end].copy()
        result._inner_bbox = bbox
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
