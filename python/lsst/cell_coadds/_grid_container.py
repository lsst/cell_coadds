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

__all__ = ("GridContainer",)

from collections.abc import Callable, Iterable, Iterator, MutableMapping
from itertools import product
from typing import TypeVar

import lsst.shoefits as shf

from ._to_upstream import CellIndex, CellShape, PixelIndex
from ._uniform_grid import UniformGrid

T = TypeVar("T")
U = TypeVar("U")


class GridContainer(MutableMapping[CellIndex, T]):
    """A container whose elements form a 2-d grid.

    Parameters
    ----------
    shape : `~lsst.skymap.Index2D`
        The number of cells in the grid in each dimension.
    offset : `~lsst.skymap.Index2D` or None, optional
        The integer offset of the grid in each dimension. If `None`, the offset
        is ``CellIndex(y=0, x=0)``.
    """

    def __init__(self, shape: CellShape, offset: CellIndex | None = None):
        super().__init__()

        self._offset = offset if offset else CellIndex(x=0, y=0)
        self._shape = shape
        self._cells: dict[CellIndex, T] = {}
        self._mapping = self._cells

    def _check(self, index: CellIndex) -> None:
        """Check if a given index belongs to the container or not."""
        if index.x < self.offset.x or index.x >= self.offset.x + self.shape.x:
            raise ValueError(
                f"x index {index.x} out of range; expected a value between {self.offset.x} and "
                f"{self.offset.x + self.shape.x -1}."
            )
        if index.y < self.offset.y or index.y >= self.offset.y + self.shape.y:
            raise ValueError(
                f"y index {index.y} out of range; expected a value between {self.offset.y} and "
                f"{self.offset.y + self.shape.y -1}."
            )

    def __repr__(self) -> str:
        return f"{type(self).__name__}(shape={self.shape}, offset={self.offset})"

    # Implement the necessary methods for MutableMapping.
    def __len__(self) -> int:
        return len(self._cells)

    def __getitem__(self, index: CellIndex) -> T:
        return self._cells[index]

    def __setitem__(self, index: CellIndex, value: T) -> None:
        self._check(index)
        self._cells[index] = value

    def __delitem__(self, index: CellIndex) -> None:
        self._check(index)
        del self._cells[index]

    def __iter__(self) -> Iterator[CellIndex]:
        return iter(self._cells)

    def indices(self) -> Iterable[CellIndex]:
        """Return an iterator over all possible indices for the container.

        Unlike `keys`, this method returns an iterator over all valid indices,
        whether the corresponding value is set or not.

        See Also
        --------
        keys
        """
        iterator = product(
            range(self.offset.y, self.offset.y + self.shape.y),
            range(self.offset.x, self.offset.x + self.shape.x),
        )
        for y, x in iterator:
            yield CellIndex(x=x, y=y)

    @property
    def shape(self) -> CellShape:
        """Number of cells in the container in each dimension."""
        return self._shape

    @property
    def offset(self) -> CellIndex:
        """Index of the first cell in the container."""
        return self._offset

    @property
    def size(self) -> int:
        """The number of cells expected in the container.

        This does not indicate the number of cells that have been filled.
        Use len() instead.
        """
        return self._shape.x * self._shape.y

    @property
    def first(self) -> T:
        """The cell at the lower left corner of the container."""
        return self._cells[self.offset]

    @property
    def last(self) -> T:
        """The cell at the upper right corner of the container."""
        return self._cells[CellIndex(x=self.offset.x + self.shape.x - 1, y=self.offset.y + self.shape.y - 1)]

    def subset_overlapping(self, grid: UniformGrid, bbox: shf.Box) -> GridContainer:
        """Return a new GridContainer with cells that overlap a bounding box.

        Parameters
        ----------
        grid : `UniformGrid`
            Grid that maps the container's cells to the coordinates used to
            define the bounding box.  May define a grid that is a superset of
            the container's cells.
        bbox : `lsst.shoefits.Box`
            Bounding box that returned cells must overlap.

        Returns
        -------
        grid_container : `GridContainer`
            GridContainer with just the cells that overlap the bounding box.
        """
        if (clipped_bbox := bbox.intersection(grid.bbox)) is None:
            return GridContainer(CellShape(x=0, y=0))
        offset = grid.index(PixelIndex(x=clipped_bbox.x.min, y=clipped_bbox.y.min))
        last = grid.index(PixelIndex(x=clipped_bbox.x.max, y=clipped_bbox.y.max))
        shape = CellShape(x=last.x + 1 - offset.x, y=last.y + 1 - offset.y)
        gc = GridContainer[T](shape, offset)
        for index in gc.indices():
            gc[index] = self[index]
        return gc

    def rebuild_transformed(self, transform: Callable[[T], U]) -> GridContainer[U]:
        """Return a GridContainer with the same shape and offset.

        The cell values are created by applying a callback function to each
        cell value in this object.

        Parameters
        ----------
        transform : Callable[[T], T]
            A callable function that takes a cell value and returns a new
        """
        gc = GridContainer[U](self.shape, self.offset)
        for key, value in self.items():
            gc[key] = transform(value)
        return gc
