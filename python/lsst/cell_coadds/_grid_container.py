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

from collections.abc import Callable, Iterable, Iterator, MutableMapping
from itertools import product
from typing import TYPE_CHECKING, TypeVar

from lsst.skymap import Index2D

from ._uniform_grid import UniformGrid

if TYPE_CHECKING:
    import lsst.geom as geom


__all__ = ("GridContainer",)

T = TypeVar("T")


class GridContainer(MutableMapping[Index2D, T]):
    """A container whose elements form a 2-d grid.

    Parameters
    ----------
    shape : `lsst.skymap.Index2D`
        The number of cells in the grid in each dimension.
    offset : `lsst.skymap.Index2D` or None, optional
        The integer offset of the grid in each dimension. If `None`, the offset
        is Index2D(0, 0).
    """

    def __init__(self, shape: Index2D, offset: Index2D | None = None) -> None:
        super().__init__()

        self._offset = offset if offset else Index2D(0, 0)
        self._shape = shape
        self._cells: dict[Index2D, T] = {}
        self._mapping = self._cells

    def _check(self, index: Index2D) -> None:
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

    def __getitem__(self, index: Index2D) -> T:
        return self._cells[index]

    def __setitem__(self, index: Index2D, value: T) -> None:
        self._check(index)
        self._cells[index] = value

    def __delitem__(self, index: Index2D) -> None:
        self._check(index)
        del self._cells[index]

    def __iter__(self) -> Iterator[T]:
        return iter(self._cells)

    # def keys(self) -> Iterable[Index2D]:  # populated_indices
    #     """Return an iterator over the indices in the container with values.

    #     See also
    #     --------
    #     indices
    #     """
    #     yield from self._cells.keys()

    def indices(self) -> Iterable[Index2D]:
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
            yield Index2D(x, y)

    @property
    def shape(self) -> Index2D:
        """Number of cells in the container in each dimension."""
        return self._shape

    @property
    def offset(self) -> Index2D:
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
        return self._cells[Index2D(x=self.offset.x + self.shape.x - 1, y=self.offset.y + self.shape.y - 1)]

    def subset_overlapping(self, grid: UniformGrid, bbox: geom.Box2I) -> GridContainer:
        """Return a new GridContainer with cells that overlap a bounding box.

        Parameters
        ----------
        grid : `~lsst.cell_coadds.UniformGrid`
            Grid that maps the container's cells to the coordinates used to
            define the bounding box.  May define a grid that is a super of the
            container's cells.
        bbox : `~lsst.geom.Box2I`
            Bounding box that returned cells must overlap.

        Returns
        -------
        grid_container : `GridContainer`
            GridContainer with just the cells that overlap the bounding box.
        """
        bbox = bbox.clippedTo(grid.bbox)
        offset = grid.index(bbox.getBegin())
        last = grid.index(bbox.getMax())
        end = Index2D(last.x + 1, last.y + 1)
        shape = Index2D(end.x - offset.x, end.y - offset.y)
        gc = GridContainer[T](shape, offset)
        for index in gc.indices():
            gc[index] = self[index]
        return gc

    def rebuild_transformed(self, transform: Callable[[T], T]) -> GridContainer:
        """Return a GridContainer with the same shape and offset.

        The cell values are created by applying a callback function to each
        cell value in this object.

        Parameters
        ----------
        transform: Callable[[T], T]
            A callable function that takes a cell value and returns a new
        """
        gc = GridContainer[T](self.shape, self.offset)
        for key in self.keys():  # noqa: SIM118
            gc[key] = transform(self[key])
        return gc
