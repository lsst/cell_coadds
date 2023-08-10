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

import lsst.geom as geom
from lsst.skymap import Index2D

__all__ = ("UniformGrid",)


class UniformGrid:
    """A 2-dimensional integer grid.

    Parameters
    ----------
    cell_size : `lsst.geom.Extent2I`
        The size of each interior grid cell.
    shape : `lsst.skymap.Index2D`
        The number of cells in the grid in each dimension.
    padding : `int`, optional
        The number of pixels to pad the grid in each dimension.
    min : `lsst.geom.Point2I` or None, optional
        The minimum (lower left) corner of the interior grid, excluding
        ``padding``. If `None`, the minimum corner is set to be (0, 0).
    """

    def __init__(
        self, cell_size: geom.Extent2I, shape: Index2D, *, padding: int = 0, min: geom.Point2I | None = None
    ) -> None:
        self._cell_size = cell_size
        self._shape = shape
        self._padding = padding
        if min is None:
            min = geom.Point2I(0, 0)
        self._bbox = geom.Box2I(
            geom.Point2I(min.x, min.y),
            geom.Extent2I(cell_size.getX() * shape.x, cell_size.getY() * shape.y),
        )

    # Factory methods for constructing a UniformGrid
    @classmethod
    def from_bbox_shape(cls, bbox: geom.Box2I, shape: Index2D, padding: int = 0) -> UniformGrid:
        """Factory method to construct from a bounding box and a shape.

        Parameters
        ----------
        bbox : `lsst.geom.Box2I`
            Bounding box of the full grid (without including ``padding``).
        shape : `lsst.skymap.Index2D`
            Number of cells in the grid in each dimension.
            Must divide the ``bbox`` width and height evenly.
        padding : `int`, optional
            The number of pixels to pad the grid in each dimension.

        Returns
        -------
        grid : `UniformGrid`
            A new UniformGrid instance.

        Raises
        ------
        LengthError
            Raised if ``shape`` dimensions do not divide the ``bbox``
            dimensions evenly.
        """
        cls._validate_bbox_shape(bbox, shape)

        cell_size = geom.Extent2I(bbox.getWidth() // shape.x, bbox.getHeight() // shape.y)
        return cls(cell_size, shape, min=bbox.getMin(), padding=padding)

    @classmethod
    def from_bbox_cell_size(cls, bbox: geom.Box2I, cell_size: geom.Extent2I, padding: int = 0) -> UniformGrid:
        """Factor method to construct from a bounding box and a cell size.

        Parameters
        ----------
        bbox : `lsst.geom.Box2I`
            Bounding box of the full grid (without including ``padding``).
        cell_size : `lsst.geom.Extent2I`
            Size of each interior grid cell.
            Must divide the ``bbox`` width and height evenly.
        padding : `int`, optional
            The number of pixels to pad the grid in each dimension.

        Returns
        -------
        grid : `UniformGrid`
            A new UniformGrid instance.

        Raises
        ------
        IndexError
            Raised if ``cell_size`` dimensions do not divide the ``bbox``
            dimensions evenly.
        """
        cls._validate_bbox_cell_size(bbox, cell_size)
        shape = Index2D(bbox.getWidth() // cell_size.x, bbox.getHeight() // cell_size.y)
        return cls(cell_size, shape, padding=padding, min=bbox.getMin())

    # Methods to validate the input parameters
    @staticmethod
    def _validate_bbox_shape(bbox: geom.Box2I, shape: Index2D) -> None:
        if bbox.getWidth() % shape.x != 0:
            raise IndexError(
                f"Bounding box width {bbox.getWidth()} is not evenly divided by x shape {shape.x}."
            )
        if bbox.getHeight() % shape.y != 0:
            raise IndexError(
                f"Bounding box height {bbox.getHeight()} is not evenly divided by y shape {shape.y}."
            )

    @staticmethod
    def _validate_bbox_cell_size(bbox: geom.Box2I, cell_size: geom.Extent2I) -> None:
        if bbox.getWidth() % cell_size.x != 0:
            raise IndexError(
                f"Bounding box width {bbox.getWidth()} is not evenly divided by x cell_size "
                f"{cell_size.getX()}."
            )

        if bbox.getHeight() % cell_size.y != 0:
            raise IndexError(
                f"Bounding box height {bbox.getHeight()} is not evenly divided by y cell_size "
                f"{cell_size.getY()}."
            )

    # Pythonic property getters
    @property
    def bbox(self) -> geom.Box2I:
        return self._bbox

    @property
    def bbox_with_padding(self) -> geom.Box2I:
        return self._bbox.dilatedBy(self._padding)

    @property
    def cell_size(self) -> geom.Extent2I:
        return self._cell_size

    @property
    def shape(self) -> Index2D:
        return self._shape

    @property
    def padding(self) -> int:
        return self._padding

    # Implement C++ like getters
    def get_bbox(self) -> geom.Box2I:
        return self._bbox

    def get_bbox_with_padding(self) -> geom.Box2I:
        return self.bbox_with_padding

    def get_cell_size(self) -> geom.Extent2I:
        return self._cell_size

    def get_shape(self) -> Index2D:
        return self._shape

    def get_padding(self) -> int:
        return self._padding

    # Dunder methods
    def __eq__(self, other: object) -> bool:
        if not isinstance(other, UniformGrid):
            return False
        return self._bbox == other._bbox and self._cell_size == other._cell_size

    def __repr__(self) -> str:
        return (
            f"UniformGrid(cell_size={repr(self.cell_size)}, shape={self.shape}, "
            f"min={repr(self.bbox.getMin())})"
        )

    # Convenience methods
    def index(self, position: geom.Point2I) -> Index2D:
        """Index of the cell that contains the given point.

        Parameters
        ----------
        position : `lsst.geom.Point2I`
            A point in the grid.

        Returns
        -------
        index : `lsst.skymap.Index2D`
            A 2D index of the cell containing ``position``.

        Raises
        ------
        ValueError
            Raised if ``position`` is not within the grid's bounding box
            including the padding.
        """
        if not self.bbox_with_padding.contains(position):
            raise ValueError(
                f"Position {position} is not within outer bounding box {self.bbox_with_padding}.s"
            )

        offset = position - self.bbox.getBegin()

        if offset.x < 0:
            x = 0
        elif offset.x >= self.shape.x * self.cell_size.x:
            x = self.shape.x - 1
        else:
            x = offset.x // self.cell_size.x

        if offset.y < 0:
            y = 0
        elif offset.y >= self.shape.y * self.cell_size.y:
            y = self.shape.y - 1
        else:
            y = offset.y // self.cell_size.y

        return Index2D(x, y)

    def min_of(self, index: Index2D) -> geom.Point2I:
        """Minimum point of a single cell's bounding box.

        Parameters
        ----------
        index : `~lsst.skymap.Index2D`
            A 2D index of the cell.

        Returns
        -------
        point : `lsst.geom.Point2I`
            The minimum point of the cell's bounding box.
        """
        if not (0 <= index.x < self._shape.x and 0 <= index.y < self._shape.y):
            raise ValueError(f"{index} is not within the grid's shape {self._shape}.")

        offset = geom.Point2I(
            -self._padding if index.x == 0 else 0,
            -self._padding if index.y == 0 else 0,
        )
        return geom.Point2I(
            index.x * self.cell_size.x + self.bbox.getBeginX() + offset.x,
            index.y * self.cell_size.y + self.bbox.getBeginY() + offset.y,
        )

    def bbox_of(self, index: Index2D) -> geom.Box2I:
        """Bounding box of the cell at the given index.

        Parameters
        ----------
        index : `~lsst.skymap.Index2D`
            A 2D index of the cell.

        Returns
        -------
        bbox : `lsst.geom.Box2I`
            The bounding box of the cell.
        """
        # Compute the buffer to add if ``index`` corresponds to the leftmost or
        # the rightmost cell or the topmost or the bottommost cell.
        buffer = geom.Extent2I(
            self.padding if index.x in {0, self.shape.x - 1} else 0,
            self.padding if index.y in {0, self.shape.y - 1} else 0,
        )
        return geom.Box2I(self.min_of(index), self.cell_size + buffer)
