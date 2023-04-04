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
        The size of each grid cell.
    shape : `lsst.skymap.Index2D`
        The number of cells in the grid in each dimension.
    min : `lsst.geom.Point2I` or None, optional
        The minimum (lower left) corner of the grid. If `None`, the minimum
        corner is set to be (0, 0).
    """

    def __init__(self, cell_size: geom.Extent2I, shape: Index2D, min: geom.Point2I | None = None) -> None:
        self._cell_size = cell_size
        self._shape = shape
        if min is None:
            min = geom.Point2I(0, 0)
        self._bbox = geom.Box2I(min, geom.Extent2I(cell_size.getX() * shape.x, cell_size.getY() * shape.y))

    # Factory methods for constructing a UniformGrid
    @classmethod
    def from_bbox_shape(cls, bbox: geom.Box2I, shape: Index2D) -> UniformGrid:
        """Factory method to construct from a bounding box and a shape.

        Parameters
        ----------
        bbox : `lsst.geom.Box2I`
            Bounding box of the full grid.
        shape : `lsst.skymap.Index2D`
            Number of cells in the grid in each dimension.
            Must divide the ``bbox`` width and height evenly.

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
        return cls(cell_size, shape, bbox.getMin())

    @classmethod
    def from_bbox_cell_size(cls, bbox: geom.Box2I, cell_size: geom.Extent2I) -> UniformGrid:
        """Factor method to construct from a bounding box and a cell size.

        Parameters
        ----------
        bbox : `lsst.geom.Box2I`
            Bounding box of the full grid.
        cell_size : `lsst.geom.Extent2I`
            Size of each grid cell.
            Must divide the ``bbox`` width and height evenly.

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
        return cls(cell_size, shape, bbox.getMin())

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
    def cell_size(self) -> geom.Extent2I:
        return self._cell_size

    @property
    def shape(self) -> Index2D:
        return self._shape

    # Implement C++ like getters
    def get_bbox(self) -> geom.Box2I:
        return self._bbox

    def get_cell_size(self) -> geom.Extent2I:
        return self._cell_size

    def get_shape(self) -> Index2D:
        return self._shape

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
            Raised if ``position`` is not within the grid's bounding box.
        """
        if not self.bbox.contains(position):
            raise ValueError(f"Position {position} is not within bounding box {self.bbox}.s")

        offset = position - self.bbox.getBegin()
        return Index2D(offset.x // self._cell_size.x, offset.y // self._cell_size.y)

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
        # TODO: Need a check to see if the index is a valid one.
        return geom.Point2I(
            index.x * self.cell_size.x + self.bbox.getBeginX(),
            index.y * self.cell_size.y + self.bbox.getBeginY(),
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
        return geom.Box2I(self.min_of(index), self.cell_size)
