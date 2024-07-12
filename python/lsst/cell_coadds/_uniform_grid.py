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

__all__ = ("UniformGrid",)

from typing import Self

import lsst.shoefits as shf
import pydantic

from ._to_upstream import CellIndex, CellShape, PixelIndex, PixelShape


class UniformGrid(pydantic.BaseModel):
    """A 2-dimensional integer grid.

    Parameters
    ----------
    cell_size : `PixelShape`
        The size of each interior grid cell.
    shape : `lsst.skymap.Index2D`
        The number of cells in the grid in each dimension.
    padding : `int`, optional
        The number of pixels to pad the grid in each dimension.  This increases
        the cell size of the first and last rows and columns of cells.
    min : `lsst.geom.Point2I` or None, optional
        The minimum (lower left) corner of the interior grid, excluding
        ``padding``. If `None`, the minimum corner is set to be (0, 0).
    """

    cell_size: PixelShape
    """The size of each interior grid cell in pixels."""

    shape: CellShape
    """The number of cells in each dimension."""

    padding: int = pydantic.Field(ge=0)
    """Number of pixels to pad the outermost cells with on all sides."""

    bbox: shf.Box
    """Bounding box of the full grid, not including `padding`."""

    # Factory methods for constructing a UniformGrid
    @classmethod
    def from_bbox_shape(cls, bbox: shf.Box, shape: CellShape, padding: int = 0) -> Self:
        """Generate a UniformGrid instance from a bounding box and a shape.

        Parameters
        ----------
        bbox : `lsst.shoefits.Box`
            Bounding box of the full grid, not including ``padding``.
        shape : `CellShape`
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
        if not all(b.size % s == 0 for b, s in zip(bbox, shape, strict=True)):
            raise IndexError(f"Sizes of bounding box {bbox} are not evenly divided by shape {shape}.")
        cell_size = PixelShape.from_yx([b.size // s for b, s in zip(bbox, shape, strict=True)])
        return cls(
            cell_size=cell_size,
            shape=shape,
            padding=padding,
            bbox=bbox,
        )

    @classmethod
    def from_bbox_cell_size(cls, bbox: shf.Box, cell_size: PixelShape, padding: int = 0) -> Self:
        """Generate a UniformGrid instance from a bounding box and a cell size.

        Parameters
        ----------
        bbox : `lsst.shoefits.Box`
            Bounding box of the full grid, not including ``padding``.
        cell_size : `PixelShape`
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
        if not all(b.size % c == 0 for b, c in zip(bbox, cell_size, strict=True)):
            raise IndexError(f"Sizes of bounding box {bbox} are not evenly divided by shape {cell_size}.")
        shape = CellShape.from_yx([b.size // c for b, c in zip(bbox, cell_size, strict=True)])
        return cls(
            cell_size=cell_size,
            shape=shape,
            padding=padding,
            bbox=bbox,
        )

    @classmethod
    def from_cell_size_shape(
        cls, cell_size: PixelShape, shape: CellShape, padding: int = 0, min: PixelIndex | None = None
    ) -> Self:
        """Generate a UniformGrid instance from a bounding box and a cell size.

        Parameters
        ----------
        cell_size : `PixelShape`
            Size of each interior grid cell.
        shape : `CellShape`
            Number of cells in the grid in each dimension.
        padding : `int`, optional
            The number of pixels to pad the grid in each dimension.
        min : `PixelIndex`, optional
            Minimum point of bbox (not including padding).  Defaults to
            ``(0, 0)``.

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
        if min is None:
            min = PixelIndex(x=0, y=0)
        return cls(
            cell_size=cell_size,
            shape=shape,
            bbox=shf.Box(
                *[
                    shf.Interval.from_size(c * s, start=m)
                    for c, s, m in zip(cell_size, shape, min, strict=True)
                ]
            ),
            padding=padding,
        )

    @pydantic.model_validator(mode="after")
    def _validate(self) -> Self:
        if any(b.size != s * c for b, s, c in zip(self.bbox, self.shape, self.cell_size, strict=True)):
            raise ValueError(f"bbox {self.bbox}, shape {self.shape}, and {self.cell_size} are not consistent")
        return self

    @property
    def bbox_with_padding(self) -> shf.Box:
        return self.bbox.dilated_by(self.padding)

    # Convenience methods
    def index(self, position: PixelIndex) -> CellIndex:
        """Index of the cell that contains the given point.

        Parameters
        ----------
        position : `PixelIndex`
            A pixel-coordinate point in the grid.

        Returns
        -------
        index : `CellIndex`
            A 2D index of the cell containing ``position``.

        Raises
        ------
        ValueError
            Raised if ``position`` is not within the grid's bounding box
            including the padding.
        """
        if position not in self.bbox_with_padding:
            raise ValueError(
                f"Position {position} is not within outer bounding box {self.bbox_with_padding}.s"
            )
        result = []
        for n in range(2):
            offset = position[n] - self.bbox[n].start
            if offset < 0:
                result.append(0)
            elif offset >= self.shape[n] * self.cell_size[n]:
                result.append(self.shape[n] - 1)
            else:
                result.append(offset // self.cell_size[n])
        return CellIndex.from_yx(result)

    def min_of(self, index: CellIndex) -> PixelIndex:
        """Minimum point of a single cell's bounding box.

        Parameters
        ----------
        index : `CellIndex`
            A 2D index of the cell.

        Returns
        -------
        point : `PixelIndex`
            The minimum point of the cell's bounding box.
        """
        if not all(0 <= i < s for i, s in zip(index, self.shape, strict=True)):
            raise ValueError(f"{index} is not within the grid's shape {self.shape}.")
        return PixelIndex.from_yx(
            [
                i * c + b.start - (self.padding if i == 0 else 0)
                for i, c, b in zip(index, self.cell_size, self.bbox, strict=True)
            ]
        )

    def bbox_of(self, index: CellIndex) -> shf.Box:
        """Bounding box of the cell at the given index.

        Parameters
        ----------
        index : `CellIndex`
            A 2D index of the cell.

        Returns
        -------
        bbox : `lsst.geom.Box2I`
            The bounding box of the cell.
        """
        cell_size = tuple(
            [c + (self.padding if i == 0 else 0) for i, c in zip(index, self.cell_size, strict=True)]
        )
        return shf.Box.from_shape(cell_size, start=self.min_of(index))
