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

import pickle
import unittest

import lsst.utils.tests
from lsst.cell_coadds import UniformGrid
from lsst.geom import Box2I, Extent2I, Point2I
from lsst.skymap import Index2D


class UniformGridTestCase(unittest.TestCase):
    """Tests for UniformGrid and GridIndex/Index2D's C++/Python
    translation.
    """

    def setUp(self) -> None:  # noqa: D102
        self.x0 = 1
        self.y0 = 2
        self.bw = 15
        self.bh = 12
        self.nx = 5
        self.ny = 6
        self.cw = 3
        self.ch = 2

        self.bbox = Box2I(Point2I(x=self.x0, y=self.y0), Extent2I(x=self.bw, y=self.bh))
        self.cell_size = Extent2I(x=self.cw, y=self.ch)
        self.shape = Index2D(x=self.nx, y=self.ny)
        self.padding = 0

    def test_ctor_bbox_cell_size(self) -> None:
        """Test UniformGrid after construction with (bbox, cell_size)."""
        grid = UniformGrid.from_bbox_cell_size(self.bbox, self.cell_size)
        self._check(grid)

    def test_ctor_bbox_shape(self) -> None:
        """Test UniformGrid after construction with (bbox, shape)."""
        grid = UniformGrid.from_bbox_shape(self.bbox, self.shape)
        self._check(grid)

    def test_ctor_cell_size_shape_min(self) -> None:
        """Test UniformGrid after construction with (cell_size, shape, min)."""
        grid = UniformGrid(self.cell_size, self.shape, min=self.bbox.getMin())
        self._check(grid)

    def _check(self, grid: UniformGrid) -> None:
        self.assertEqual(grid.bbox, self.bbox)
        self.assertEqual(grid.cell_size, self.cell_size)
        self.assertEqual(grid.shape, self.shape)
        self.assertEqual(grid.padding, self.padding)
        self.assertIsInstance(grid.shape, Index2D)
        self.assertEqual(grid.bbox_of(Index2D(x=3, y=4)), Box2I(Point2I(x=10, y=10), self.cell_size))
        index = grid.index(Point2I(x=11, y=9))
        self.assertEqual(index, Index2D(x=3, y=3))
        self.assertIsInstance(index, Index2D)

    def test_index(self):
        """Test various inputs to UniformGrid.index."""
        grid = UniformGrid.from_bbox_cell_size(self.bbox, self.cell_size)
        self.assertEqual(grid.index(self.bbox.getMin()), Index2D(x=0, y=0))
        self.assertEqual(grid.index(self.bbox.getMax()), Index2D(x=4, y=5))
        self.assertEqual(grid.index(Point2I(x=9, y=5)), Index2D(x=2, y=1))
        self.assertEqual(grid.index(Point2I(x=9, y=6)), Index2D(x=2, y=2))
        self.assertEqual(grid.index(Point2I(x=10, y=5)), Index2D(x=3, y=1))
        self.assertEqual(grid.index(Point2I(x=10, y=6)), Index2D(x=3, y=2))
        with self.assertRaises(ValueError):
            grid.index(self.bbox.getMin() - Extent2I(x=0, y=1))
        with self.assertRaises(ValueError):
            grid.index(self.bbox.getMin() - Extent2I(x=1, y=0))
        with self.assertRaises(ValueError):
            grid.index(self.bbox.getMin() - Extent2I(x=1, y=1))
        with self.assertRaises(ValueError):
            grid.index(self.bbox.getMax() + Extent2I(x=0, y=1))
        with self.assertRaises(ValueError):
            grid.index(self.bbox.getMax() + Extent2I(x=1, y=0))
        with self.assertRaises(ValueError):
            grid.index(self.bbox.getMax() + Extent2I(x=1, y=1))

    def test_repr(self):
        """Test that UniformGrid.__repr__ round-trips through eval."""
        for grid in (
            UniformGrid.from_bbox_cell_size(self.bbox, self.cell_size),
            UniformGrid.from_bbox_shape(self.bbox, self.shape, padding=1),
            UniformGrid(self.cell_size, self.shape, min=self.bbox.getMin(), padding=3),
        ):
            self.assertEqual(eval(repr(grid)), grid, msg=repr(grid))

    @lsst.utils.tests.methodParameters(padding=[0, 3])
    def test_pickle(self, padding: int):
        """Test that UniformGrid objects are pickleable."""
        grid1 = UniformGrid.from_bbox_cell_size(self.bbox, self.cell_size, padding=padding)
        grid2 = pickle.loads(pickle.dumps(grid1, pickle.HIGHEST_PROTOCOL))
        self.assertIsInstance(grid2, UniformGrid)
        self.assertEqual(grid1, grid2)

    @lsst.utils.tests.methodParameters(padding=(1, 2, 3, 4, 7, 10))
    def test_padding(self, padding: int):
        """Test that bbox_of and index methods work with padding > 0."""
        grid = UniformGrid.from_bbox_cell_size(self.bbox, self.cell_size, padding=padding)

        # Test all interior cells
        for x in range(1, self.nx - 1):
            for y in range(1, self.ny - 1):
                bbox = grid.bbox_of(Index2D(x=x, y=y))
                self.assertEqual(bbox.getDimensions(), self.cell_size)
                self.assertEqual(bbox.getMin(), Point2I(x=x * self.cw + self.x0, y=y * self.ch + self.y0))

        # Test the four corners
        for x in (0, self.nx - 1):
            for y in (0, self.ny - 1):
                bbox = grid.bbox_of(Index2D(x=x, y=y))
                self.assertEqual(
                    bbox.getDimensions(),
                    self.cell_size + Extent2I(padding, padding),
                )
                self.assertEqual(
                    bbox.getMin(),
                    Point2I(
                        x=x * self.cw + self.x0 - padding * (x == 0),
                        y=y * self.ch + self.y0 - padding * (y == 0),
                    ),
                )

        # Test along the two horizontal edges
        for x in (0, self.nx - 1):
            for y in range(1, self.ny - 1):
                bbox = grid.bbox_of(Index2D(x=x, y=y))
                self.assertEqual(bbox.getDimensions(), self.cell_size + Extent2I(padding, 0))
                self.assertEqual(
                    bbox.getMin(),
                    Point2I(x=x * self.cw + self.x0 - padding * (x == 0), y=y * self.ch + self.y0),
                )

        # Test along the two vertical edges
        for x in range(1, self.nx - 1):
            for y in (0, self.ny - 1):
                bbox = grid.bbox_of(Index2D(x=x, y=y))
                self.assertEqual(bbox.getDimensions(), self.cell_size + Extent2I(0, padding))
                self.assertEqual(
                    bbox.getMin(),
                    Point2I(x=x * self.cw + self.x0, y=y * self.ch + self.y0 - padding * (y == 0)),
                )

        # Check the mapping between positions and indices
        positions_index = (
            # Check the four interior corners
            (Point2I(x=self.x0, y=self.y0), Index2D(x=0, y=0)),  # Lower-left corner, in the x and y buffer
            (
                Point2I(x=self.x0 + self.bw - 1, y=self.y0),
                Index2D(x=self.nx - 1, y=0),
            ),  # Lower-right corner, in the x buffer
            (
                Point2I(x=self.x0, y=self.y0 + self.bh - 1),
                Index2D(x=0, y=self.ny - 1),
            ),  # Upper-left corner, in the y buffer
            (
                Point2I(x=self.x0 + self.bw - 1, y=self.y0 + self.bh - 1),
                Index2D(x=self.nx - 1, y=self.ny - 1),
            ),  # Upper-right corner, in the x and y buffer
            # Check for one cell from the origin in either directions
            (Point2I(x=self.x0, y=self.y0 + self.ch), Index2D(x=0, y=1)),  # Left edge, in the x buffer
            (Point2I(x=self.x0 + self.cw, y=self.y0), Index2D(x=1, y=0)),
            # Check the four corners of the lower left buffer region
            (Point2I(x=self.x0 - padding, y=self.y0), Index2D(0, 0)),
            (Point2I(x=self.x0, y=self.y0 - padding), Index2D(0, 0)),
            (Point2I(x=self.x0 - padding, y=self.y0), Index2D(0, 0)),
            (Point2I(x=self.x0 - padding, y=self.y0 - padding), Index2D(0, 0)),
            # Check the four corners of the lower right buffer region
            (Point2I(x=self.x0 + self.bw, y=self.y0), Index2D(self.nx - 1, 0)),
            (Point2I(x=self.x0 + self.bw + padding - 1, y=self.y0), Index2D(self.nx - 1, 0)),
            (Point2I(x=self.x0 + self.bw + padding - 1, y=self.y0 - padding), Index2D(self.nx - 1, 0)),
            (Point2I(x=self.x0 + self.bw, y=self.y0 - padding), Index2D(self.nx - 1, 0)),
            # Check the four corners of the upper left buffer region
            (Point2I(x=self.x0, y=self.y0 + self.bh), Index2D(0, self.ny - 1)),
            (Point2I(x=self.x0, y=self.y0 + self.bh + padding - 1), Index2D(0, self.ny - 1)),
            (Point2I(x=self.x0 - padding, y=self.y0 + self.bh + padding - 1), Index2D(0, self.ny - 1)),
            (Point2I(x=self.x0 - padding, y=self.y0 + self.bh), Index2D(0, self.ny - 1)),
            # Check the four corners of the upper right buffer region
            (Point2I(x=self.x0 + self.bw, y=self.y0 + self.bh), Index2D(self.nx - 1, self.ny - 1)),
            (
                Point2I(x=self.x0 + self.bw, y=self.y0 + self.bh + padding - 1),
                Index2D(self.nx - 1, self.ny - 1),
            ),
            (
                Point2I(x=self.x0 + self.bw + padding - 1, y=self.y0 + self.bh),
                Index2D(self.nx - 1, self.ny - 1),
            ),
            (
                Point2I(x=self.x0 + self.bw + padding - 1, y=self.y0 + self.bh + padding - 1),
                Index2D(self.nx - 1, self.ny - 1),
            ),
        )

        for position, index in positions_index:
            self.assertEqual(grid.index(position), index, msg=f"{position} does not map to {index}")

        # Check that positions outside the grid raise an exception
        self.assertRaises(ValueError, grid.index, Point2I(x=self.x0 - padding - 1, y=self.y0))
        self.assertRaises(ValueError, grid.index, Point2I(x=self.x0 + self.bw + padding, y=self.y0))
        self.assertRaises(ValueError, grid.index, Point2I(x=self.x0, y=self.y0 - padding - 1))
        self.assertRaises(ValueError, grid.index, Point2I(x=self.x0, y=self.y0 + self.bh + padding))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    """Test for memory/resource leaks."""


def setup_module(module):  # noqa: D103
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
