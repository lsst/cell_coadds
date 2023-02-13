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

from lsst.cell_coadds import UniformGrid
from lsst.geom import Box2I, Extent2I, Point2I
from lsst.skymap import Index2D


class UniformGridTestCase(unittest.TestCase):
    """Tests for UniformGrid and GridIndex/Index2D's C++/Python
    translation."""

    def setUp(self) -> None:
        self.bbox = Box2I(Point2I(x=1, y=2), Extent2I(x=15, y=12))
        self.cell_size = Extent2I(x=3, y=2)
        self.shape = Index2D(x=5, y=6)

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
        grid = UniformGrid(self.cell_size, self.shape, self.bbox.getMin())
        self._check(grid)

    def _check(self, grid: UniformGrid) -> None:
        self.assertEqual(grid.bbox, self.bbox)
        self.assertEqual(grid.cell_size, self.cell_size)
        self.assertEqual(grid.shape, self.shape)
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
        grid = UniformGrid.from_bbox_cell_size(self.bbox, self.cell_size)
        self.assertEqual(eval(repr(grid)), grid, msg=repr(grid))
        # TODO: Add more test cases, especially with non-zero min if this
        # is a good test.

    def test_index_overloads(self):
        """Test methods that accept either a single (x, y) object argument or
        kw-only x and y args.
        """
        return
        grid = UniformGrid.from_bbox_cell_size(self.bbox, self.cell_size)
        self.assertEqual(grid.index(Point2I(x=9, y=5)), grid.index(x=9, y=5))
        self.assertEqual(grid.min_of(Index2D(x=1, y=3)), grid.min_of(x=1, y=3))
        self.assertEqual(grid.bbox_of(Index2D(x=1, y=3)), grid.bbox_of(x=1, y=3))
        with self.assertRaises(TypeError):
            grid.index(9, 5)
        with self.assertRaises(TypeError):
            grid.min_of(1, 3)
        with self.assertRaises(TypeError):
            grid.bbox_of(1, 3)
        # Unlike C++, Python does not allow overloading methods with different
        # signatures. Therefore, explicitly test that invalid calls to the
        # methods raise the error that we expect.
        with self.assertRaises(TypeError):
            grid.index(position=Point2I(x=9, y=5), x=9, y=5)
        with self.assertRaises(TypeError):
            grid.min_of(index=Index2D(x=1, y=3), x=1)
        with self.assertRaises(TypeError, msg="sdfs"):
            grid.bbox_of(x=1)
        with self.assertRaises(TypeError):
            grid.min_of(Index2D(x=1, y=3), x=1, y=3)
        with self.assertRaises(TypeError):
            grid.bbox_of(Index2D(x=1, y=3), 1, 3)

    def test_pickle(self):
        """Test that UniformGrid objects are pickleable."""
        grid1 = UniformGrid.from_bbox_cell_size(self.bbox, self.cell_size)
        grid2 = pickle.loads(pickle.dumps(grid1, pickle.HIGHEST_PROTOCOL))
        self.assertIsInstance(grid2, UniformGrid)
        self.assertEqual(grid1, grid2)


if __name__ == "__main__":
    unittest.main()
