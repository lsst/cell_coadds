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

import unittest

from lsst.cell_coadds import UniformGrid
from lsst.geom import Box2I, Extent2I, Point2I
from lsst.skymap import Index2D


class UniformGridTestCase(unittest.TestCase):
    """Tests for UniformGrid and GridIndex/Index2D's C++/Python
    translation."""

    def setUp(self) -> None:
        self.bbox = Box2I(Point2I(1, 2), Extent2I(15, 12))
        self.cell_size = Extent2I(3, 2)
        self.shape = (5, 6)

    def test_ctor_bbox_cell_size(self) -> None:
        """Test UniformGrid after construction with (bbox, cell_size)."""
        self._check(UniformGrid(self.bbox, self.cell_size))

    def test_ctor_bbox_shape(self) -> None:
        """Test UniformGrid after construction with (bbox, shape)."""
        self._check(UniformGrid(self.bbox, self.shape))

    def test_ctor_cell_size_shape_min(self) -> None:
        """Test UniformGrid after construction with (cell_size, shape, min)."""
        self._check(UniformGrid(self.cell_size, self.shape, self.bbox.getMin()))

    def _check(self, grid: UniformGrid) -> None:
        self.assertEqual(grid.bbox, self.bbox)
        self.assertEqual(grid.cell_size, self.cell_size)
        self.assertEqual(grid.shape, self.shape)
        self.assertIsInstance(grid.shape, Index2D)
        self.assertEqual(grid.bbox_of((3, 4)), Box2I(Point2I(10, 10), self.cell_size))
        index = grid.index(Point2I(11, 9))
        self.assertEqual(index, (3, 3))
        self.assertIsInstance(index, Index2D)


if __name__ == "__main__":
    unittest.main()
