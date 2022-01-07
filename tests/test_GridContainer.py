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

import copy
import unittest
from typing import Dict, List

from lsst.cell_coadds import GridContainer, GridContainerBuilder
from lsst.skymap import Index2D


class GridContainerTestCase(unittest.TestCase):
    """Tests for GridContainer and GridIndex/Index2D's C++/Python
    translation."""

    def _fill(self, builder: GridContainerBuilder[Dict[str, int]]) -> None:
        """Populate a GridContainerBuilder with dicts that map "x" or "y" to
        the cell index in that dimension.
        """
        for y in range(builder.offset.y, builder.offset.y + builder.shape.y):
            for x in range(builder.offset.x, builder.offset.x + builder.shape.x):
                builder.set(x=x, y=y, value={"x": x, "y": y})

    def _check(self, container: GridContainer[Dict[str, int]]) -> None:
        """Perform a complete battery of tests on a GridContainer instance."""
        for value in container:
            self.assertEqual(container.get(x=value["x"], y=value["y"]), value)
            self.assertEqual(container.get(Index2D(x=value["x"], y=value["y"])), value)
            self.assertEqual(container[Index2D(**value)], value)
        self.assertEqual(container.first, {"x": container.offset.x, "y": container.offset.y})
        self.assertEqual(
            container.last,
            {
                "x": container.offset.x + container.shape.x - 1,
                "y": container.offset.y + container.shape.y - 1,
            },
        )
        transformed_builder = container.rebuild_transformed(lambda cell: list(cell.values()))
        transformed: GridContainer[List[int]] = transformed_builder.finish()
        self.assertEqual(transformed.shape, container.shape)
        self.assertEqual(transformed.offset, container.offset)
        self.assertEqual(list(transformed), [list(cell.values()) for cell in container])
        emptied_builder: GridContainerBuilder[Dict[str, int]] = transformed.rebuild_empty()
        self._fill(emptied_builder)
        emptied: GridContainer[Dict[str, int]] = emptied_builder.finish()
        self.assertEqual(emptied.shape, container.shape)
        self.assertEqual(emptied.offset, container.offset)
        self.assertEqual(list(emptied), list(container))
        copied = copy.copy(container)
        self.assertEqual(copied.shape, container.shape)
        self.assertEqual(copied.offset, container.offset)
        self.assertEqual(list(copied), list(container))
        self.assertTrue(
            all(copied_cell is original_cell for copied_cell, original_cell in zip(copied, container))
        )
        deep_copied = copy.deepcopy(container)
        self.assertEqual(deep_copied.shape, container.shape)
        self.assertEqual(deep_copied.offset, container.offset)
        self.assertEqual(list(deep_copied), list(container))
        self.assertTrue(
            all(
                deep_copied_cell is not original_cell
                for deep_copied_cell, original_cell in zip(deep_copied, container)
            )
        )

    def test_simple_ctor(self):
        """Test a GridContainer built with the shape-only builder
        constructor."""
        shape = Index2D(x=3, y=2)
        builder: GridContainerBuilder[Dict[str, int]] = GridContainerBuilder(shape)
        self.assertEqual(builder.shape, shape)
        self.assertIsInstance(builder.shape, Index2D)
        self.assertEqual(builder.offset, Index2D(x=0, y=0))
        self.assertIsInstance(builder.offset, Index2D)
        self.assertEqual(len(builder), shape[0] * shape[1])
        with self.assertRaises(Exception):
            builder.finish()
        self._fill(builder)
        self._check(builder.finish())

    def test_complex_ctor(self):
        """Test a GridContainer built with the shape-and-offset builder
        constructor."""
        shape = Index2D(x=3, y=2)
        offset = Index2D(x=1, y=2)
        builder: GridContainerBuilder[Dict[str, int]] = GridContainerBuilder(shape, offset)
        self.assertEqual(builder.shape, shape)
        self.assertIsInstance(builder.shape, Index2D)
        self.assertEqual(builder.offset, offset)
        self.assertIsInstance(builder.offset, Index2D)
        self.assertEqual(len(builder), shape[0] * shape[1])
        with self.assertRaises(Exception):
            builder.finish()
        self._fill(builder)
        self._check(builder.finish())


if __name__ == "__main__":
    unittest.main()
