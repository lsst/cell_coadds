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
import pickle
import unittest

from lsst.cell_coadds import GridContainer, UniformGrid
from lsst.geom import Box2I, Extent2I, Point2I
from lsst.skymap import Index2D
from lsst.utils.tests import methodParameters


class GridContainerTestCase(unittest.TestCase):
    """Tests for GridContainer and GridIndex/Index2D's C++/Python
    translation.
    """

    def _fill(self, container: GridContainer[dict[str, int]]) -> None:
        """Populate a GridContainer with dicts that map "x" or "y" to
        the cell index in that dimension.
        """
        for y in range(container.offset.y, container.offset.y + container.shape.y):
            for x in range(container.offset.x, container.offset.x + container.shape.x):
                container[Index2D(x, y)] = {"x": x, "y": y}

    def _check(self, container: GridContainer[dict[str, int]]) -> None:
        """Perform a complete battery of tests on a GridContainer instance."""
        for value in container.values():
            self.assertEqual(container[Index2D(x=value["x"], y=value["y"])], value)
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
        transformed: GridContainer[list[int]] = container.rebuild_transformed(
            lambda cell: list(cell.values())  # type: ignore
        )
        self.assertEqual(transformed.shape, container.shape)
        self.assertEqual(transformed.offset, container.offset)
        self.assertEqual(list(transformed.keys()), list(container.keys()))
        self.assertEqual(list(transformed.values()), [list(v.values()) for v in container.values()])
        emptied = GridContainer[dict[str, int]](container.shape, container.offset)
        self._fill(emptied)
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
                for deep_copied_cell, original_cell in zip(deep_copied, container, strict=True)
            )
        )

    def test_simple_ctor(self) -> None:
        """Test a GridContainer built with the shape-only GridContainer
        constructor.
        """
        shape = Index2D(x=3, y=2)
        gc = GridContainer(shape)
        self.assertEqual(gc.shape, shape)
        self.assertIsInstance(gc.shape, Index2D)
        self.assertEqual(gc.offset, Index2D(x=0, y=0))
        self.assertIsInstance(gc.offset, Index2D)
        self.assertEqual(gc.size, shape[0] * shape[1])
        self.assertEqual(len(gc), 0)  # unfilled container has length 0.
        self._fill(gc)
        self.assertEqual(len(gc), shape[0] * shape[1])
        self._check(gc)

    def test_complex_ctor(self) -> None:
        """Test a GridContainer built with the shape-and-offset GridContainer
        constructor.
        """
        shape = Index2D(x=3, y=2)
        offset = Index2D(x=1, y=2)
        gc = GridContainer(shape, offset)
        self.assertEqual(gc.shape, shape)
        self.assertIsInstance(gc.shape, Index2D)
        self.assertEqual(gc.offset, offset)
        self.assertIsInstance(gc.offset, Index2D)
        self.assertEqual(len(gc), 0)
        self.assertEqual(gc.size, shape[0] * shape[1])
        self._fill(gc)
        self.assertEqual(len(gc), shape[0] * shape[1])
        self._check(gc)

    def test_subset_overlapping(self) -> None:
        """Test various inputs to GridContainer.subset_overlapping."""
        cell_size = Extent2I(x=3, y=2)
        full_bbox = Box2I(Point2I(x=1, y=2), Extent2I(x=15, y=12))
        full_shape = Index2D(x=5, y=6)
        full_container = GridContainer(shape=full_shape)
        self._fill(full_container)
        grid = UniformGrid.from_bbox_cell_size(full_bbox, cell_size)
        self.assertEqual(grid.shape, full_shape)

        # Subset with the orignal bounding box; should behave like a deepcopy.
        subset_container_full = full_container.subset_overlapping(grid, full_bbox)
        self.assertEqual(subset_container_full.shape, full_container.shape)
        self.assertEqual(subset_container_full.offset, full_container.offset)
        self.assertEqual(list(subset_container_full), list(full_container))
        subset_container_full[Index2D(x=2, y=2)] = {"x": -1, "y": -1}
        self.assertEqual(list(subset_container_full.keys()), list(full_container.keys()))
        self.assertNotEqual(list(subset_container_full.values()), list(full_container.values()))

        # Subset the full container with a nontrivial bbox.
        bbox_1 = Box2I(Point2I(x=6, y=4), Point2I(x=10, y=7))
        subset_container_1 = full_container.subset_overlapping(grid, bbox_1)
        self.assertEqual(subset_container_1.offset, Index2D(x=1, y=1))
        self.assertEqual(subset_container_1.shape, Index2D(x=3, y=2))
        union_bbox_1 = Box2I()
        for v in subset_container_1.keys():  # noqa: SIM118
            cell_bbox = grid.bbox_of(v)
            self.assertTrue(cell_bbox.overlaps(bbox_1))
            union_bbox_1.include(cell_bbox)
        self.assertTrue(union_bbox_1.contains(bbox_1))
        self._check(subset_container_1)

        # Subset the subset container with an even smaller bbox, to check the
        # case where the original offset is nonzero.
        bbox_2 = Box2I(Point2I(x=6, y=5), Point2I(x=7, y=5))
        subset_container_2 = subset_container_1.subset_overlapping(grid, bbox_2)
        self.assertEqual(subset_container_2.offset, Index2D(x=1, y=1))
        self.assertEqual(subset_container_2.shape, Index2D(x=2, y=1))
        union_bbox_2 = Box2I()
        for v in subset_container_2.keys():  # noqa: SIM118
            cell_bbox = grid.bbox_of(v)
            self.assertTrue(cell_bbox.overlaps(bbox_1))
            union_bbox_2.include(cell_bbox)
        self.assertTrue(union_bbox_2.contains(bbox_2))
        self._check(subset_container_2)

        # Subsetting the container by a bbox that isn't contained by the cells
        # raise a KeyError.
        with self.assertRaises(KeyError):
            subset_container_1.subset_overlapping(grid, full_bbox)

    def test_rebuild_transformed(self) -> None:
        """Test that rebuild_transformed method works by squaring
        transformation.
        """
        container = GridContainer(shape=Index2D(x=3, y=2))
        self._fill(container)
        container = container.rebuild_transformed(
            transform=lambda cell: {key: value**2 for key, value in cell.items()}
        )
        for index in container.indices():
            self.assertTrue(container[index]["x"] == index.x**2)
            self.assertTrue(container[index]["y"] == index.y**2)

    @methodParameters(offset=(Index2D(x=1, y=2), None))
    def test_pickle(self, offset) -> None:
        """Test that we can serialize GridContainer with pickle."""
        shape = Index2D(x=3, y=2)
        gc = GridContainer(shape, offset)
        self._fill(gc)
        grid_container = gc
        pickled_container = pickle.loads(pickle.dumps(grid_container, pickle.HIGHEST_PROTOCOL))
        self.assertIsInstance(pickled_container, GridContainer)
        # This line below is failing, so compare internals.
        # self.assertEqual(pickled_container, grid_container)
        self.assertEqual(pickled_container.shape, grid_container.shape)
        self.assertEqual(pickled_container.offset, grid_container.offset)
        self.assertEqual(pickled_container.first, grid_container.first)
        self.assertEqual(pickled_container.last, grid_container.last)
        self.assertListEqual(list(pickled_container.__iter__()), list(grid_container.__iter__()))
        self._check(pickled_container)


if __name__ == "__main__":
    unittest.main()
