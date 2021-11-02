// -*- LSST-C++ -*-
/*
 * This file is part of cell_coadds.
 *
 * Developed for the LSST Data Management System.
 * This product includes software developed by the LSST Project
 * (https://www.lsst.org).
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#ifndef LSST_CELL_COADDS_UniformGrid_h
#define LSST_CELL_COADDS_UniformGrid_h

#include <vector>
#include "lsst/geom/Box.h"

namespace lsst {
namespace cell_coadds {

/**
 * A 2-d integer grid.
 */
class UniformGrid {
public:
    /**
     * A 2-d index or shape in a grid.
     *
     * This class is mapped to a (y, x) tuple in Python rather than being
     * wrapped directly.
     */
    struct Index {
        int x;
        int y;
    };

    /**
     * Construct from bounding box and stride.
     *
     * @param bbox    Bounding box of the full grid.
     *
     * @param stride  Size of each grid cell.  Must divide the bbox width
     *                and height evenly.
     */
    UniformGrid(geom::Box2I const& bbox, geom::Extent2I const& stride);

    /**
     * Construct from bounding box and stride.
     *
     * @param bbox    Bounding box of the full grid.
     *
     * @param shape   Number of cells in the grid in each dimension.  Must
     *                divide the bbox width and height evenly.
     */
    UniformGrid(geom::Box2I const& bbox, Index const& shape);

    /**
     * Find the index of the cell that contains the given point.
     */
    Index index(geom::Point2I const& position) const;

    /**
     * Flatten a 2-d index into an integer suitable for addressing a 1-d array.
     *
     * This is guaranteed to be row-major order (all cells in a row, and then
     * cells in the next row); this method is a convenience, not an attempt
     * to fully encapsulate the ordering.
     */
    std::size_t flatten(Index const& index) const { return _shape.x * index.y + index.x; }

    /**
     * Return a new grid with the same spacing containing only the given
     * index range.
     *
     * @param min   Minimum index to include, inclusive.
     * @param max   Maximum index to include, inclusive.
     *
     * This method is mapped to `__getitem__` in Python, with
     * `(y_slice, x_slice)` arguments.
     */
    UniformGrid subset(Index const& min, Index const& max) const;

    /**
     * Return the bounding box of a single cell.
     */
    geom::Box2I bbox_of(Index const& index) const;

    /**
     * Return the bounding box of the full grid.
     */
    geom::Box2I get_bbox() const { return _bbox; }

    /**
     * Return the dimensions of each cell (also the separation between them).
     */
    geom::Extent2I get_cell_size() const { return _cell_size; }

    /**
     * Return the number of cells in each dimension.
     */
    Index get_shape() const { return _shape; }

private:
    geom::Box2I _bbox;
    geom::Extent2I _cell_size;
    Index _shape;
};

}  // namespace cell_coadds
}  // namespace lsst

#endif  // !LSST_CELL_COADDS_UniformGrid_h
