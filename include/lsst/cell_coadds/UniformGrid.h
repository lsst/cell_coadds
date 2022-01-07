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

#include "lsst/cell_coadds/GridIndex.h"
#include "lsst/geom/Box.h"

namespace lsst {
namespace cell_coadds {

/**
 * A 2-d integer grid.
 */
class UniformGrid final {
public:
    using Index = GridIndex;

    /**
     * Construct from bounding box and cell size.
     *
     * @param bbox       Bounding box of the full grid.
     *
     * @param cell_size  Size of each grid cell.  Must divide the bbox width
     *                   and height evenly.
     */
    UniformGrid(geom::Box2I const& bbox, geom::Extent2I const& cell_size);

    /**
     * Construct from bounding box and shape.
     *
     * @param bbox    Bounding box of the full grid.
     *
     * @param shape   Number of cells in the grid in each dimension.  Must
     *                divide the bbox width and height evenly.
     */
    UniformGrid(geom::Box2I const& bbox, Index const& shape);

    /**
     * Construct from cell size and shape.
     *
     * @param cell_size  Size of each grid cell.
     *
     * @param shape      Number of cells in the grid in each dimension.
     *
     * @param min     Minimum x and y coordinates of the bounding box.
     */
    UniformGrid(
        geom::Extent2I const& cell_size, Index const& shape, geom::Point2I const& min = geom::Point2I());

    //@{
    /**
     * Equality comparison.
     */
    bool operator==(UniformGrid const& other) const;
    bool operator!=(UniformGrid const& other) const { return !(*this == other); }
    //@}

    /**
     * Find the index of the cell that contains the given point.
     */
    Index index(geom::Point2I const& position) const;

    /**
     * Return the minimum point of a single cell's bounding box.
     */
    geom::Point2I min_of(Index const& index) const;

    /**
     * Return the bounding box of a single cell.
     */
    geom::Box2I bbox_of(Index const& index) const { return geom::Box2I(min_of(index), _cell_size); }

    /**
     * Return the bounding box of the full grid.
     */
    geom::Box2I const& get_bbox() const { return _bbox; }

    /**
     * Return the dimensions of each cell (also the separation between them).
     */
    geom::Extent2I const& get_cell_size() const { return _cell_size; }

    /**
     * Return the number of cells in each dimension.
     */
    Index const& get_shape() const { return _shape; }

private:
    geom::Box2I _bbox;
    geom::Extent2I _cell_size;
    Index _shape;
};

}  // namespace cell_coadds
}  // namespace lsst

#endif  // !LSST_CELL_COADDS_UniformGrid_h
