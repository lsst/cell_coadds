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

#include "lsst/cell_coadds/UniformGrid.h"

namespace lsst {
namespace cell_coadds {

UniformGrid::UniformGrid(geom::Box2I const& bbox, geom::Extent2I const& cell_size)
        : _bbox(bbox),
          _cell_size(cell_size),
          _shape{bbox.getWidth() / cell_size.getX(), bbox.getHeight() / cell_size.getY()} {
    if (_bbox.getWidth() % _cell_size.getX() != 0) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthError,
            (boost::format("Bounding box width %s is not evenly divided by x cell_size %s.") %
             _bbox.getWidth() % _cell_size.getX())
                .str());
    }
    if (_bbox.getHeight() % _cell_size.getY() != 0) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthError,
            (boost::format("Bounding box height %s is not evenly divided by y cell_size %s.") %
             _bbox.getHeight() % _cell_size.getY())
                .str());
    }
}

UniformGrid::UniformGrid(geom::Box2I const& bbox, Index const& shape)
        : _bbox(bbox), _cell_size(bbox.getWidth() / shape.x, bbox.getHeight() / shape.y), _shape(shape) {
    if (_bbox.getWidth() % _shape.x != 0) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthError,
            (boost::format("Bounding box width %s is not evenly divided by x shape %s.") % _bbox.getWidth() %
             _shape.x)
                .str());
    }
    if (_bbox.getHeight() % _shape.y != 0) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthError,
            (boost::format("Bounding box height %s is not evenly divided by y shape %s.") %
             _bbox.getHeight() % _shape.y)
                .str());
    }
}

UniformGrid::UniformGrid(geom::Extent2I const& cell_size, Index const& shape, geom::Point2I const& min)
        : _bbox(min, geom::Extent2I(cell_size.getX() * shape.x, cell_size.getY() * shape.y)),
          _cell_size(cell_size),
          _shape(shape) {}

UniformGrid::Index UniformGrid::index(geom::Point2I const& position) const {
    geom::Extent2I offset = position - _bbox.getBegin();
    Index result = {offset.getX() / _cell_size.getX(), offset.getY() / _cell_size.getY()};
    if (result.x < 0 || result.x >= _shape.x) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthError,
            (boost::format("Position %s is not within bounding box %s.") % position % _bbox).str());
    }
    return result;
}

geom::Point2I UniformGrid::min_of(Index const& index) const {
    return geom::Point2I(
        index.x * _cell_size.getX() + _bbox.getBeginX(), index.y * _cell_size.getY() + _bbox.getBeginY());
}

}  // namespace cell_coadds
}  // namespace lsst
