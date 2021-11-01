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

UniformGrid::UniformGrid(geom::Box2I const& bbox, geom::Extent2I const& stride)
        : _bbox(bbox),
          _stride(stride),
          _shape{bbox.getWidth() / stride.getX(), bbox.getHeight() / stride.getY()} {
    if (_bbox.getWidth() % _stride.getX() != 0) {
        throw LSST_EXCEPT(pex::exceptions::LengthError,
                          (boost::format("Bounding box width %s is not evenly divided by x stride %s.") %
                           _bbox.getWidth() % _stride.getX())
                                  .str());
    }
    if (_bbox.getHeight() % _stride.getY() != 0) {
        throw LSST_EXCEPT(pex::exceptions::LengthError,
                          (boost::format("Bounding box height %s is not evenly divided by y stride %s.") %
                           _bbox.getHeight() % _stride.getY())
                                  .str());
    }
}

UniformGrid::UniformGrid(geom::Box2I const& bbox, Index const& shape)
        : _bbox(bbox), _stride(bbox.getWidth() / shape.x, bbox.getHeight() / shape.y), _shape(shape) {
    if (_bbox.getWidth() % _shape.x != 0) {
        throw LSST_EXCEPT(pex::exceptions::LengthError,
                          (boost::format("Bounding box width %s is not evenly divided by x shape %s.") %
                           _bbox.getWidth() % _shape.x)
                                  .str());
    }
    if (_bbox.getHeight() % _shape.y != 0) {
        throw LSST_EXCEPT(pex::exceptions::LengthError,
                          (boost::format("Bounding box height %s is not evenly divided by y shape %s.") %
                           _bbox.getHeight() % _shape.y)
                                  .str());
    }
}

UniformGrid::Index UniformGrid::index(geom::Point2I const& position) const {
    geom::Extent2I offset = position - _bbox.getBegin();
    Index result = {offset.getX() / _stride.getX(), offset.getY() / _stride.getY()};
    if (result.x < 0 || result.x >= _shape.x) {
        throw LSST_EXCEPT(
                pex::exceptions::LengthError,
                (boost::format("Position %s is not within bounding box %s.") % position % _bbox).str());
    }
    return result;
}

UniformGrid UniformGrid::subset(Index const& min, Index const& max) const {
    return UniformGrid(geom::Box2I(geom::Point2I(min.x * _stride.getX() + _bbox.getBeginX(),
                                                 min.y * _stride.getY() * _bbox.getBeginY()),
                                   geom::Extent2I((1 + max.x - min.x) * _stride.getX(),
                                                  (1 + max.y - min.y) * _stride.getY())),
                       _stride);
}

geom::Box2I UniformGrid::bbox_of(Index const& index) const {
    return geom::Box2I(geom::Point2I(index.x * _stride.getX() + _bbox.getBeginX(),
                                     index.y * _stride.getY() + _bbox.getBeginY()),
                       _stride);
}

}  // namespace cell_coadds
}  // namespace lsst
