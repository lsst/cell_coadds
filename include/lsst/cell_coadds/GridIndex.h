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
#ifndef LSST_CELL_COADDS_GridIndex_h
#define LSST_CELL_COADDS_GridIndex_h

namespace lsst {
namespace cell_coadds {

/**
 * A 2-d index or shape in a grid.
 *
 * This class is mapped to the lsst.skymap.Index2D (x, y) namedtuple in Python
 * rather than being wrapped directly.
 */
struct GridIndex final {
    int x;
    int y;

    bool operator==(GridIndex const& other) const { return x == other.x && y == other.y; }
    bool operator!=(GridIndex const& other) const { return !(*this == other); }
};

}  // namespace cell_coadds
}  // namespace lsst

#endif  // !LSST_CELL_COADDS_GridIndex_h
