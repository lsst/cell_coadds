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
#ifndef LSST_CELL_COADDS_python_h
#define LSST_CELL_COADDS_python_h

#include "lsst/cell_coadds/GridIndex.h"
#include "pybind11/pybind11.h"

namespace pybind11 {
namespace detail {

/*
 *  Custom type-caster that translates the C++ GridIndex struct into the Python
 *  lsst.skymap.Index2D named tuple, and converts arbitrary Python 2-tuples
 *  into C++ GridIndex.  If we add a named tuple for cells to the skymap
 *  package, we can switch to using that here instead.
 */
template <>
struct type_caster<lsst::cell_coadds::GridIndex> {
    PYBIND11_TYPE_CASTER(lsst::cell_coadds::GridIndex, _("Index2D"));

public:
    bool load(handle src, bool);
    static handle
    cast(lsst::cell_coadds::GridIndex src, return_value_policy /* policy */, handle /* parent */);
};

}  // namespace detail
}  // namespace pybind11

#endif  // !LSST_CELL_COADDS_python_h
