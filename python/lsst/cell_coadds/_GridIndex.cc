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
#include "pybind11/pybind11.h"
#include "lsst/cell_coadds/python.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace pybind11 {
namespace detail {

bool type_caster<lsst::cell_coadds::GridIndex>::load(handle src, bool) {
    if (PyArg_ParseTuple(src.ptr(), "ii", &value.x, &value.y)) {
        return true;
    } else {
        PyErr_Clear();
        return false;
    }
};

handle type_caster<lsst::cell_coadds::GridIndex>::cast(
    lsst::cell_coadds::GridIndex src, return_value_policy /* policy */, handle /* parent */) {
    // Static variable to hold the named tuple's type object.
    static PyObject* py_type = nullptr;
    if (!py_type) {
        // Attempt to import that type object.
        PyObject* module = PyImport_ImportModule("lsst.skymap");
        if (!module) {
            throw error_already_set();
        }
        py_type = PyObject_GetAttrString(module, "Index2D");
        if (!py_type) {
            throw error_already_set();
        }
    }
    return PyObject_CallFunction(py_type, "ii", src.x, src.y);
}

}  // namespace detail
}  // namespace pybind11
