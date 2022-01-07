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
#include "lsst/cpputils/python.h"
#include "lsst/cell_coadds/UniformGrid.h"
#include "lsst/cell_coadds/python.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace cell_coadds {

void wrapUniformGrid(utils::python::WrapperCollection& wrappers) {
    wrappers.wrapType(py::class_<UniformGrid>(wrappers.module, "UniformGrid"), [](auto& mod, auto& cls) {
        cls.def(py::init<geom::Box2I const&, geom::Extent2I const&>(), "bbox"_a, "cell_size"_a);
        cls.def(py::init<geom::Box2I const&, GridIndex const&>(), "bbox"_a, "shape"_a);
        cls.def(
            py::init<geom::Extent2I const&, GridIndex const&, geom::Point2I const&>(),
            "cell_size"_a,
            "shape"_a,
            "min"_a = geom::Point2I());
        cls.def("__eq__", &UniformGrid::operator==, py::is_operator());
        cls.def("index", &UniformGrid::index, "position"_a);
        cls.def("min_of", &UniformGrid::min_of, "index"_a);
        cls.def("bbox_of", &UniformGrid::bbox_of, "index"_a);
        cls.def_property_readonly("bbox", &UniformGrid::get_bbox, py::return_value_policy::copy);
        cls.def_property_readonly("cell_size", &UniformGrid::get_cell_size, py::return_value_policy::copy);
        cls.def_property_readonly("shape", &UniformGrid::get_shape, py::return_value_policy::copy);
    });
}

}  // namespace cell_coadds
}  // namespace lsst
