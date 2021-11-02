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
#include "pybind11/stl.h"
#include "lsst/cpputils/python.h"
#include "lsst/cell_coadds/StitchedPsf.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace cell_coadds {

void wrapStitchedPsf(utils::python::WrapperCollection& wrappers) {
    wrappers.wrapType(
        py::class_<StitchedPsf, std::shared_ptr<StitchedPsf>, meas::algorithms::ImagePsf>(
            wrappers.module, "StitchedPsf"),
        [](auto& mod, auto& cls) {
            cls.def(
                py::init([](GridContainer<py::object> container, UniformGrid const& grid) {
                    auto builder = std::move(container).rebuild_transformed([](py::object obj) {
                        return py::cast<std::shared_ptr<afw::detection::Psf::Image>>(obj);
                    });
                    return StitchedPsf(std::move(builder).finish(), grid);
                }),
                "images"_a,
                "grid"_a);
            // All other methods are implementations of base
            // class methods, and can use their wrappers.
        });
}

}  // namespace cell_coadds
}  // namespace lsst
