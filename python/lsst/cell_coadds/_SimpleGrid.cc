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
#include "lsst/cell_coadds/SimpleGrid.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace cell_coadds {

namespace {

// Instead of wrapping SimpleGrid::Index directly, we map it to a (y, x) tuple,
// since that's what's natural for 2-d indexing in Python.  We go through
// std::pair since pybind11 will conver that to/from Python tuple
// automatically.

inline SimpleGrid::Index pair_to_index(std::pair<int, int> const& pair) {
    return SimpleGrid::Index{pair.second, pair.first};
}

inline std::pair<int, int> index_to_pair(SimpleGrid::Index const& index) {
    return std::make_pair(index.y, index.x);
}

}  // namespace

void wrapSimpleGrid(utils::python::WrapperCollection& wrappers) {
    wrappers.wrapType(py::class_<SimpleGrid>(wrappers.module, "SimpleGrid"), [](auto& mod, auto& cls) {
        cls.def(py::init<geom::Box2I const&, geom::Extent2I const&>(), "bbox"_a, "strides"_a);
        cls.def(py::init([](geom::Box2I const& bbox, std::pair<int, int> const& shape) {
                    return SimpleGrid(bbox, pair_to_index(shape));
                }),
                "bbox"_a, "shape"_a);
        cls.def(
                "index",
                [](SimpleGrid const& self, geom::Point2I const& position) {
                    return index_to_pair(self.index(position));
                },
                "position"_a);
        cls.def(
                "flatten",
                [](SimpleGrid const& self, std::pair<int, int> const& index) {
                    return self.flatten(pair_to_index(index));
                },
                "index"_a);
        cls.def("subset",
                [](SimpleGrid const& self, std::pair<int, int> const& max, std::pair<int, int> const& max) {
                    return self.subset(pair_to_index(min), pair_to_index(max));
                });
        cls.def(
                "bbox_of",
                [](SimpleGrid const& self, std::pair<int, int> const& index) {
                    return self.bbox_of(pair_to_index(index));
                },
                "position"_a);
        cls.def_property_readonly("bbox", &SimpleGrid::get_bbox);
        cls.def_property_readonly("stride", &SimpleGrid::get_stride);
        cls.def_property_readonly("shape",
                                  [](SimpleGrid const& self) { return index_to_pair(self.get_shape()); });
    });
}

}  // namespace cell_coadds
}  // namespace lsst
