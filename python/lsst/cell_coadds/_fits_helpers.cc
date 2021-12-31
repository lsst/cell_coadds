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
#include "lsst/cell_coadds/fits_helpers.h"
#include "lsst/cell_coadds/python.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace cell_coadds {

void wrap_fits_helpers(utils::python::WrapperCollection& wrappers) {
    wrappers.addSignatureDependency("lsst.afw.image");
    wrappers.addSignatureDependency("lsst.afw.table.io");
    wrappers.wrapType(
        py::class_<FitsWriteHelper, std::shared_ptr<FitsWriteHelper>>(wrappers.module, "FitsWriteHelper"),
        [](auto& mod, auto& cls) {
            cls.def(
                py::init<std::string const&, std::shared_ptr<daf::base::PropertyList>>(),
                "filename"_a,
                "primary_metadata"_a = nullptr);
            cls.def(
                "write_image",
                &FitsWriteHelper::write_image,
                "image"_a,
                "metadata"_a = nullptr,
                "options"_a = nullptr);
            cls.def(
                "write_mask",
                &FitsWriteHelper::write_mask,
                "mask"_a,
                "metadata"_a = nullptr,
                "options"_a = nullptr);
            cls.def("write_catalog", &FitsWriteHelper::write_catalog, "catalog"_a);
            cls.def("put_persistable", &FitsWriteHelper::put_persistable, "obj"_a, "permissive"_a = false);
            cls.def("finish", &FitsWriteHelper::finish);
        });
    wrappers.wrapType(
        py::class_<FitsReadHelper, std::shared_ptr<FitsReadHelper>>(wrappers.module, "FitsReadHelper"),
        [](auto& mod, auto& cls) {
            cls.def(py::init<std::string const&>(), "filename"_a);
            cls.def_property_readonly("primary_metadata", &FitsReadHelper::get_primary_metadata);
            cls.def_property_readonly("filename", &FitsReadHelper::get_filename);
            // Wrap lambdas that convert result to shared_ptr for read_image
            // and read_mask; pybind11 isn't smart enough to do this
            // automatically, and the wrapper declaration says that they're
            // held by shared_ptr.
            // We also use keep_alive<0, 1>() to tell pybind11 to keep the
            // FitsReadHelper instance alive as long as any returned
            // ImageFitsReader or MaskFitsReader is alive; this is necessary
            // because the returned readers utilize the parent helper's
            // underlying afw::fits::Fits pointer.
            cls.def(
                "read_image",
                [](FitsReadHelper& self, int hdu) {
                    return std::shared_ptr<afw::image::ImageFitsReader>(self.read_image(hdu));
                },
                "hdu"_a,
                py::keep_alive<0, 1>());
            cls.def(
                "read_mask",
                [](FitsReadHelper& self, int hdu) {
                    return std::shared_ptr<afw::image::MaskFitsReader>(self.read_mask(hdu));
                },
                "hdu"_a,
                py::keep_alive<0, 1>());
            cls.def("read_catalog", &FitsReadHelper::read_catalog, "hdu"_a);
            cls.def("read_archive", &FitsReadHelper::read_archive, "hdu"_a);
        });
}

}  // namespace cell_coadds
}  // namespace lsst
