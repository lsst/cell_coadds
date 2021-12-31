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
#ifndef LSST_CELL_COADDS_fits_helpers_h
#define LSST_CELL_COADDS_fits_helpers_h

#include <memory>
#include <optional>
#include <vector>

#include "lsst/daf/base/PropertyList.h"
#include "lsst/afw/fits.h"
#include "lsst/afw/image/Image.h"
#include "lsst/afw/image/Mask.h"
#include "lsst/afw/table.h"
#include "lsst/afw/image/ImageFitsReader.h"
#include "lsst/afw/image/MaskFitsReader.h"
#include "lsst/afw/table/io/OutputArchive.h"
#include "lsst/afw/table/io/Persistable.h"
#include "lsst/afw/table/io/InputArchive.h"

namespace lsst {
namespace cell_coadds {

class FitsWriteHelper {
public:
    FitsWriteHelper(
        std::string const& filename,
        std::shared_ptr<daf::base::PropertyList const> primary_metadata = nullptr);

    FitsWriteHelper(FitsWriteHelper const&) = delete;
    FitsWriteHelper(FitsWriteHelper&&) = delete;

    FitsWriteHelper& operator=(FitsWriteHelper const&) = delete;
    FitsWriteHelper& operator=(FitsWriteHelper&&) = delete;

    void write_image(
        afw::image::Image<float> const& image,
        daf::base::PropertyList const* metadata = nullptr,
        afw::fits::ImageWriteOptions const* options = nullptr);

    void write_mask(
        afw::image::Mask<> const& mask,
        daf::base::PropertyList const* metadata = nullptr,
        afw::fits::ImageWriteOptions const* options = nullptr);

    void write_catalog(afw::table::BaseCatalog const& catalog);

    int put_persistable(std::shared_ptr<afw::table::io::Persistable const> obj, bool permissive = false);

    void finish();

private:
    afw::fits::Fits _fits;
    std::optional<afw::table::io::OutputArchive> _archive;
};

class FitsReadHelper {
public:
    explicit FitsReadHelper(std::string const& filename);

    FitsReadHelper(FitsReadHelper const&) = delete;
    FitsReadHelper(FitsReadHelper&&) = delete;

    FitsReadHelper& operator=(FitsReadHelper const&) = delete;
    FitsReadHelper& operator=(FitsReadHelper&&) = delete;

    std::shared_ptr<daf::base::PropertyList> get_primary_metadata() const { return _primary_metadata; }

    std::string get_filename() const { return _fits.getFileName(); }

    std::unique_ptr<afw::image::ImageFitsReader> read_image(int hdu);

    std::unique_ptr<afw::image::MaskFitsReader> read_mask(int hdu);

    afw::table::BaseCatalog read_catalog(int hdu);

    afw::table::io::InputArchive read_archive(int hdu);

private:
    afw::fits::Fits _fits;
    std::shared_ptr<daf::base::PropertyList> _primary_metadata;
};

}  // namespace cell_coadds
}  // namespace lsst

#endif  // !LSST_CELL_COADDS_fits_helpers_h
