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

#include "lsst/cell_coadds/fits_helpers.h"

namespace lsst {
namespace cell_coadds {

FitsWriteHelper::FitsWriteHelper(
    std::string const& filename, std::shared_ptr<daf::base::PropertyList const> primary_metadata)
        : _fits(filename, "w", afw::fits::Fits::AUTO_CLOSE | afw::fits::Fits::AUTO_CHECK), _archive() {
    _fits.createEmpty();
    _fits.writeMetadata(*primary_metadata);
}

void FitsWriteHelper::write_image(
    afw::image::Image<float> const& image,
    daf::base::PropertyList const* metadata,
    afw::fits::ImageWriteOptions const* options) {
    if (!options) {
        image.writeFits(_fits, metadata);
    } else {
        image.writeFits(_fits, *options, metadata);
    }
}

void FitsWriteHelper::write_mask(
    afw::image::Mask<> const& mask,
    daf::base::PropertyList const* metadata,
    afw::fits::ImageWriteOptions const* options) {
    if (!options) {
        mask.writeFits(_fits, metadata);
    } else {
        mask.writeFits(_fits, *options, metadata);
    }
}

void FitsWriteHelper::write_catalog(afw::table::BaseCatalog const& catalog) { catalog.writeFits(_fits); }

int FitsWriteHelper::put_persistable(
    std::shared_ptr<afw::table::io::Persistable const> obj, bool permissive) {
    if (!_archive) {
        _archive.emplace();
    }
    return _archive->put(obj, permissive);
}

void FitsWriteHelper::finish() {
    if (_archive) {
        _archive->writeFits(_fits);
    }
    _fits.closeFile();
}

FitsReadHelper::FitsReadHelper(std::string const& filename)
        : _fits(filename, "r", afw::fits::Fits::AUTO_CLOSE | afw::fits::Fits::AUTO_CHECK),
          _primary_metadata(std::make_shared<daf::base::PropertyList>()) {
    _fits.readMetadata(*_primary_metadata, true);
}

std::unique_ptr<afw::image::ImageFitsReader> FitsReadHelper::read_image(int hdu) {
    _fits.setHdu(hdu, false);
    return std::make_unique<afw::image::ImageFitsReader>(&_fits);
}

std::unique_ptr<afw::image::MaskFitsReader> FitsReadHelper::read_mask(int hdu) {
    _fits.setHdu(hdu, false);
    return std::make_unique<afw::image::MaskFitsReader>(&_fits);
}

afw::table::BaseCatalog FitsReadHelper::read_catalog(int hdu) {
    _fits.setHdu(hdu, false);
    return afw::table::BaseCatalog::readFits(_fits);
}

afw::table::io::InputArchive FitsReadHelper::read_archive(int hdu) {
    _fits.setHdu(hdu, false);
    return afw::table::io::InputArchive::readFits(_fits);
}

}  // namespace cell_coadds
}  // namespace lsst
