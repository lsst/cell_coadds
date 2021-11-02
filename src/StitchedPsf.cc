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

#include "lsst/pex/exceptions.h"
#include "lsst/cell_coadds/StitchedPsf.h"

namespace lsst {
namespace cell_coadds {

StitchedPsf::StitchedPsf(
    GridContainer<std::shared_ptr<afw::detection::Psf::Image>> const& images, UniformGrid const& grid)
        : _images(images), _grid(grid) {}

std::shared_ptr<afw::detection::Psf> StitchedPsf::resized(int width, int height) const {
    if (width % 2 != 1) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthError,
            (boost::format("resized width must be a positive odd integer; got %d.") % width).str());
    }
    if (height % 2 != 1) {
        throw LSST_EXCEPT(
            pex::exceptions::LengthError,
            (boost::format("resized height must be a positive odd integer; got %d.") % height).str());
    }
    geom::Box2I bbox(geom::Point2I(-width / 2, -height / 2), geom::Point2I(width / 2, height / 2));
    auto new_images =
        GridContainer(_images).rebuild_transformed([bbox](std::shared_ptr<afw::detection::Psf::Image> image) {
            if (image->getBBox().contains(bbox)) {
                return std::make_shared<afw::detection::Psf::Image>(image->subset(bbox));
            } else {
                // Make a new image big enough to fit current bbox and new bbox,
                // copy current image into it, then subset that for the returned
                // PSF.
                afw::detection::Psf::Image bigger_image(bbox.expandedTo(image->getBBox()));
                bigger_image = 0.0f;
                bigger_image.subset(image->getBBox()).assign(*image);
                return std::make_shared<afw::detection::Psf::Image>(bigger_image.subset(bbox));
            }
        });
    return std::make_shared<StitchedPsf>(std::move(new_images).finish(), _grid);
}

std::shared_ptr<afw::detection::Psf::Image> StitchedPsf::doComputeKernelImage(
    geom::Point2D const& position, afw::image::Color const& color) const {
    return _images[_grid.index(geom::Point2I(position))];
}

geom::Box2I StitchedPsf::doComputeBBox(geom::Point2D const& position, afw::image::Color const& color) const {
    return _images[_grid.index(geom::Point2I(position))]->getBBox();
}

}  // namespace cell_coadds
}  // namespace lsst
