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
#ifndef LSST_CELL_COADDS_StitchedPsf_h
#define LSST_CELL_COADDS_StitchedPsf_h

#include <cstdint>
#include <memory>
#include <vector>

#include "lsst/meas/algorithms/ImagePsf.h"

namespace lsst {
namespace cell_coadds {

/**
 * A piecewise Psf implementation backed by a 2-d grid of images.
 */
class StitchedPsf : public meas::algorithms::ImagePsf {
public:
    /**
     * Construct by taking ownership of a strided vector of images.
     *
     * @param images   PSF images, with index in the vector given by
     *                 `row*n_cols + col`.  Images need not have the same
     *                 dimensions, but must have odd dimensions.
     * @param bbox     Bounding box of the complete piecewise PSF.
     * @param strides  Distances in pixels between cells.  These must subdivide
     *                 the width and height of `bbox` exactly, yielding the
     *                 number of rows and columns, which must in term multiply
     *                 to equal the length of `images`.
     */
    StitchedPsf(std::vector<std::shared_ptr<afw::detection::Psf::Image>> images, geom::Box2I const& bbox,
                geom::Extent2I const& strides);

    StitchedPsf(StitchedPsf const&) = default;
    StitchedPsf(StitchedPsf&&) = default;

    std::shared_ptr<afw::detection::Psf> clone() const override {
        return std::make_shared<StitchedPsf>(*this);
    }

    std::shared_ptr<afw::detection::Psf> resized(int width, int height) const override;

protected:
    std::shared_ptr<afw::detection::Psf::Image> doComputeKernelImage(
            geom::Point2D const& position, afw::image::Color const& color) const override;
    geom::Box2I doComputeBBox(geom::Point2D const& position, afw::image::Color const& color) const override;

private:
    // All reasoning about how to transform the flat vector of images into
    // a grid of images associated with cells lives here.
    std::shared_ptr<afw::detection::Psf::Image> _image_at(geom::Point2I const& position) const;

    // Indexes are `row*n_cols + col`, but use _image_at instead of encoding
    // that convention elsehwere.
    //
    // We use shared_ptr here because that's what appears in the interfaces we
    // inherit from Psf.  Using shared_ptr there was really a mistake - Image
    // already holds its pixels via an internal shared_ptr, so passing an Image
    // around is already cheap, and adding more indirection is probably a
    // slight pessimization rather than the optimization it was intended to be.
    // But it's not worth fighting that here.
    std::vector<std::shared_ptr<afw::detection::Psf::Image>> _images;

    // Bounding box of full piecewise PSF.  Has nothing to do with the bounding
    // boxes of the elements of _images, but must satisfy
    //   _bbox.getWidth() = _strides.getX()*_n_cols
    //   _bbox.getHeight() = _strides.getY()*_n_rows
    geom::Box2I _bbox;

    // Offsets between cells that define the piecewise behavior.
    geom::Extent2I _strides;

    // Number of cells in each dimension.
    int _n_rows;
    int _n_cols;
};

}  // namespace cell_coadds
}  // namespace lsst

#endif  // !LSST_CELL_COADDS_StitchedPsf_h
