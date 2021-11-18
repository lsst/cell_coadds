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

#include "lsst/cell_coadds/GridContainer.h"
#include "lsst/cell_coadds/UniformGrid.h"
#include "lsst/meas/algorithms/ImagePsf.h"

namespace lsst {
namespace cell_coadds {

/**
 * A piecewise Psf implementation backed by a 2-d grid of images.
 */
class StitchedPsf final : public meas::algorithms::ImagePsf {
public:
    /**
     * Construct by taking ownership of a strided vector of images.
     *
     * @param images   PSF images.  Images need not have the same dimensions,
     *                 but must have odd dimensions.
     *
     * @param grid     Object that defines the geometry of the piecewise image
     *                 this PSF corresponds to.
     */
    StitchedPsf(
        GridContainer<std::shared_ptr<afw::detection::Psf::Image>> const& images, UniformGrid const& grid);

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
    // We use shared_ptr here because that's what appears in the interfaces we
    // inherit from Psf.  Using shared_ptr there was really a mistake - Image
    // already holds its pixels via an internal shared_ptr, so passing an Image
    // around is already cheap, and adding more indirection is probably a
    // slight pessimization rather than the optimization it was intended to be.
    // But it's not worth fighting that here.
    GridContainer<std::shared_ptr<afw::detection::Psf::Image>> _images;

    UniformGrid _grid;
};

}  // namespace cell_coadds
}  // namespace lsst

#endif  // !LSST_CELL_COADDS_StitchedPsf_h
