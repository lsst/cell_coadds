# This file is part of cell_coadds.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from __future__ import annotations

__all__ = ("StitchedCoadd",)

from typing import TYPE_CHECKING, Self

import lsst.shoefits as shf

from ._common_components import CommonComponents
from ._grid_container import GridContainer
from ._image_planes import ImagePlanes
from ._uniform_grid import UniformGrid

if TYPE_CHECKING:
    from ._multiple_cell_coadd import MultipleCellCoadd


class StitchedCoadd(ImagePlanes, CommonComponents):
    """A coadd that stitches together images from adjacent cells.

    Notes
    -----
    A `StitchedCoadd` coadd never shares image data with a `MultipleCellCoadd`;
    pixel values must always be copied.  PSF models images may be shared,
    however.
    """

    grid: UniformGrid
    """Object that defines the piecewise grid (of inner cell regions) that
    this object stitches together.

    This may include cells outside the region covered by these image planes.
    """

    psf_images: GridContainer[shf.Image]
    """The piecewise PSF of this image as a grid of PSF model images."""

    @classmethod
    def from_cell_coadd(cls, cell_coadd: MultipleCellCoadd, bbox: shf.Box | None = None) -> Self:
        raise NotImplementedError("TODO DM-45189")
