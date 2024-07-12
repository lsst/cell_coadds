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

__all__ = ("ExplodedCoadd",)

from typing import TYPE_CHECKING, Self

import lsst.shoefits as shf

from ._image_planes import ImagePlanes
from ._uniform_grid import UniformGrid

if TYPE_CHECKING:
    from ._multiple_cell_coadd import MultipleCellCoadd


class ExplodedCoadd(ImagePlanes):
    """A coadd that stitches together the outer regions of the cells in a
    `MultipleCellCoadd` (including multiple values for most pixels).

    Notes
    -----
    `ExplodedCoadd` is intended primarily for serialization and visualization,
    and with this in mind it stores its PSF as a single stitched image as well.

    An `ExplodedCoadd` and a `MultipleCellCoadd` may share data, but generally
    only if the `MultipleCellCoadd` is constructed from the `ExplodedCoadd`,
    not the other way around.
    """

    grid: UniformGrid
    """Object that defines the piecewise grid (of outer cell regions) that
    this object stitches together.

    This grid always starts at ``(0, 0)``; because there is no way to align
    this "exploded" grid with the inner cell grid over more than one cell, no
    attempt is made to align it with the inner cell grid's overall offset.
    """

    psf_grid: UniformGrid
    """Object that describes the grid on which PSF model images are stitched
    together.
    """

    psf_image: shf.Image
    """A stitched-together image of the PSF models for each cell."""

    @classmethod
    def from_cell_coadd(cls, cell_coadd: MultipleCellCoadd, pad_psfs_with: float | None) -> Self:
        """Build an `ExplodedCoadd` from a `MultipleCellCoadd`.

        Parameters
        ----------
        cell_coadd : `MultipleCellCoadd`
            Cell-based coadd to stitch together.
        pad_psfs_with : `float` or `None`, optional
            A floating-point value to pad PSF images with so each PSF-image
            cell has the same dimensions as the image (outer) cell it
            corresponds to. If `None`, PSF images will not be padded and the
            full PSF image will generally be smaller than the exploded image it
            corresponds to.

        Returns
        -------
        exploded : `ExplodedCoadd`
            Stitched outer-cell coadd.
        """
        raise NotImplementedError("TODO DM-45189")
