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

__all__ = ("MultipleCellCoadd",)

from collections.abc import Iterable
from typing import TYPE_CHECKING

import lsst.shoefits as shf

from ._common_components import CommonComponents, CommonComponentsProperties
from ._exploded_coadd import ExplodedCoadd
from ._grid_container import GridContainer
from ._single_cell_coadd import SingleCellCoadd
from ._stitched_coadd import StitchedCoadd
from ._to_upstream import PixelShape
from ._uniform_grid import UniformGrid

if TYPE_CHECKING:
    import astropy.io.fits


class MultipleCellCoadd(CommonComponentsProperties):
    """A data structure for coadds built from many overlapping cells.

    Notes
    -----
    `MultipleCellCoadd` is designed to be used both by measurement algorithms
    that are able to take advantage of cell boundaries and overlap regions
    (which can use the ``.cells`` attribute to access `SingleCellCoadd` objects
    directly) and measurement algorithms that just want one image and don't
    care (or don't care much) about discontinuities (which can use `stitch` to
    obtain such an image).
    """

    def __init__(
        self,
        cells: Iterable[SingleCellCoadd],
        grid: UniformGrid,
        outer_cell_size: PixelShape,
        psf_image_size: PixelShape,
        *,
        common: CommonComponents,
    ):
        self._grid = grid
        self._outer_cell_size = outer_cell_size
        self._psf_image_size = psf_image_size
        self._common = common
        self._cells = GridContainer[SingleCellCoadd](self._grid.shape)

        for cell in cells:
            index = cell.identifiers.cell
            self._cells[index] = cell
            if cell.inner.bbox != self._grid.bbox_of(index):
                raise ValueError(
                    f"Cell at index {index} has inner bbox {cell.inner.bbox}, "
                    f"but grid expects {self._grid.bbox_of(index)}."
                )
            if cell.outer.bbox.shape != self._outer_cell_size:
                raise ValueError(
                    f"Cell at index {index} has outer dimensions {cell.outer.bbox.shape}, "
                    f"but coadd expects {self._outer_cell_size}."
                )
            if cell.psf_image.bbox.shape != self._psf_image_size:
                raise ValueError(
                    f"Cell at index {index} has PSF image with dimensions {cell.psf_image.bbox.shape}, "
                    f"but coadd expects {self._psf_image_size}."
                )

        self._mask_fractions = {cell.outer.mask_fractions is not None for cell in self._cells.values()}
        self._has_mask_fractions = self._mask_fractions.pop()
        if self._has_mask_fractions:
            raise ValueError("Cells are inconsistent in presence/absence of mask fractions.")
        n_noise_realizations = {len(cell.outer.noise_realizations) for cell in self._cells.values()}
        self._n_noise_realizations = n_noise_realizations.pop()
        if n_noise_realizations:
            n_noise_realizations.add(self._n_noise_realizations)
            raise ValueError(
                f"Inconsistent number of noise realizations ({n_noise_realizations}) between cells."
            )

    @property
    def cells(self) -> GridContainer[SingleCellCoadd]:
        """The grid of single-cell coadds, indexed by (y, x)."""
        return self._cells

    @property
    def n_noise_realizations(self) -> int:
        """The number of noise realizations cells are guaranteed to have."""
        return self._n_noise_realizations

    @property
    def has_mask_fractions(self) -> bool:
        """Whether cells have ``mask_fractions`` planes that record the
        fraction of contributing pixels that were interpolated bad pixels.
        """
        return self._has_mask_fractions

    @property
    def grid(self) -> UniformGrid:
        """Object that defines the inner geometry for all cells.

        The padding in this grid reflects the outer cell boundaries of the
        outermost cells in the grid.
        """
        return self._grid

    @property
    def outer_cell_size(self) -> PixelShape:
        """Dimensions of the outer region of each cell."""
        return self._outer_cell_size

    @property
    def psf_image_size(self) -> PixelShape:
        """Dimensions of PSF model images."""
        return self._psf_image_size

    @property
    def outer_bbox(self) -> shf.Box:
        """The rectangular region fully covered by all cell outer bounding
        boxes.

        This is an alias for ``self.grid.bbox_with_padding``.
        """
        return self.grid.bbox_with_padding

    @property
    def inner_bbox(self) -> shf.Box:
        """The rectangular region fully covered by all cell inner bounding
        boxes.

        This is an alias for ``self.grid.bbox``.
        """
        return self.grid.bbox

    @property
    def common(self) -> CommonComponents:
        # Docstring inherited.
        return self._common

    def stitch(self, bbox: shf.Box | None = None) -> StitchedCoadd:
        """Return a contiguous (but in general discontinuous) coadd by
        stitching together inner cells.

        Parameters
        ----------
        bbox : `Box2I`, optional
            Region for the returned coadd; default is ``self.inner_bbox``.

        Returns
        -------
        stitched : `StitchedCellCoadd`
            Contiguous coadd covering the given area.  Each image plane is
            actually constructed when first accessed, not when this method
            is called.
        """
        # In the future, stitching algorithms that apply ramps to smooth
        # discontinuities may also be provided; we'd implement that by having
        # this return different types (from a common ABC), perhaps dispatched
        # by an enum.
        return StitchedCoadd.from_cell_coadd(self, bbox=bbox)

    def explode(self, pad_psfs_with: float | None = None) -> ExplodedCoadd:
        """Return a coadd whose image planes stitch together the outer regions
        of each cell, duplicating pixels in the overlap regions.

        Parameters
        ----------
        pad_psfs_with : `float` or None, optional
            A floating-point value to pad PSF images with so each PSF-image
            cell has the same dimensions as the image (outer) cell it
            corresponds to.  If `None`, PSF images will not be padded and the
            full PSF image will generally be smaller than the exploded image it
            corresponds to.

        Returns
        -------
        exploded : `ExplodedCoadd`
            Exploded version of the coadd.
        """
        return ExplodedCoadd.from_cell_coadd(self, pad_psfs_with=pad_psfs_with)

    @classmethod
    def read_fits(cls, filename: str) -> MultipleCellCoadd:
        """Read a MultipleCellCoadd from a FITS file.

        Parameters
        ----------
        filename : `str`
            The path to the FITS file to read.

        Returns
        -------
        cell_coadd : `MultipleCellCoadd`
            The MultipleCellCoadd object read from the FITS file.
        """
        # Avoid circular import, keep afw dependence optional.
        from ._fits import CellCoaddFitsReader

        reader = CellCoaddFitsReader(filename)
        return reader.readAsMultipleCellCoadd()

    @classmethod
    def readFits(cls, *args, **kwargs) -> MultipleCellCoadd:  # type: ignore[no-untyped-def]
        """Alias to `read_fits` method.

        Notes
        -----
        This method exists for compatability with the rest of the codebase.
        The presence of this method allows for reading in via
        `lsst.obs.base.formatters.FitsGenericFormatter`.
        Whenever possible, use `read_fits` instead, since this method may be
        deprecated in the near future.
        """
        return cls.read_fits(*args, **kwargs)

    def write_fits(
        self, filename: str, overwrite: bool = False, metadata: astropy.io.fits.Header | None = None
    ) -> None:
        """Write the coadd as a FITS file.

        Parameters
        ----------
        filename : `str`
            The path to the FITS file to write.
        overwrite : `bool`, optional
            Whether to overwrite an existing file?
        metadata : `astropy.io.fits.Header`, optional
            Additional metadata to write to the FITS header.
        """
        # Avoid circular import, keep afw dependence optional.
        from ._fits import writeMultipleCellCoaddAsFits

        writeMultipleCellCoaddAsFits(self, filename, overwrite=overwrite, metadata=metadata)

    def writeFits(self, *args, **kwargs) -> None:  # type: ignore[no-untyped-def]
        """Alias to `write_fits` method.

        Notes
        -----
        This method exists for compatability with the rest of the codebase.
        The presence of this method allows for persistence via
        `lsst.obs.base.formatters.FitsGenericFormatter`.
        Whenever possible, use `write_fits` instead, since this method may be
        deprecated in the near future.
        """
        self.write_fits(*args, **kwargs)
