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

__all__ = (
    "CellCoaddFitsFormatter",
    "CellCoaddFitsReader",
    "writeMultipleCellCoaddAsFits",
)

import os
from typing import Any, Mapping

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import numpy as np
from astropy.io import fits
from lsst.afw.image import ImageD, ImageF
from lsst.daf.base import PropertySet
from lsst.geom import Box2I, Extent2I, Point2I
from lsst.obs.base.formatters.fitsGeneric import FitsGenericFormatter
from lsst.skymap import Index2D

from ._common_components import CoaddUnits, CommonComponents
from ._identifiers import CellIdentifiers, PatchIdentifiers
from ._image_planes import OwnedImagePlanes
from ._multiple_cell_coadd import MultipleCellCoadd, SingleCellCoadd
from ._uniform_grid import UniformGrid


class CellCoaddFitsFormatter(FitsGenericFormatter):
    """Interface for writing and reading cell coadds to/from FITS files.

    This assumes the existence of readFits and writeFits methods (for now).
    """


class CellCoaddFitsReader:
    """A reader class to read from a FITS file and produce cell-based coadds.

    This reader class has read methods that can either return a single
    component without reading the entire file (e.g., readBBox, readWcs)
    and read methods that return a full coadd (e.g.,
    readAsMultipleCellCoadd, readAsExplodedCellCoadd, readAsStitchedCoadd).

    Parameters
    ----------
    filename : `str`
        The name of the FITS file to read.
    """

    def __init__(self, filename: str) -> None:
        if not os.path.exists(filename):
            raise FileNotFoundError(f"File {filename} not found")

        self.filename = filename

    def readAsMultipleCellCoadd(self) -> MultipleCellCoadd:
        """Read the FITS file as a MultipleCellCoadd object."""
        with fits.open(self.filename) as hdu_list:
            data = hdu_list[1].data
            header = hdu_list[1].header

            # Read in WCS
            ps = PropertySet()
            ps.update(hdu_list[0].header)
            wcs = afwGeom.makeSkyWcs(ps)

            # Build the quantities needed to construct a MultipleCellCoadd.
            common = CommonComponents(
                units=CoaddUnits(1),  # TODO: read from FITS TUNIT1
                wcs=wcs,
                band=header["BAND"],
                identifiers=PatchIdentifiers(
                    skymap=header["SKYMAP"],
                    tract=header["TRACT"],
                    patch=Index2D(x=header["PATCH_X"], y=header["PATCH_Y"]),
                    band=header["BAND"],
                ),
            )

            grid_cell_size = Extent2I(header["GRCELL1"], header["GRCELL2"])  # Inner size of a single cell.
            grid_shape = Extent2I(header["GRSHAPE1"], header["GRSHAPE2"])
            grid_min = Point2I(header["GRMIN1"], header["GRMIN2"])
            grid = UniformGrid(cell_size=grid_cell_size, shape=grid_shape, min=grid_min)

            # This is the inner bounding box for the multiple cell coadd
            inner_bbox = Box2I(
                Point2I(header["INBBOX11"], header["INBBOX12"]),
                Point2I(header["INBBOX21"], header["INBBOX22"]),
            )

            outer_cell_size = Extent2I(header["OCELL1"], header["OCELL2"])
            psf_image_size = Extent2I(header["PSFSIZE1"], header["PSFSIZE2"])

            coadd = MultipleCellCoadd(
                (
                    self._readSingleCellCoadd(
                        data=row,
                        header=header,
                        common=common,
                        outer_cell_size=outer_cell_size,
                        psf_image_size=psf_image_size,
                        inner_cell_size=grid_cell_size,
                    )
                    for row in data
                ),
                grid=grid,
                outer_cell_size=outer_cell_size,
                psf_image_size=psf_image_size,
                inner_bbox=inner_bbox,
                common=common,
            )

        return coadd

    @staticmethod
    def _readSingleCellCoadd(
        data: Mapping[str, Any],
        common: CommonComponents,
        header: Mapping[str, Any],
        *,
        outer_cell_size: Extent2I,
        inner_cell_size: Extent2I,
        psf_image_size: Extent2I,
    ) -> SingleCellCoadd:
        """Read a coadd from a FITS file.

        Parameters
        ----------
        data : `Mapping`
            The data from the FITS file. Usually, a single row from the binary
            table representation.
        common : `CommonComponents`
            The common components of the coadd.
        outer_cell_size : `Extent2I`
            The size of the outer cell.
        psf_image_size : `Extent2I`
            The size of the PSF image.
        inner_cell_size : `Extent2I`
            The size of the inner cell.

        Returns
        -------
        coadd : `SingleCellCoadd`
            The coadd read from the file.
        """
        buffer = (outer_cell_size - inner_cell_size) // 2

        psf = ImageD(
            array=data["psf"].astype(np.float64),
            xy0=(-(psf_image_size // 2)).asPoint(),  # integer division and negation do not commute.
        )  # use the variable
        xy0 = Point2I(
            inner_cell_size.x * data["cell_id"][0] - buffer.x + header["GRMIN1"],
            inner_cell_size.y * data["cell_id"][1] - buffer.y + header["GRMIN2"],
        )
        mask = afwImage.Mask(data["mask"].astype(np.int32), xy0=xy0)
        image_planes = OwnedImagePlanes(
            image=ImageF(
                data["image"].astype(np.float32),
                xy0=xy0,
            ),
            mask=mask,
            variance=ImageF(data["variance"].astype(np.float32), xy0=xy0),
            noise_realizations=[],
            mask_fractions={},
        )

        identifiers = CellIdentifiers(
            cell=Index2D(data["cell_id"][0], data["cell_id"][1]),
            skymap=common.identifiers.skymap,
            tract=common.identifiers.tract,
            patch=common.identifiers.patch,
            band=common.identifiers.band,
        )

        return SingleCellCoadd(
            outer=image_planes,
            psf=psf,
            inner_bbox=Box2I(
                corner=Point2I(
                    inner_cell_size.x * data["cell_id"][0] + header["GRMIN1"],
                    inner_cell_size.y * data["cell_id"][1] + header["GRMIN2"],
                ),
                dimensions=inner_cell_size,
            ),
            common=common,
            identifiers=identifiers,
            inputs=None,  # type: ignore[arg-type]  # TODO: Pass a sensible value here.
        )

    # TODO: Make this reader a context manager that handles file closures.

    def readWcs(self) -> afwGeom.SkyWcs:
        """Read the WCS information from the FITS file.

        Returns
        -------
        wcs : `~lsst.afw.geom.SkyWcs`
            The WCS information read from the FITS file.
        """
        # Read in WCS
        ps = PropertySet()
        with fits.open(self.filename) as hdu_list:
            ps.update(hdu_list[0].header)
        wcs = afwGeom.makeSkyWcs(ps)
        return wcs


def writeMultipleCellCoaddAsFits(
    multiple_cell_coadd: MultipleCellCoadd,
    filename: str,
    overwrite: bool = False,
    metadata: PropertySet | None = None,
) -> None:
    """Write a MultipleCellCoadd object to a FITS file.

    Parameters
    ----------
    multiple_cell_coadd : `MultipleCellCoadd`
        The multiple cell coadd to write to a FITS file.
    filename : `str`
        The name of the file to write to.
    overwrite : `bool`, optional
        Whether to overwrite the file if it already exists?
    metadata : `~lsst.daf.base.PropertySet`, optional
        Additional metadata to write to the FITS file.
    """
    cell_id = fits.Column(
        name="cell_id",
        format="2I",
        array=[cell.identifiers.cell for cell in multiple_cell_coadd.cells.values()],
    )

    image_array = [cell.outer.image.array for cell in multiple_cell_coadd.cells.values()]
    unit_array = [cell.common.units.name for cell in multiple_cell_coadd.cells.values()]
    image = fits.Column(
        name="image",
        unit=unit_array[0],
        format=f"{image_array[0].size}E",
        dim=f"({image_array[0].shape[1]}, {image_array[0].shape[0]})",
        array=image_array,
    )

    mask_array = [cell.outer.mask.array for cell in multiple_cell_coadd.cells.values()]
    mask = fits.Column(
        name="mask",
        format=f"{mask_array[0].size}I",
        dim=f"({mask_array[0].shape[1]}, {mask_array[0].shape[0]})",
        array=mask_array,
    )

    variance_array = [cell.outer.variance.array for cell in multiple_cell_coadd.cells.values()]
    variance = fits.Column(
        name="variance",
        format=f"{variance_array[0].size}E",
        dim=f"({variance_array[0].shape[1]}, {variance_array[0].shape[0]})",
        array=variance_array,
    )

    psf_array = [cell.psf_image.array for cell in multiple_cell_coadd.cells.values()]
    psf = fits.Column(
        name="psf",
        format=f"{psf_array[0].size}D",
        dim=f"({psf_array[0].shape[1]}, {psf_array[0].shape[0]})",
        array=[cell.psf_image.array for cell in multiple_cell_coadd.cells.values()],
    )

    col_defs = fits.ColDefs([cell_id, image, mask, variance, psf])
    hdu = fits.BinTableHDU.from_columns(col_defs)

    grid_cell_size = multiple_cell_coadd.grid.cell_size
    grid_shape = multiple_cell_coadd.grid.shape
    grid_min = multiple_cell_coadd.grid.bbox.getMin()
    grid_cards = {
        "GRCELL1": grid_cell_size.x,
        "GRCELL2": grid_cell_size.y,
        "GRSHAPE1": grid_shape.x,
        "GRSHAPE2": grid_shape.y,
        "GRMIN1": grid_min.x,
        "GRMIN2": grid_min.y,
    }
    hdu.header.extend(grid_cards)

    outer_cell_size_cards = {
        "OCELL1": multiple_cell_coadd.outer_cell_size.x,
        "OCELL2": multiple_cell_coadd.outer_cell_size.y,
    }
    hdu.header.extend(outer_cell_size_cards)

    psf_image_size_cards = {
        "PSFSIZE1": multiple_cell_coadd.psf_image_size.x,
        "PSFSIZE2": multiple_cell_coadd.psf_image_size.y,
    }
    hdu.header.extend(psf_image_size_cards)

    inner_bbox_cards = {
        "INBBOX11": multiple_cell_coadd.inner_bbox.minX,
        "INBBOX12": multiple_cell_coadd.inner_bbox.minY,
        "INBBOX21": multiple_cell_coadd.inner_bbox.maxX,
        "INBBOX22": multiple_cell_coadd.inner_bbox.maxY,
    }
    hdu.header.extend(inner_bbox_cards)

    wcs = multiple_cell_coadd.common.wcs
    wcs_cards = wcs.getFitsMetadata().toDict()
    primary_hdu = fits.PrimaryHDU()
    primary_hdu.header.extend(wcs_cards)

    hdu.header["TUNIT1"] = multiple_cell_coadd.common.units.name
    # This assumed to be the same as multiple_cell_coadd.common.identifers.band
    # See DM-38843.
    hdu.header["BAND"] = multiple_cell_coadd.common.band
    hdu.header["SKYMAP"] = multiple_cell_coadd.common.identifiers.skymap
    hdu.header["TRACT"] = multiple_cell_coadd.common.identifiers.tract
    hdu.header["PATCH_X"] = multiple_cell_coadd.common.identifiers.patch.x
    hdu.header["PATCH_Y"] = multiple_cell_coadd.common.identifiers.patch.y

    if metadata is not None:
        hdu.header.extend(metadata.toDict())

    hdu_list = fits.HDUList([primary_hdu, hdu])
    hdu_list.writeto(filename, overwrite=overwrite)
