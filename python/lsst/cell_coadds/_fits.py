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

"""Module to handle FITS serialization and de-serialization.

The routines to write and read the files are in the same module, as a change to
one is typically accompanied by a corresponding change to another. Code changes
relating to writing the file must bump to the version number denoted by the
module constant FILE_FORMAT_VERSION.

Although the typical use case is for newer versions of the code to read files
written by an older version, for the purposes of deciding the newer version
string, it is helpful to think about an older version of the reader attempting
to read a newer version of the file on disk. The policy for bumping the version
is as follows:

1. When the on-disk file format written by this module changes such that the
previous version of the reader can still read files written by the newer
version, then there should be a minor bump.

2. When the on-disk format written by this module changes in a way that will
prevent the previous version of the reader from reading a file produced by the
current version of the module, then there should be a major bump. This usually
means that the new version of the reader cannot read older file either,
save the temporary support with deprecation warnings, possibly until a new
release of the Science Pipelines is made.

Examples
--------
1. A file with VERSION=1.3 should still be readable by the reader in
this module when the module-level constant FILE_FORMAT_VERSION=1.4. A file
written with VERSION=1.4 will typically be readable by a reader when the
module-level FILE_FORMAT_VERSION=1.3, although such a use case is not expected.
A concrete example of change
that requires only a minor bump is adding another BinTable that keeps track of
the input visits.

2. An example of major change would be migrating from using
BinTableHDU to ImageHDU to save data. Even if the reader supports reading
either of this formats based on the value of VERSION from the header, it should
be a major change because the previous version of the reader cannot read data
from ImageHDUs.

Unit tests only check that a file written can be read by the concurrent version
of the module, but not by any of the previous ones. Hence, bumping
FILE_FORMAT_VERSION to the appropriate value is ultimately at the discretion of
the developers.

A major bump must also be recorded in the `isCompatibleWith` method.
It is plausible that different (non-consequent) major format versions can be
read by the same reader (due to reverting back to an earlier format, or to
something very similar). `isCompatibleWith` method offers the convenience of
checking if a particular format version can be read by the current reader.

Note that major version 0 is considered unstable and experimental and none of
the guarantee above applies.
"""

from __future__ import annotations

__all__ = (
    "CellCoaddFitsFormatter",
    "CellCoaddFitsReader",
    "IncompatibleVersionError",
    "writeMultipleCellCoaddAsFits",
)

import logging
import os
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from typing import Any

import lsst.shoefits as shf
import numpy as np
from astropy.io import fits
from lsst.obs.base.formatters.fitsGeneric import FitsGenericFormatter
from lsst.skymap import Index2D
from packaging import version

from ._common_components import CoaddUnits, CommonComponents
from ._grid_container import GridContainer
from ._identifiers import CellIdentifiers, ObservationIdentifiers, PatchIdentifiers
from ._image_planes import ImagePlanes
from ._multiple_cell_coadd import MultipleCellCoadd, SingleCellCoadd
from ._to_upstream import CellIndex, CellShape, PixelIndex, PixelShape
from ._uniform_grid import UniformGrid

FILE_FORMAT_VERSION = "0.3"
"""Version number for the file format as persisted, presented as a string of
the form M.m, where M is the major version, m is the minor version.
"""

logger = logging.getLogger(__name__)


class IncompatibleVersionError(RuntimeError):
    """Exception raised when the CellCoaddFitsReader version is not compatible
    with the FITS file attempted to read.
    """


@dataclass
class VisitRecord:
    """A dataclass to hold relevant info about a visit.

    This is intended for use with this module.
    """

    visit: int
    day_obs: int
    physical_filter: str


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

    # Minimum and maximum compatible file format versions are listed as
    # iterables so as to allow for discontiguous intervals.
    MINIMUM_FILE_FORMAT_VERSIONS = ("0.1",)
    MAXIMUM_FILE_FORMAT_VERSIONS = ("1.0",)

    def __init__(self, filename: str) -> None:
        if not os.path.exists(filename):
            raise FileNotFoundError(f"File {filename} not found")

        self.filename = filename

    @classmethod
    def isCompatibleWith(cls, written_version: str, /) -> bool:
        """Check if the serialization version is compatible with the reader.

        This is a convenience method to ask if the current version of this
        class can read a file, based on the VERSION in its header.

        Parameters
        ----------
        written_version: `str`
            The VERSION of the file to be read.

        Returns
        -------
        compatible : `bool`
            Whether the reader can read a file whose VERSION is
            ``written_version``.

        Notes
        -----
        This accepts the other version as a positional argument only.
        """
        written_version_object = version.parse(written_version)
        for min_version, max_version in zip(
            cls.MINIMUM_FILE_FORMAT_VERSIONS,
            cls.MAXIMUM_FILE_FORMAT_VERSIONS,
            strict=True,
        ):
            if version.parse(min_version) <= written_version_object < version.parse(max_version):
                return True

        return False

    def readAsMultipleCellCoadd(self) -> MultipleCellCoadd:
        """Read the FITS file as a MultipleCellCoadd object.

        Raises
        ------
        IncompatibleError
            Raised if the version of this module that wrote the file is
            incompatible with this module that is reading it in.
        """
        with fits.open(self.filename) as hdu_list:
            header = hdu_list[1].header
            written_version = header.get("VERSION", "0.1")
            if not self.isCompatibleWith(written_version):
                raise IncompatibleVersionError(
                    f"{self.filename} was written with version {written_version}"
                    f"but attempting to read it with a reader designed for {FILE_FORMAT_VERSION}"
                )
            if written_version != FILE_FORMAT_VERSION:
                logger.info(
                    "Reading %s having version %s with reader designed for %s",
                    self.filename,
                    written_version,
                    FILE_FORMAT_VERSION,
                )

            written_version = version.parse(written_version)

            # TODO: Remove this when FILE_FORMAT_VERSION is bumped to 1.0
            if written_version < version.parse("0.3"):
                header.rename_keyword("BAND", "FILTER")

            data = hdu_list[1].data

            # TODO: DM-45189: read WCS from header
            wcs = None

            # Build the quantities needed to construct a MultipleCellCoadd.
            common = CommonComponents(
                units=CoaddUnits(header["TUNIT1"]),
                wcs=wcs,
                band=header["FILTER"],
                identifiers=PatchIdentifiers(
                    skymap=header["SKYMAP"],
                    tract=header["TRACT"],
                    patch=Index2D(x=header["PATCH_X"], y=header["PATCH_Y"]),
                    band=header["FILTER"],
                ),
            )

            # Inner size of a single cell.
            grid_cell_size = PixelShape(x=header["GRCELL1"], y=header["GRCELL2"])
            # Number of cells in the grid in each dimension.
            grid_shape = CellShape(x=header["GRSHAPE1"], y=header["GRSHAPE2"])
            # Coordinate of first pixel in the unpadded grid.
            grid_min = PixelIndex(x=header["GRMIN1"], y=header["GRMIN2"])
            # This format version doesn't save the grid padding directly; it's
            # derivable from the inner cell size and the outer cell size.
            outer_cell_size = PixelShape(x=header["OCELL1"], y=header["OCELL2"])
            padding_x, bad_remainder_x = divmod(outer_cell_size.x - grid_cell_size.x, 2)
            padding_y, bad_remainder_y = divmod(outer_cell_size.y - grid_cell_size.y, 2)
            assert not bad_remainder_x, "Cell padding in x cannot be symmetric."
            assert not bad_remainder_y, "Cell padding in y cannot be symmetric."
            assert padding_x == padding_y, "Cell padding should be the same in x and y."
            grid = UniformGrid.from_cell_size_shape(
                grid_cell_size, grid_shape, min=grid_min, padding=padding_x
            )

            inner_bbox = shf.Box.factory[
                header["INBBOX12"] : header["INBBOX22"],
                header["INBBOX11"] : header["INBBOX21"],
            ]
            assert inner_bbox == grid.bbox, "INBBOX and GR* keys are actually redundant."

            psf_image_size = PixelShape(x=header["PSFSIZE1"], y=header["PSFSIZE2"])

            # Attempt to get inputs for each cell.
            inputs = GridContainer[list[ObservationIdentifiers]](shape=grid.shape)
            if written_version >= version.parse("0.3"):
                visit_dict = {
                    row["visit"]: VisitRecord(
                        visit=row["visit"],
                        physical_filter=row["physical_filter"],
                        day_obs=row["day_obs"],
                    )
                    for row in hdu_list[hdu_list.index_of("VISIT")].data
                }
                link_table = hdu_list[hdu_list.index_of("CELL")].data
                for link_row in link_table:
                    cell_id = Index2D(link_row["cell_x"], link_row["cell_y"])
                    visit = link_row["visit"]
                    obs_id = ObservationIdentifiers(
                        instrument=header["INSTRUME"],
                        visit=visit,
                        detector=link_row["detector"],
                        day_obs=visit_dict[visit].day_obs,
                        physical_filter=visit_dict[visit].physical_filter,
                    )
                    if cell_id in inputs:
                        inputs[cell_id] += [obs_id]
                    else:
                        inputs[cell_id] = [obs_id]
            else:
                logger.info(
                    "Cell inputs are available for VERSION=0.3 or later. The file provided has ",
                    "VERSION = %s",
                    written_version,
                )

            # TODO DM-45189: We don't seem to be saving the mask plane
            # definitions at all, so for now we just set them all to be
            # undefined.  Should probably hard-code the afw global default mask
            # plane here instead.
            mask_schema = shf.MaskSchema([None] * 32, dtype=np.uint8)

            coadd = MultipleCellCoadd(
                (
                    self._readSingleCellCoadd(
                        data=row,
                        common=common,
                        inputs=inputs[CellIndex(x=row["cell_id"][0], y=row["cell_id"][1])],
                        grid=grid,
                        outer_cell_size=outer_cell_size,
                        psf_image_size=psf_image_size,
                        mask_schema=mask_schema,
                    )
                    for row in data
                ),
                grid=grid,
                outer_cell_size=outer_cell_size,
                psf_image_size=psf_image_size,
                common=common,
                mask_schema=mask_schema,
            )

        return coadd

    @staticmethod
    def _readSingleCellCoadd(
        data: Mapping[str, Any],
        common: CommonComponents,
        *,
        inputs: Iterable[ObservationIdentifiers],
        grid: UniformGrid,
        outer_cell_size: PixelShape,
        psf_image_size: PixelShape,
        mask_schema: shf.MaskSchema,
    ) -> SingleCellCoadd:
        """Read a coadd from a FITS file.

        Parameters
        ----------
        data : `Mapping`
            The data from the FITS file. Usually, a single row from the binary
            table representation.
        common : `CommonComponents`
            The common components of the coadd.
        inputs : `Iterable` [`ObservationIdentifiers`]
            Any iterable of ObservationIdentifiers instances that contributed
            to this cell.
        grid : `UniformGrid`
            Cell geometry grid.
        outer_cell_size : `PixelShape`
            The size of the outer cell.
        psf_image_size : `PixelShape`
            The size of the PSF image.
        mask_schema : `lsst.shoefits.MaskSchema`
            Definitions of coadded mask planes.

        Returns
        -------
        coadd : `SingleCellCoadd`
            The coadd read from the file.
        """
        psf = shf.Image(
            data["psf"].astype(np.float64),
            # integer division and negation do not commute.
            start=(-(psf_image_size.y // 2), -(psf_image_size.x // 2)),
        )
        cell_index = CellIndex.from_xy(data["cell_id"])
        inner_bbox = grid.bbox_of(cell_index)
        outer_start = PixelIndex(x=inner_bbox.x.start - grid.padding, y=inner_bbox.y.start - grid.padding)
        mask_array = np.frombuffer(data["mask"].tobytes(), dtype=np.uint8).reshape(
            *tuple(outer_cell_size) + (4,)
        )
        mask = shf.Mask(mask_array, start=outer_start, schema=mask_schema)
        unit = common.units.to_astropy()
        image_planes = ImagePlanes(
            image=shf.Image(
                data["image"].astype(np.float32),
                start=outer_start,
                unit=unit,
            ),
            mask=mask,
            variance=shf.Image(data["variance"].astype(np.float32), start=outer_start, unit=unit),
            noise_realizations=[],
            mask_fractions=None,
        )

        identifiers = CellIdentifiers(
            cell=CellIndex(x=data["cell_id"][0], y=data["cell_id"][1]),
            skymap=common.identifiers.skymap,
            tract=common.identifiers.tract,
            patch=common.identifiers.patch,
            band=common.identifiers.band,
        )

        return SingleCellCoadd(
            outer=image_planes,
            psf=psf,
            inner_bbox=inner_bbox,
            common=common,
            identifiers=identifiers,
            inputs=inputs,
        )


def writeMultipleCellCoaddAsFits(
    multiple_cell_coadd: MultipleCellCoadd,
    filename: str,
    overwrite: bool = False,
    metadata: fits.Header | None = None,
) -> fits.HDUList:
    """Write a MultipleCellCoadd object to a FITS file.

    Parameters
    ----------
    multiple_cell_coadd : `MultipleCellCoadd`
        The multiple cell coadd to write to a FITS file.
    filename : `str`
        The name of the file to write to.
    overwrite : `bool`, optional
        Whether to overwrite the file if it already exists?
    metadata : `~astropy.io.fits.Header`, optional
        Additional metadata to write to the FITS file.

    Returns
    -------
    hdu_list : `~astropy.io.fits.HDUList`
        The FITS file as an HDUList.

    Notes
    -----
    Changes to this function that modify the way the file is written to disk
    must be accompanied with a change to FILE_FORMAT_VERSION.
    """
    # Create metadata tables:
    # 1. Visit table containing information about the visits.
    # 2. Cell table containing info about the visit+detector for each cell.
    visit_records: list[Any] = []
    cell_records: list[Any] = []
    instrument_set = set()
    for cell_id, single_cell_coadd in multiple_cell_coadd.cells.items():
        for observation_id in single_cell_coadd.inputs:
            visit_records.append(
                (observation_id.visit, observation_id.physical_filter, observation_id.day_obs)
            )
            cell_records.append((cell_id.x, cell_id.y, observation_id.visit, observation_id.detector))
            instrument_set.add(observation_id.instrument)

    assert len(instrument_set) == 1, "All cells must have the same instrument."
    instrument = instrument_set.pop()

    visit_recarray = np.rec.fromrecords(
        recList=sorted(set(visit_records), key=lambda x: x[0]),  # Sort by visit.
        formats=None,  # formats has specified to please mypy. See numpy#26376.
        names=(
            "visit",
            "physical_filter",
            "day_obs",
        ),
    )
    cell_recarray = np.rec.fromrecords(
        recList=cell_records,
        formats=None,  # formats has specified to please mypy. See numpy#26376.
        names=(
            "cell_x",
            "cell_y",
            "visit",
            "detector",
        ),
    )

    visit_hdu = fits.BinTableHDU.from_columns(visit_recarray, name="VISIT")
    cell_hdu = fits.BinTableHDU.from_columns(cell_recarray, name="CELL")

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

    grid = multiple_cell_coadd.grid
    grid_cell_size = grid.cell_size
    grid_shape = grid.shape
    grid_cards = {
        "GRCELL1": grid_cell_size.x,
        "GRCELL2": grid_cell_size.y,
        "GRSHAPE1": grid_shape.x,
        "GRSHAPE2": grid_shape.y,
        "GRMIN1": grid.bbox.x.min,
        "GRMIN2": grid.bbox.y.min,
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
        "INBBOX11": grid.bbox.x.min,
        "INBBOX12": grid.bbox.y.min,
        "INBBOX21": grid.bbox.x.max,
        "INBBOX22": grid.bbox.y.max,
    }
    hdu.header.extend(inner_bbox_cards)

    primary_hdu = fits.PrimaryHDU()
    # TODO DM-45189: save WCS
    hdu.header["VERSION"] = FILE_FORMAT_VERSION
    hdu.header["TUNIT1"] = multiple_cell_coadd.common.units.name
    # This assumed to be the same as multiple_cell_coadd.common.identifers.band
    # See DM-38843.
    hdu.header["INSTRUME"] = instrument
    hdu.header["FILTER"] = multiple_cell_coadd.common.band
    hdu.header["SKYMAP"] = multiple_cell_coadd.common.identifiers.skymap
    hdu.header["TRACT"] = multiple_cell_coadd.common.identifiers.tract
    hdu.header["PATCH_X"] = multiple_cell_coadd.common.identifiers.patch.x
    hdu.header["PATCH_Y"] = multiple_cell_coadd.common.identifiers.patch.y

    if metadata is not None:
        hdu.header.extend(metadata)

    hdu_list = fits.HDUList([primary_hdu, hdu, cell_hdu, visit_hdu])
    hdu_list.writeto(filename, overwrite=overwrite)

    return hdu_list
