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

import lsst.afw.geom as afwGeom
import lsst.geom as geom
from lsst.daf.butler import DataCoordinate, DimensionUniverse

"""Collection of utility functions to be used for unit tests."""


def generate_data_id(
    *,
    tract: int = 9813,
    patch: int = 42,
    cell_x: int = 4,
    cell_y: int = 2,
    band: str = "r",
    detector_id: int = 9,
    visit_id: int = 1234,
    detector_max: int = 109,
    visit_max: int = 10000
) -> DataCoordinate:
    """Generate a DataCoordinate instance to use as data_id.

    Parameters
    ----------
    tract : `int`, optional
        Tract ID for the data_id
    patch : `int`, optional
        Patch ID for the data_id
    cell_x : `int`, optional
        X index of the cell this patch corresponds to
    cell_y : `int`, optional
        Y index of the cell this patch corresponds to
    band : `str`, optional
        Band for the data_id
    detector_id : `int`, optional
        Detector ID for the data_id
    visit_id : `int`, optional
        Visit ID for the data_id
    detector_max : `int`, optional
        Maximum detector ID for the data_id
    visit_max : `int`, optional
        Maximum visit ID for the data_id

    Returns
    -------
    data_id : `lsst.daf.butler.DataCoordinate`
        An expanded data_id instance.
    """
    universe = DimensionUniverse()

    instrument = universe["instrument"]
    instrument_record = instrument.RecordClass(
        name="DummyCam",
        class_name="lsst.obs.base.instrument_tests.DummyCam",
        detector_max=detector_max,
        visit_max=visit_max,
    )

    skymap = universe["skymap"]
    skymap_record = skymap.RecordClass(name="test_skymap")

    band_element = universe["band"]
    band_record = band_element.RecordClass(name=band)

    visit = universe["visit"]
    visit_record = visit.RecordClass(id=visit_id, instrument="test")

    detector = universe["detector"]
    detector_record = detector.RecordClass(id=detector_id, instrument="test")

    physical_filter = universe["physical_filter"]
    physical_filter_record = physical_filter.RecordClass(name=band, instrument="test", band=band)

    patch_element = universe["patch"]
    patch_record = patch_element.RecordClass(
        skymap="test_skymap", tract=tract, patch=patch, cell_x=cell_x, cell_y=cell_y
    )

    # A dictionary with all the relevant records.
    record = dict(
        instrument=instrument_record,
        visit=visit_record,
        detector=detector_record,
        patch=patch_record,
        tract=9813,
        band=band_record.name,
        skymap=skymap_record.name,
        physical_filter=physical_filter_record,
    )

    # A dictionary with all the relevant recordIds.
    record_id = record.copy()
    for key in ("visit", "detector"):
        record_id[key] = record_id[key].id

    # TODO: Catching mypy failures on Github Actions should be made easier,
    # perhaps in DM-36873. Igroring these for now.
    data_id = DataCoordinate.standardize(record_id, universe=universe)
    return data_id.expanded(record)


def generate_wcs(*, scale: float = 0.168, flipX: bool = True) -> afwGeom.SkyWcs:
    """Generate a SkyWcs instant with a given pixel scale.

    Parameters
    ----------
    scale : `float`, optional
        Pixel scale in arcseconds.
    flipX : `bool`, optional
        Flip the X axis.

    Returns
    -------
    wcs : `lsst.afw.geom.SkyWcs`
        A SkyWcs instance.
    """
    orientation = -45 * geom.degrees
    scale = scale * geom.arcseconds
    crpix = geom.Point2D(100, 100)
    crval = geom.SpherePoint(30, 60, geom.degrees)
    cdMatrix = afwGeom.makeCdMatrix(scale=scale, orientation=orientation, flipX=flipX)
    return afwGeom.makeSkyWcs(crpix=crpix, crval=crval, cdMatrix=cdMatrix)
