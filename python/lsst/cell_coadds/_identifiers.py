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
    "PatchIdentifiers",
    "CellIdentifiers",
    "ObservationIdentifiers",
)


from dataclasses import dataclass
from typing import Any, Self, cast

from lsst.daf.butler import DataCoordinate, DimensionRecord
from lsst.skymap import Index2D


class BaseIdentifiers:
    """Base class for identifiers. This acts like a mixin."""

    def __getitem__(self, key: str) -> Any:
        """Get an attribute by name.

        Parameters
        ----------
        key : `str`
            Name of the attribute to get.

        Returns
        -------
        value : `int`, `float`, `str`, `~lsst.skymap.Index2D` or `None`
            Value of the attribute.
        """
        return getattr(self, key)


@dataclass(frozen=True)
class PatchIdentifiers(BaseIdentifiers):
    """Struct of identifiers for a coadd patch."""

    skymap: str
    """The name of the skymap this patch belongs to.
    """

    tract: int
    """The name of the tract this patch belongs to.
    """

    patch: Index2D
    """Identifiers for the patch itself.
    """

    band: str | None
    """Name of the band, if any.
    """

    @classmethod
    def from_data_id(cls, data_id: DataCoordinate) -> PatchIdentifiers:
        """Construct from a data ID.

        Parameters
        ----------
        data_id : `~lsst.daf.butler.DataCoordinate`
            Fully-expanded data ID that includes the 'patch' dimension and
            optionally the `band` dimension.

        Returns
        -------
        identifiers : `PatchIdentifiers`
            Struct of identifiers for this patch.
        """
        patch_record = cast(DimensionRecord, data_id.records["patch"])
        return cls(
            skymap=cast(str, data_id["skymap"]),
            tract=cast(int, data_id["tract"]),
            patch=Index2D(x=patch_record.cell_x, y=patch_record.cell_y),
            band=cast(str, data_id.get("band")),
        )


@dataclass(frozen=True)
class CellIdentifiers(PatchIdentifiers):
    """Struct of identifiers for a coadd cell."""

    cell: Index2D
    """Identifiers for the cell itself."""

    @classmethod
    def from_data_id(  # type: ignore [override]
        cls, data_id: DataCoordinate, cell: Index2D
    ) -> CellIdentifiers:
        """Construct from a data ID and a cell index.

        Parameters
        ----------
        data_id : `~lsst.daf.butler.DataCoordinate`
            Fully-expanded data ID that includes the 'patch' dimension and
            optionally the `band` dimension.
        cell : `~lsst.skymap.Index2D`
            Index of the cell within the patch.

        Returns
        -------
        identifiers : `CellIdentifiers`
            Struct of identifiers for this cell within a patch.
        """
        patch_record = cast(DimensionRecord, data_id.records["patch"])
        return cls(
            skymap=cast(str, data_id["skymap"]),
            tract=cast(int, data_id["tract"]),
            patch=Index2D(x=patch_record.cell_x, y=patch_record.cell_y),
            band=cast(str, data_id.get("band")),
            cell=cell,
        )


@dataclass(frozen=True)
class ObservationIdentifiers(BaseIdentifiers):
    """Struct of identifiers for an observation that contributed to a coadd
    cell.
    """

    instrument: str
    """Name of the instrument that this observation was taken with.
    """

    physical_filter: str
    """Name of the physical filter that this observation was taken with.
    """

    visit: int
    """Unique identifier for the visit.

    A visit may be comprised of more than one exposure only if all were
    observed back-to-back with no dithers, allowing them to be combined early
    the processing with no resampling.  All detector-level images in a visit
    share the same visit ID.
    """

    day_obs: int
    """A day and night of observations that rolls over during daylight hours.
      The identifier is an decimal integer-concatenated date, i.e. YYYYMMDD,
      with the exact rollover time observatory-dependent.
    """

    detector: int
    """Unique identifier for the detector.
    """

    @property
    def ccd(self) -> int:
        """Alias for the detector.

        This is provided for compatibility with the older API.
        """
        return self.detector

    @classmethod
    def from_data_id(cls, data_id: DataCoordinate, *, backup_detector: int = -1) -> ObservationIdentifiers:
        """Construct from a data ID.

        Parameters
        ----------
        data_id : `~lsst.daf.butler.DataCoordinate`
            Fully-expanded data ID that includes the 'visit', 'detector' and
            'day_obs' dimensions.
        backup_detector : `int`, optional
            Detector ID to use as a backup if not present in ``data_id``.
            This is not used if detector information is available in
            ``data_id`` and does not override it.

        Returns
        -------
        identifiers : `ObservationIdentifiers`
            Struct of identifiers for this observation.
        """
        detector = data_id.get("detector", backup_detector)
        day_obs = data_id.get("day_obs")
        return cls(
            instrument=cast(str, data_id["instrument"]),
            physical_filter=cast(str, data_id["physical_filter"]),
            visit=cast(int, data_id["visit"]),
            day_obs=cast(int, day_obs),
            detector=detector,
        )

    def __lt__(self, other: Self, /) -> bool:
        return (self.visit, self.detector) < (other.visit, other.detector)
