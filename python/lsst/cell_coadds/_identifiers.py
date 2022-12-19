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
from typing import TYPE_CHECKING

from lsst.skymap import Index2D

if TYPE_CHECKING:
    from lsst.daf.butler import DataCoordinate


@dataclass(frozen=True)
class PatchIdentifiers:
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
        # The 'type: ignore' directives below are present because we require a
        # a fully-expanded data ID, but this isn't embedded in the types.
        return cls(
            skymap=data_id["skymap"],  # type: ignore
            tract=data_id["tract"],  # type: ignore
            patch=Index2D(
                x=data_id.records["patch"].cell_x, y=data_id.records["patch"].cell_y  # type: ignore
            ),
            band=data_id.get("band"),
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
        return cls(
            skymap=data_id["skymap"],  # type: ignore
            tract=data_id["tract"],  # type: ignore
            patch=Index2D(
                x=data_id.records["patch"].cell_x, y=data_id.records["patch"].cell_y  # type: ignore
            ),
            band=data_id.get("band"),
            cell=cell,
        )


@dataclass(frozen=True)
class ObservationIdentifiers:
    """Struct of identifiers for an observation that contributed to a coadd
    cell."""

    instrument: str
    """Name of the instrument that this observation was taken with.
    """

    packed: int
    """ID that uniquely identifies both the visit and detector by packing
    together their IDs.
    """

    visit: int
    """Unique identifier for the visit.

    A visit may be comprised of more than one exposure only if all were
    observed back-to-back with no dithers, allowing them to be combined early
    the processing with no resampling.  All detector-level images in a visit
    share the same visit ID.
    """

    detector: int
    """Unique identifier for the detector.
    """

    @classmethod
    def from_data_id(cls, data_id: DataCoordinate) -> ObservationIdentifiers:
        """Construct from a data ID.

        Parameters
        ----------
        data_id : `~lsst.daf.butler.DataCoordinate`
            Fully-expanded data ID that includes the 'visit' and 'detector'
            dimensions.

        Returns
        -------
        identifiers : `ObservationIdentifiers`
            Struct of identifiers for this observation.
        """
        # The 'type: ignore' directives below are present because we require a
        # a fully-expanded data ID, but this isn't embedded in the types.
        return cls(
            instrument=data_id["instrument"],  # type: ignore
            packed=data_id.pack("visit_detector"),  # type: ignore
            visit=data_id["visit"],  # type: ignore
            detector=data_id["detector"],  # type: ignore
        )
