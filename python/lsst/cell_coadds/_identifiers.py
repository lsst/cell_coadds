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
    "GridIdentifiers",
    "PatchIdentifiers",
    "CellIdentifiers",
    "ObservationIdentifiers",
)


from dataclasses import dataclass
from typing import Tuple


@dataclass(frozen=True)
class GridIdentifiers:
    """Struct of identifiers that identify an element in a grid.

    This may be used to represent patches within a tract, or cells within a
    patch.
    """

    packed: int
    """ID that uniquely identifies both this element and grid in which it lives
    by packing together the IDs of each.
    """

    sequential: int
    """Unique sequential integer ID for this element in the grid."""

    x: int
    """The column of this element in its grid.
    """

    y: int
    """The row of this element in its grid.
    """

    @property
    def index(self) -> Tuple[int, int]:
        """The index as a 2-d (x, y) tuple."""
        return (self.x, self.y)


@dataclass(frozen=True)
class PatchIdentifiers:
    """Struct of identifiers for a coadd patch."""

    skymap: str
    """The name of the skymap this patch belongs to.
    """

    tract: int
    """The name of the tract this patch belongs to.
    """

    patch: GridIdentifiers
    """Identifiers for the patch itself.
    """


@dataclass(frozen=True)
class CellIdentifiers(PatchIdentifiers):
    """Struct of identifiers for a coadd cell."""

    cell: GridIdentifiers
    """Identifiers for the cell itself."""


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
