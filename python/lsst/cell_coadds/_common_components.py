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

__all__ = ("CoaddUnits", "CommonComponents", "CommonComponentsProperties")

import enum
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from lsst.afw.geom import SkyWcs

    from ._identifiers import PatchIdentifiers


class CoaddUnits(enum.Enum):
    """Enumeration of units a coadd's pixels may have.

    Notes
    -----
    This list is intentionally limited to the units we know we need rather than
    units we think somewhat might want (which is also why we do not support any
    flux unit other than nJy).
    """

    nJy = enum.auto()
    """Pixels represent flux in nanojanskies.
    """

    chi = enum.auto()
    """Pixels represent a signal-to-noise ratio.
    """


@dataclass(frozen=True, eq=False, repr=False)
class CommonComponents:
    """Struct containing image attributes that are common to all cells in a
    patch.
    """

    units: CoaddUnits
    """Units of the coadd's data pixels.
    """

    wcs: SkyWcs
    """World Coordinate System object that maps the pixel grid to sky
    coordinates.
    """

    band: str | None
    """String label for the filter bandpass.

    May be `None` only for coadds that represent a combination of multiple
    passbands (e.g. chi^2 detection coadds), not just to indicate absence of
    knowledge.
    """

    identifiers: PatchIdentifiers
    """Struct of unique identifiers for this coadd's patch."""


class CommonComponentsProperties(ABC):
    """A mix-in class that provides properties for all common components to
    any derived class that provides a property for the common components struct
    itself.
    """

    @property
    def units(self) -> CoaddUnits:
        """Units of the coadd's data pixels."""
        return self.common.units

    @property
    def wcs(self) -> SkyWcs:
        """World Coordinate System object that maps the pixel grid to sky
        coordinates.
        """
        return self.common.wcs

    @property
    def band(self) -> str | None:
        """String label for the filter bandpass.

        May be `None` only for coadds that represent a combination of multiple
        passbands (e.g. chi^2 detection coadds), not just to indicate absence
        of knowledge.
        """
        return self.common.band

    @property
    def identifiers(self) -> PatchIdentifiers:
        """Struct of unique identifiers for this coadd's patch."""
        return self.common.identifiers

    @property
    @abstractmethod
    def common(self) -> CommonComponents:
        """Struct of image components common to all cells in a patch."""
        raise NotImplementedError()
