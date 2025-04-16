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

"""Collection of type hint stubs for use in type checking."""


from __future__ import annotations

from collections.abc import Mapping
from typing import Protocol, Self, overload

import numpy as np
from frozendict import frozendict

from lsst.geom import Box2I, Point2I

from ._stitched_aperture_correction import StitchedApertureCorrection


class ImageLike(Protocol):
    """Interface for objects that behave like `lsst.afw.image.Image` with
    respect to subimage slicing, bounding box access, and XY0 offset.
    """

    @overload
    def __getitem__(self, slices: tuple[slice, slice]) -> Self:
        pass

    @overload
    def __getitem__(self, bbox: Box2I) -> Self:
        pass

    def __setitem__(self, bbox: Box2I, other: Self | int) -> None:  # noqa: D105
        pass

    def getBBox(self) -> Box2I:  # noqa: D102
        pass

    def getXY0(self) -> Point2I:  # noqa: D102
        pass

    def setXY0(self, xy0: Point2I) -> None:  # noqa: D102
        pass

    @property
    def array(self) -> np.ndarray:  # noqa: D102
        pass


SingleCellCoaddApCorrMap: type = frozendict[str, float]
"""A type alias for aperture correction maps for single cell coadds."""

StitchedCoaddApCorrMap: type = Mapping[str, StitchedApertureCorrection]
"""A type alias for aperture correction maps for stitched coadds."""
