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
import logging
from collections.abc import Iterable, Mapping

import logging
from collections.abc import Iterable, Mapping
from typing import overload

import numpy as np

import lsst.geom as geom
from lsst.skymap import Index2D

from ._uniform_grid import UniformGrid

__all__ = ("StitchedApertureCorrection",)

logger = logging.getLogger(__name__)


class StitchedApertureCorrection:
    """A class that quacks like BoundedField and defined on a grid piecewise.

    Parameters
    ----------
    ugrid : `~lsst.cell_coadds.UniformGrid`
        The uniform grid that defines the bounding boxes for each piece.
    gc : `Mapping[Index2D, float]`
        The grid container that holds the values for each cell.
    """

    def __init__(self, ugrid: UniformGrid, gc: Mapping[Index2D, float]):
        self.ugrid = ugrid
        self.gc = gc

    @overload
    def evaluate(self, x: geom.Point2D | geom.Point2I, y: None) -> float: ...  # noqa: E704

    @overload
    def evaluate(self, x: Iterable[float], y: Iterable[float]) -> np.ndarray: ...  # noqa: E704

    def evaluate(self, x, y=None):  # type: ignore
        """Evaluate the BoundedField at a given point on the image.

        Parameters
        ----------
        x : `~lsst.geom.Point2D` or `~lsst.geom.Point2I`, or `Iterable[float]`
            The point at which to evaluate the BoundedField, or an iterable of
            x-coordinates.
        y : `Iterable[float]`, optional
            The y-coordinates of the points at which to evaluate the aperture
            correction.
            If `None`, `x` is assumed to be a single point.

        Returns
        -------
        value: `float`, or `numpy.ndarray`
            The value of the BoundedField at the specified point.
        """
        if y is None:
            eval_point = geom.Point2I(x)
            idx = self.ugrid.index(eval_point)
            try:
                return self.gc[idx]
            except KeyError:
                logger.info("No aperture correction found for %s", idx)
                return 1.0

        else:
            return np.array([self.evaluate(geom.Point2I(xx, yy)) for xx, yy in zip(x, y, strict=True)])
