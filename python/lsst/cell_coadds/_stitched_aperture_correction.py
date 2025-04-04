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

from collections.abc import Mapping

import lsst.geom as geom
from lsst.skymap import Index2D

from ._uniform_grid import UniformGrid

__all__ = ("StitchedApertureCorrection",)


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

    def evaluate(self, point: geom.Point2D | geom.Point2I) -> float:
        """Evaluate the BoundedField at a given point on the image.

        Parameters
        ----------
        point : `~lsst.geom.Point2D` or `~lsst.geom.Point2I`
            The point at which to evaluate the BoundedField.

        Returns
        -------
        value: `float`
            The value of the BoundedField at the specified point.
        """
        eval_point = geom.Point2I(point)
        idx = self.ugrid.index(eval_point)
        return self.gc[idx]
