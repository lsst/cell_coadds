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
from collections.abc import Generator, Iterable, Mapping

import lsst.geom as geom
from lsst.skymap import Index2D

from ._identifiers import ObservationIdentifiers
from ._uniform_grid import UniformGrid

__all__ = ("StitchedExposureCatalog",)

logger = logging.getLogger(__name__)


class StitchedExposureCatalog:
    """A class that quacks like an `lsst.afw.table.ExposureCatalog`.

    Parameters
    ----------
    ugrid : `~lsst.cell_coadds.UniformGrid`
        The uniform grid that defines the bounding boxes for each piece.
    gc : `Mapping[Index2D, tuple[ObservationIdentifiers, ...]`
        The grid container that holds the values for each cell.
    """

    def __init__(self, ugrid: UniformGrid, gc: Mapping[Index2D, tuple[ObservationIdentifiers, ...]]):
        self.ugrid = ugrid
        self.gc = gc
        self._inputs: list[ObservationIdentifiers] = []

    def __iter__(self) -> Generator[ObservationIdentifiers]:
        """Iterate over the unique set of `ObservationIdentifiers`
        corresponding to the visits that overlap the full bounding box.
        """
        if not self._inputs:
            inputs: list = []
            for idx in self.gc:
                inputs += self.gc[idx]
            self._inputs = sorted(set(inputs))
        yield from self._inputs

    def __len__(self) -> int:
        if not self._inputs:
            return len([self.__iter__()])
        return len(self._inputs)

    def subsetContaining(self, point: geom.Point2D) -> Iterable[ObservationIdentifiers]:
        """Obtain the set of `ObservationIdentifiers` corresponding to the
        visits at a given point on the image.

        Parameters
        ----------
        point : `~lsst.geom.Point2D` or `~lsst.geom.Point2I`
            The point at which to evaluate the BoundedField.

        Returns
        -------
        value: `~collections.abc.Iterable[ObservationIdentifiers]`
            The set of visits overlapping a given point on the image.
        """
        eval_point = geom.Point2I(point)
        idx = self.ugrid.index(eval_point)
        try:
            return self.gc[idx]
        except KeyError:
            logger.info("No inputs found for %s", idx)
            return ()
