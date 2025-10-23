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

import pickle
from functools import partial
from typing import TYPE_CHECKING, Any, ClassVar

import numpy as np

import lsst.geom as geom
from lsst.afw.detection import InvalidPsfError
from lsst.afw.image import ImageD
from lsst.afw.typehandling import StorableHelperFactory
from lsst.meas.algorithms import ImagePsf

from ._grid_container import GridContainer
from ._uniform_grid import UniformGrid

if TYPE_CHECKING:
    from lsst.afw.image import Color

__all__ = ("StitchedPsf",)


class StitchedPsf(ImagePsf):
    """A piecewise PSF implementation backed by a 2-d grid of images."""

    # We need to ensure a C++ StorableHelperFactory is constructed and
    # available before any unpersists of this class.  Placing this "private"
    # class attribute here accomplishes that.
    _factory: ClassVar[type[StorableHelperFactory]] = StorableHelperFactory("lsst.cell_coadds", "StitchedPsf")

    def __init__(self, images: GridContainer[ImageD], grid: UniformGrid) -> None:
        self._validate_args(images, grid)

        super().__init__()
        self._images = images
        self._grid = grid
        self._averagePosition = None

    @staticmethod
    def _validate_args(images: GridContainer[ImageD], grid: UniformGrid) -> None:
        """Validate the images and grid.

        Parameters
        ----------
        images : `GridContainer`
            The images to validate.
        grid : `UniformGrid`
            The grid to validate.

        Raises
        ------
        ValueError
            Raised if the images and grid are incompatible.
        """
        min_x = min(index.x for index in images.indices())
        min_y = min(index.y for index in images.indices())
        max_x = max(index.x for index in images.indices())
        max_y = max(index.y for index in images.indices())

        if ((max_x - min_x + 1) > grid.shape.x) or ((max_y - min_y + 1) > grid.shape.y):
            raise ValueError("Images do not fit on grid.")

    @property
    def images(self) -> GridContainer[ImageD]:
        """The images that make up this PSF."""
        return self._images

    @property
    def grid(self) -> UniformGrid:
        """The grid on which the images are placed."""
        return self._grid

    def getAveragePosition(self) -> geom.Point2D:
        """Get a position where PSF can be evaluated on a patch.

        This defaults to the center of the patch bounding box, unless there are
        no inputs there. In that case, it switches to find an arbitrary cell,
        typically at a corner that has inputs and returns the center position
        of the cell.
        """
        if self._averagePosition is None:
            center = self._grid.bbox.getCenter()
            if self.grid.index(geom.Point2I(center)) not in self._images:
                arbitrary_index = next(iter(self._images))
                bbox = self._grid.bbox_of(arbitrary_index)
                center = bbox.getCenter()

            self._averagePosition = center

        return self._averagePosition

    # The _do* methods make use of the ImagePsf trampoline.
    def _doComputeBBox(self, position: geom.Point2D | geom.Point2I, color: Color = None) -> geom.Box2I:
        try:
            return self._images[self._grid.index(geom.Point2I(position))].getBBox()
        except (KeyError, ValueError):
            raise InvalidPsfError("No inputs exists at position.") from None

    def _doComputeKernelImage(self, position: geom.Point2D | geom.Point2I, color: Color = None) -> ImageD:
        try:
            return self._images[self._grid.index(geom.Point2I(position))]
        except (KeyError, ValueError):
            raise InvalidPsfError("No inputs exists at position.") from None

    def clone(self) -> StitchedPsf:
        """Return a deep copy of this object."""
        return StitchedPsf(self.images, self.grid)

    def __deepcopy__(self, memo: dict[int, Any] | None = None) -> StitchedPsf:
        """Return a deep copy of this object."""
        return StitchedPsf(self.images, self.grid)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, StitchedPsf):
            return False

        if not (self.grid == other.grid):
            return False

        for index in self.images.indices():
            if (
                not (self.images[index].array.shape == other.images[index].array.shape)
                or not np.equal(self.images[index].array, other.images[index].array).all()
            ):
                return False

        return True

    @staticmethod
    def _callback(image: ImageD, bbox: geom.Box2I) -> ImageD:
        if image.getBBox().contains(bbox):
            return image[bbox]
        else:
            # Make a new image big enough to fit current bbox and new bbox,
            # copy current image into it, then subset that for the returned
            # PSF.
            bigger_image = ImageD(bbox=bbox.expandedTo(image.getBBox()), initialValue=0.0)
            bigger_image[image.getBBox()] = image
            return bigger_image[bbox]

    def resized(self, width: int, height: int) -> StitchedPsf:
        if not (width % 2 == 1 and width > 0):
            raise ValueError("resized width must be a positive odd integer; got {width}.")
        if not (height % 2 == 1 and height > 0):
            raise ValueError("resized height must be a positive odd integer; got {height}.")

        bbox = geom.Box2I(geom.Point2I(-(width // 2), -(height // 2)), geom.Extent2I(width, height))
        gc = self._images.rebuild_transformed(transform=partial(self._callback, bbox=bbox))
        return StitchedPsf(gc, self.grid)

    @staticmethod
    def isPersistable() -> bool:
        return True

    @staticmethod
    def _getPersistenceName() -> str:
        return "StitchedPsf"

    @staticmethod
    def _getPythonModule() -> str:
        return __name__

    # The get/set state methods are needed to support pickle.
    def __getstate__(self) -> dict:
        return {"images": self.images, "grid": self.grid}

    def __setstate__(self, state: dict) -> None:
        StitchedPsf.__init__(self, state["images"], state["grid"])

    def _write(self) -> bytes:
        return pickle.dumps((self._images, self._grid))

    @staticmethod
    def _read(pkl: bytes) -> StitchedPsf:
        return StitchedPsf(*pickle.loads(pkl))

    def writeFits(self, name: str) -> None:
        """Persist the PSF as a FITS file."""
        raise NotImplementedError("FITS persistence not implemented for StitchedPsf.")

    def readFits(self, name: str) -> None:
        """Persist the PSF as a FITS file."""
        raise NotImplementedError("FITS persistence not implemented for StitchedPsf.")
