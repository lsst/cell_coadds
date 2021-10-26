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
    "ImagePlanes",
    "OwnedImagePlanes",
    "ViewImagePlanes",
)

from abc import ABC, abstractmethod
from typing import Sequence

from lsst.afw.image import ImageF, Mask, MaskedImageF
from lsst.geom import Box2I


class ImagePlanes(ABC):
    """Struct interface for the image-like planes we coadd.

    Notes
    -----
    The extends the set of planes in `lsst.afw.image.MaskedImage` by adding
    noise realizations and an "interpolation fraction" image.
    """

    @property
    @abstractmethod
    def bbox(self) -> Box2I:
        """The bounding box common to all image planes."""
        raise NotImplementedError()

    @property
    @abstractmethod
    def image(self) -> ImageF:
        """The data image itself."""
        raise NotImplementedError()

    @property
    @abstractmethod
    def mask(self) -> Mask:
        """An integer bitmask."""
        raise NotImplementedError()

    @property
    @abstractmethod
    def variance(self) -> ImageF:
        """Per-pixel variances for the image."""
        raise NotImplementedError()

    @property
    @abstractmethod
    def missing_fraction(self) -> ImageF:
        """An image of the fraction of input pixels that had to be interpolated
        due to missing or bad data."""
        raise NotImplementedError()

    @property
    @abstractmethod
    def noise_realizations(self) -> Sequence[ImageF]:
        """A sequence of noise realizations that were coadded with the same
        operations that were appled to the data image.
        """
        raise NotImplementedError()

    def asMaskedImage(self) -> MaskedImageF:
        """An `lsst.afw.image.MaskedImage` view of the image, mask, and
        variance planes.
        """
        return MaskedImageF(self.image, self.mask, self.variance)


class OwnedImagePlanes(ImagePlanes):
    """An implementation of the `ImagePlanes` interface backed by actual
    afw objects.
    """

    def __init__(
        self,
        *,
        image: ImageF,
        mask: Mask,
        variance: ImageF,
        missing_fraction: ImageF,
        noise_realizations: Sequence[ImageF] = (),
    ):
        self._image = image
        self._mask = mask
        self._variance = variance
        self._missing_fraction = missing_fraction
        self._noise_realizations = tuple(noise_realizations)

    @property
    def bbox(self) -> Box2I:
        # Docstring inherited.
        return self._image.getBBox()

    @property
    def image(self) -> ImageF:
        # Docstring inherited.
        return self._image

    @property
    def mask(self) -> Mask:
        # Docstring inherited.
        return self._mask

    @property
    def variance(self) -> ImageF:
        # Docstring inherited.
        return self._variance

    @property
    def missing_fraction(self) -> ImageF:
        # Docstring inherited.
        return self._missing_fraction

    @property
    def noise_realizations(self) -> Sequence[ImageF]:
        # Docstring inherited.
        return self._noise_realizations


class ViewImagePlanes(ImagePlanes):
    """A lazy implementation of the `ImagePlanes` interface that just takes
    subimages of another target `ImagePlanes` instance.
    """

    def __init__(self, target: ImagePlanes, bbox: Box2I):
        self._target = target
        self._bbox = bbox

    @property
    def bbox(self) -> Box2I:
        # Docstring inherited.
        return self._bbox

    @property
    def image(self) -> ImageF:
        # Docstring inherited.
        return self._target.image[self._bbox]

    @property
    def mask(self) -> Mask:
        # Docstring inherited.
        return self._target.mask[self._bbox]

    @property
    def variance(self) -> ImageF:
        # Docstring inherited.
        return self._target.variance[self._bbox]

    @property
    def missing_fraction(self) -> ImageF:
        # Docstring inherited.
        return self._target.missing_fraction[self._bbox]

    @property
    def noise_realizations(self) -> Sequence[ImageF]:
        # Docstring inherited.
        # We could make this even lazier with a custom Sequence class, but it
        # doesn't seem worthwhile.
        return tuple(r[self._bbox] for r in self._target.noise_realizations)
