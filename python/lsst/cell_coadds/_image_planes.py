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
from collections.abc import Callable, Sequence
from typing import TYPE_CHECKING

from lsst.afw.image import MaskedImageF

if TYPE_CHECKING:
    from lsst.afw.image import Mask
    from lsst.geom import Box2I

    from .typing_helpers import ImageLike


class ImagePlanes(ABC):
    """Struct interface for the image-like planes we coadd.

    Notes
    -----
    The extends the set of planes in `lsst.afw.image.MaskedImage` by adding
    noise realizations and "mask fraction" images.
    """

    @property
    @abstractmethod
    def bbox(self) -> Box2I:
        """The bounding box common to all image planes."""
        raise NotImplementedError()

    @property
    @abstractmethod
    def image(self) -> ImageLike:
        """The data image itself."""
        raise NotImplementedError()

    @property
    @abstractmethod
    def mask(self) -> Mask:
        """An integer bitmask."""
        raise NotImplementedError()

    @property
    @abstractmethod
    def variance(self) -> ImageLike:
        """Per-pixel variances for the image."""
        raise NotImplementedError()

    @property
    @abstractmethod
    def mask_fractions(self) -> ImageLike | None:
        """The (weighted) fraction of masked pixels that contribute to each
        pixel.
        """
        raise NotImplementedError()

    @property
    @abstractmethod
    def noise_realizations(self) -> Sequence[ImageLike]:
        """A sequence of noise realizations that were coadded with the same
        operations that were appled to the data image.
        """
        raise NotImplementedError()

    def asMaskedImage(self) -> MaskedImageF:
        """Return an `lsst.afw.image.MaskedImage` view of the image, mask, and
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
        image: ImageLike,
        mask: Mask,
        variance: ImageLike,
        mask_fractions: ImageLike | None = None,
        noise_realizations: Sequence[ImageLike] = (),
    ):
        self._image = image
        self._mask = mask
        self._variance = variance
        self._mask_fractions = mask_fractions
        self._noise_realizations = tuple(noise_realizations)

    @property
    def bbox(self) -> Box2I:
        # Docstring inherited.
        return self._image.getBBox()

    @property
    def image(self) -> ImageLike:
        # Docstring inherited.
        return self._image

    @property
    def mask(self) -> Mask:
        # Docstring inherited.
        return self._mask

    @property
    def variance(self) -> ImageLike:
        # Docstring inherited.
        return self._variance

    @property
    def mask_fractions(self) -> ImageLike | None:
        # Docstring inherited.
        return self._mask_fractions

    @property
    def noise_realizations(self) -> Sequence[ImageLike]:
        # Docstring inherited.
        return self._noise_realizations


class ViewImagePlanes(ImagePlanes):
    """An implementation of the `ImagePlanes` interface that extracts views
    from another target `ImagePlanes` instance.

    Parameters
    ----------
    target : `ImagePlanes`
        Planes to construct views of.
    make_view : `Callable`
        Callable that takes an original image plane and returns a view into it.
    bbox : `Box2I`, optional
        Bounding box of the new image plane.  Defaults to ``target.bbox``.
    """

    def __init__(
        self, target: ImagePlanes, make_view: Callable[[ImageLike], ImageLike], bbox: Box2I | None = None
    ):
        self._target = target
        self._bbox = bbox if bbox is not None else self._target.bbox
        self._make_view = make_view

    @property
    def bbox(self) -> Box2I:
        # Docstring inherited.
        return self._bbox

    @property
    def image(self) -> ImageLike:
        # Docstring inherited.
        return self._make_view(self._target.image)

    @property
    def mask(self) -> Mask:
        # Docstring inherited.
        return self._make_view(self._target.mask)

    @property
    def variance(self) -> ImageLike:
        # Docstring inherited.
        return self._make_view(self._target.variance)

    @property
    def mask_fractions(self) -> ImageLike | None:
        # Docstring inherited.
        if self._target.mask_fractions is not None:
            return self._make_view(self._target.mask_fractions)

        return None

    @property
    def noise_realizations(self) -> Sequence[ImageLike]:
        # Docstring inherited.
        # We could make this even lazier with a custom Sequence class, but it
        # doesn't seem worthwhile.
        return tuple(self._make_view(r) for r in self._target.noise_realizations)
