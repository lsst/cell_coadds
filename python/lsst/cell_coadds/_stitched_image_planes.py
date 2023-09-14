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

__all__ = ("StitchedImagePlanes",)

from abc import abstractmethod
from collections.abc import Callable, Iterator, Mapping, Sequence, Set
from functools import partial
from typing import TYPE_CHECKING, TypeVar

from lsst.afw.image import ImageF, Mask

from . import typing_helpers
from ._image_planes import ImagePlanes
from ._uniform_grid import UniformGrid

_T = TypeVar("_T", bound=typing_helpers.ImageLike)


if TYPE_CHECKING:
    from .typing_helpers import ImageLike


class StitchedImagePlanes(ImagePlanes):
    """An ImagePlanes intermediate base class that stitches together per-cell
    images.

    Parameters
    ----------
    bbox : `Box2I`
        The region over which contiguous piecewise images are desired.

    Notes
    -----
    This class simply inserts subimages from each cell into the full image,
    doing so when an attribute is first accessed to avoid stitching together
    planes that may never be accessed.
    """

    def __init__(self) -> None:
        self._image: ImageLike | None = None
        self._mask: Mask | None = None
        self._variance: ImageLike | None = None
        self._mask_fractions: Mapping[str, ImageLike] | None = None
        self._noise_realizations: Sequence[ImageLike] | None = None

    @property
    @abstractmethod
    def grid(self) -> UniformGrid:
        """Object that defines the piecewise grid this object stitches
        together.

        This may include cells outside the region covered by these image
        planes.
        """
        raise NotImplementedError()

    @property
    @abstractmethod
    def n_noise_realizations(self) -> int:
        """The number of noise realizations cells are guaranteed to have."""
        raise NotImplementedError()

    @property
    @abstractmethod
    def mask_fraction_names(self) -> Set[str]:
        """The names of all mask planes whose fractions were propagated in any
        cell.

        Cells that do not have a mask fraction for a particular name may be
        assumed to have the fraction for that mask plane uniformly zero.
        """
        raise NotImplementedError()

    @property
    def image(self) -> ImageLike:
        # Docstring inherited.
        if self._image is None:
            self._image = self._make_plane(ImageF(self.bbox), lambda planes: planes.image)
        return self._image

    def uncache_image(self) -> None:
        """Remove any cached `image` plane."""
        self._image = None

    @property
    def mask(self) -> Mask:
        # Docstring inherited.
        if self._mask is None:
            self._mask = self._make_plane(Mask(self.bbox), lambda planes: planes.mask)
        return self._mask

    def uncache_mask(self) -> None:
        """Remove any cached `mask` plane."""
        self._mask = None

    @property
    def variance(self) -> ImageLike:
        # Docstring inherited.
        if self._variance is None:
            self._variance = self._make_plane(ImageF(self.bbox), lambda planes: planes.variance)
        return self._variance

    def uncache_variance(self) -> None:
        """Remove any cached `variance` plane."""
        self._variance = None

    @staticmethod
    def _mask_getter(planes: ImagePlanes, name: str) -> Mask:
        return planes.mask.get(name, None)

    @property
    def mask_fractions(self) -> Mapping[str, ImageLike]:
        # Docstring inherited.
        if self._mask_fractions is None:
            # Could make this lazier with a custom Mapping class (only stitch a
            # mask fraction plane if that plane is requested), but not clear
            # it's worth the effort.
            self._mask_fractions = {
                name: self._make_plane(ImageF(self.bbox), partial(self._mask_getter, name=name))
                for name in self.mask_fraction_names
            }
        return self._mask_fractions

    def uncache_mask_fraction(self) -> None:
        """Remove any cached `mask_fraction` planes."""
        self._mask_fractions = None

    @staticmethod
    def _noise_plane_getter(planes: ImagePlanes, i: int) -> ImageF:
        return planes.noise_realizations[i]

    @property
    def noise_realizations(self) -> Sequence[ImageF]:
        # Docstring inherited.
        if self._noise_realizations is None:
            # Could make this lazier with a custom Sequence class (only stitch
            # a noise plane if that plane is requested), but not clear it's
            # worth the effort.
            self._noise_realizations = tuple(
                self._make_plane(ImageF(self.bbox), partial(self._noise_plane_getter, i=i))
                for i in range(self.n_noise_realizations)
            )
        return self._noise_realizations

    def uncache_noise_realizations(self) -> None:
        """Remove any cached `noise_realization` planes."""
        self._noise_realizations = None

    def _make_plane(self, result: _T, getter: Callable[[ImagePlanes], _T | None]) -> _T:
        """Stitch together a single image plane.

        Parameters
        ----------
        result : ImageLike
            The out `~lsst.afw.image.Image` or `~lsst.afw.image.Mask` instance
            covering the full area, to be assigned to.
        getter : `Callable`
            Callable that obtains the appropriate image-like object to assign
            a subimage from, given an `ImagePlanes` instance from a cell inner
            region.

        Returns
        -------
        result : image-like
            The same result object passed in.
        """
        for cell_planes in self._iter_cell_planes():
            common_bbox = cell_planes.bbox.clippedTo(self.bbox)
            if not common_bbox.isEmpty():
                input_plane = getter(cell_planes)
                if input_plane is None:
                    result[common_bbox] = 0
                else:
                    result[common_bbox] = input_plane[common_bbox]
        return result

    @abstractmethod
    def _iter_cell_planes(self) -> Iterator[ImagePlanes]:
        """Iterate over all cell image planes."""
        raise NotImplementedError()
