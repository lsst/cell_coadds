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

import dataclasses
from abc import abstractmethod
from typing import AbstractSet, Callable, Dict, Iterable, Iterator, List, Optional, Type, cast

from lsst.afw.image import ImageF, Mask

from . import typing_helpers
from ._cell_coadds import UniformGrid
from ._image_planes import ImagePlanes, SingleImagePlane


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

    def __init__(self, mask_fraction_names: Iterable[str], n_noise_realizations: int) -> None:
        self._image: Optional[ImageF] = None
        self._mask: Optional[Mask] = None
        self._variance: Optional[ImageF] = None
        self._mask_fractions: Dict[str, Optional[ImageF]] = dict.fromkeys(mask_fraction_names)
        self._noise_realizations: List[Optional[ImageF]] = [None] * n_noise_realizations

    @property
    @abstractmethod
    def grid(self) -> UniformGrid:
        """Object that defines the piecewise grid this object stitches
        together.

        This may include cells outside the region covered by these image
        planes."""
        raise NotImplementedError()

    @property
    def image(self) -> ImageF:
        # Docstring inherited.
        if self._image is None:
            self._image = self._stitch_plane(ImageF(self.bbox), ImagePlanes.get_image)
        return self._image

    def uncache_image(self) -> None:
        # Docstring inherited.
        self._image = None

    @property
    def mask(self) -> Mask:
        # Docstring inherited.
        if self._mask is None:
            self._mask = self._stitch_plane(Mask(self.bbox), ImagePlanes.get_mask)
        return self._mask

    def uncache_mask(self) -> None:
        # Docstring inherited.
        self._mask = None

    @property
    def variance(self) -> ImageF:
        # Docstring inherited.
        if self._variance is None:
            self._variance = self._stitch_plane(ImageF(self.bbox), ImagePlanes.get_variance)
        return self._variance

    def uncache_variance(self) -> None:
        # Docstring inherited.
        self._variance = None

    @property
    def mask_fraction_names(self) -> AbstractSet[str]:
        # Docstring inherited.
        return self._mask_fractions.keys()

    def get_mask_fraction(self, name: str) -> ImageF:
        # Docstring inherited.
        if (result := self._mask_fractions[name]) is None:
            result = self._stitch_plane(ImageF(self.bbox), ImagePlanes.mask_fraction_getter(name))
            self._mask_fractions[name] = result
        return result

    def uncache_mask_fraction(self, name: str) -> None:
        # Docstring inherited.
        self._mask_fractions[name] = None

    @property
    def n_noise_realizations(self) -> int:
        # Docstring inherited.
        return len(self._noise_realizations)

    def get_noise_realization(self, index: int) -> ImageF:
        # Docstring inherited.
        if (result := self._noise_realizations[index]) is None:
            result = self._stitch_plane(ImageF(self.bbox), ImagePlanes.noise_realization_getter(index))
            self._noise_realizations[index] = result
        return result

    def uncache_noise_realization(self, index: int) -> None:
        # Docstring inherited.
        self._noise_realizations[index] = None

    def _stitch_plane(
        self,
        result: typing_helpers.ImageLikeType,
        getter: Callable[[ImagePlanes], Optional[typing_helpers.ImageLikeType]],
    ) -> typing_helpers.ImageLikeType:
        """Stitch together a single image plane.

        Parameters
        ----------
        result : image-like
            The out `lsst.afw.image.Image` or `lsst.afw.image.Mask` instance
            covering the full area, to be assigned to.
        getter : `Callable`
            Callable that obtains the appropriate image-like object to assign
            a subimage from, given an `ImagePlanes` instance from a cell inner
            region.  May return `None` to represent a zero subimage.

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

    def _make_single_plane(
        self,
        name: str,
        description: str,
        image_type: Type[typing_helpers.ImageLikeType],
        get: Callable[[], typing_helpers.ImageLikeType],
        uncache: Callable[[], None],
    ) -> StitchedSingleImagePlane[typing_helpers.ImageLikeType]:
        return StitchedSingleImagePlane(
            name=name,
            description=description,
            image_type=image_type,
            get=get,
            uncache=uncache,
            grid=self.grid,
        )

    def __iter__(self) -> Iterator[StitchedSingleImagePlane]:
        yield from cast(Iterator[StitchedSingleImagePlane], super().__iter__())


@dataclasses.dataclass(frozen=True, eq=False)
class StitchedSingleImagePlane(SingleImagePlane[typing_helpers.ImageLikeType]):
    grid: UniformGrid
