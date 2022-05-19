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
    "ExtractSubimage",
    "FixedImagePlaneTag",
    "ImagePlaneTag",
    "ImagePlanes",
    "MaskFractionImagePlaneTag",
    "NoiseRealizationImagePlaneTag",
    "ResetOrigin",
    "parse_image_plane_tag_str",
)

import dataclasses
import enum
import re
from typing import AbstractSet, Any, Callable, ClassVar, Iterable, Iterator, Mapping, TypeVar, Union, final

from lsst.afw.image import ImageF, Mask, MaskedImageF
from lsst.geom import Box2I, Point2I

from . import typing_helpers


_S = TypeVar("_S", bound="ImagePlanes")


@dataclasses.dataclass(eq=False)
class ImagePlanes:

    __slots__ = ("image", "mask", "variance", "mask_fractions", "noise_realizations")

    image: ImageF
    """The data image itself."""

    mask: Mask
    """An integer bitmask."""

    variance: ImageF
    """Per-pixel variances for the image."""

    mask_fractions: dict[str, ImageF] = dataclasses.field(default_factory=dict)
    """A mapping from mask plane name to an image of the weighted fraction of
    input pixels with that mask bit set."""

    noise_realizations: list[ImageF] = dataclasses.field(default_factory=list)
    """A sequence of noise realizations that were coadded with the same
    operations that were appled to the data image.
    """

    @property
    def bbox(self) -> Box2I:
        """The bounding box common to all image planes."""
        return self.image.getBBox()

    @property
    def mask_fraction_names(self) -> AbstractSet[str]:
        """The names of all mask planes whose fractions are tracked."""
        return self.mask_fractions.keys()

    @property
    def n_noise_realizations(self) -> int:
        """The number of noise realizations present."""
        return len(self.noise_realizations)

    def iter_tags(self) -> Iterator[ImagePlaneTag]:
        yield from FixedImagePlaneTag.__members__.values()
        for name in self.mask_fraction_names:
            yield MaskFractionImagePlaneTag(name)
        for index in range(self.n_noise_realizations):
            yield NoiseRealizationImagePlaneTag(index)

    @classmethod
    def from_mapping(
        cls: type[_S], planes: Mapping[ImagePlaneTag, typing_helpers.ImageLike], **kwargs: Any
    ) -> _S:
        result = cls(
            planes[FixedImagePlaneTag.IMAGE],
            planes[FixedImagePlaneTag.MASK],
            planes[FixedImagePlaneTag.VARIANCE],
            **kwargs,
        )
        for tag, image in planes.items():
            # Perfect candidate for a Python 3.10 match/destructuring block.
            if type(tag) is MaskFractionImagePlaneTag:
                result.mask_fractions[tag.name] = image
            elif type(tag) is NoiseRealizationImagePlaneTag:
                if tag.index != result.n_noise_realizations:
                    # It's just a pain to deal with these being out of order,
                    # and there's no reason they should ever be.
                    raise ValueError("Missing or out-of-order noise realizations in mapping.")
                result.noise_realizations.append(image)
        return result

    @classmethod
    def from_callback(
        cls: type[_S],
        callback: Callable[[ImagePlaneTag], typing_helpers.ImageLike],
        mask_fraction_names: Iterable[str] = (),
        n_noise_realizations: int = 0,
        **kwargs: Any,
    ) -> _S:
        result = cls(
            callback(FixedImagePlaneTag.IMAGE),
            callback(FixedImagePlaneTag.MASK),
            callback(FixedImagePlaneTag.VARIANCE),
            **kwargs,
        )
        for name in mask_fraction_names:
            result.mask_fractions[name] = callback(MaskFractionImagePlaneTag(name))
        for index in range(n_noise_realizations):
            result.noise_realizations.append(callback(NoiseRealizationImagePlaneTag(index)))
        return result

    def make_mapping(self) -> dict[ImagePlaneTag, typing_helpers.ImageLike]:
        return {tag: tag.get(self) for tag in self.iter_tags()}

    def map(self, func: Callable[[typing_helpers.ImageLike], typing_helpers.ImageLike]) -> ImagePlanes:
        result = ImagePlanes(func(self.image), func(self.mask), func(self.variance))
        for name, image in self.mask_fractions.items():
            result.mask_fractions[name] = func(image)
        for image in self.noise_realizations:
            result.noise_realizations.append(func(image))
        return result

    def __getitem__(self, bbox: Box2I) -> ImagePlanes:
        return self.map(ExtractSubimage(bbox))

    def asMaskedImage(self) -> MaskedImageF:
        """Return an `lsst.afw.image.MaskedImage` view of the image, mask, and
        variance planes.
        """
        return MaskedImageF(self.image, self.mask, self.variance)


class FixedImagePlaneTag(enum.Enum):
    IMAGE = "image"
    MASK = "mask"
    VARIANCE = "variance"

    def __str__(self) -> str:
        return self.value

    def get(self, planes: ImagePlanes) -> typing_helpers.ImageLike:
        return getattr(planes, self.value)

    @property
    def image_type(self) -> type[typing_helpers.ImageLike]:
        return Mask if self is FixedImagePlaneTag.MASK else ImageF


@final
@dataclasses.dataclass(frozen=True)
class MaskFractionImagePlaneTag:
    __slots__ = ("name",)
    name: str

    _REGEX: ClassVar[re.Pattern] = re.compile(r"mf\:(.+)")
    """A regular expression that can be used to extract (as the first unnamed
    group) the name of a mask fraction image from a `SingleImagePlane.name`
    value (`str`).
    """

    def __str__(self) -> str:
        return f"mf:{self.name}"

    def get(self, planes: ImagePlanes) -> typing_helpers.ImageLike:
        return planes.mask_fractions[self.name]

    @property
    def image_type(self) -> type[typing_helpers.ImageLike]:
        return ImageF


@final
@dataclasses.dataclass(frozen=True)
class NoiseRealizationImagePlaneTag:
    __slots__ = ("index",)
    index: int

    _REGEX: ClassVar[re.Pattern] = re.compile(r"noise\:(\d+)")
    """A regular expression that can be used to extract (as the first unnamed
    group, which must be coerced to `int`) the index of a noise realization
    image from a `SingleImagePlane.name` value (`str`).
    """

    def __str__(self) -> str:
        return f"noise:{self.index}"

    def get(self, planes: ImagePlanes) -> typing_helpers.ImageLike:
        return planes.noise_realizations[self.index]

    @property
    def image_type(self) -> type[typing_helpers.ImageLike]:
        return ImageF


ImagePlaneTag = Union[FixedImagePlaneTag, MaskFractionImagePlaneTag, NoiseRealizationImagePlaneTag]


def parse_image_plane_tag_str(string: str) -> ImagePlaneTag:
    if (m := MaskFractionImagePlaneTag._REGEX.fullmatch(string)) is not None:
        return MaskFractionImagePlaneTag(m.group(1))
    if (m := NoiseRealizationImagePlaneTag._REGEX.fullmatch(string)) is not None:
        return NoiseRealizationImagePlaneTag(int(m.group(1)))
    return FixedImagePlaneTag(string)


class ExtractSubimage:

    __slots__ = ("_bbox",)

    def __init__(self, bbox: Box2I):
        self._bbox = bbox

    def __call__(self, original: typing_helpers.ImageLike) -> typing_helpers.ImageLike:
        return original[self._bbox]


class ResetOrigin:

    __slots__ = ("_xy0",)

    def __init__(self, xy0: Point2I):
        self._xy0 = xy0

    def __call__(self, original: typing_helpers.ImageLike) -> typing_helpers.ImageLike:
        result: typing_helpers.ImageLike = original.clone(False)
        result.setXY0(self._xy0)
        return result


class Clone:

    __slots__ = ("_deep",)

    def __init__(self, deep: bool):
        self._deep = deep

    def __call__(self, original: typing_helpers.ImageLike) -> typing_helpers.ImageLike:
        return original.clone(self._deep)
