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
    "SingleImagePlane",
    "ViewImagePlanes",
)

from abc import ABC, abstractmethod
from dataclasses import dataclass
from functools import partial
from typing import (
    AbstractSet,
    Any,
    Callable,
    Dict,
    Generic,
    Iterator,
    List,
    Mapping,
    Optional,
    Sequence,
    Type,
    overload,
)

from lsst.afw.image import ImageF, Mask, MaskedImageF
from lsst.geom import Box2I

from . import typing_helpers


class ImagePlanes(ABC):
    """Struct interface for the image-like planes we coadd.

    Notes
    -----
    The extends the set of planes in `lsst.afw.image.MaskedImage` by adding
    noise realizations and "mask fraction" images.

    The `ImagePlanes` is designed to be used both by code that knows the
    specific planes it wants (via various properties and other direct
    accessors) and code that works generically over all image planes
    (via iteration, which yields `SingleImagePlane` instances).

    In addition to the methods and properties decorated with `abstractmethod`,
    derived classes should override the various ``uncache_`` methods if they
    lazily compute and cache image planes on-the-fly, and override ``__iter__``
    (delegating to `super`) if they add any planes not included in the base
    class interface.  Derived classes that add new attributes to each plane
    should subclass `SingleImagePlane`, override `_make_single_plane`, and
    override `__iter__` (even if just to update the type annotation).
    """

    MASK_FRACTION_NAME_TEMPLATE = "mf.{}"
    """The `str.format` template used to generate `SingleImagePlane.name`
    values for mask fraction images (`str`).
    """

    MASK_FRACTION_NAME_REGEX = r"mf\.(.+)"
    """A regular expression that can be used to extract (as the first unnamed
    group) the name of a mask fraction image from a `SingleImagePlane.name`
    value (`str`).
    """

    NOISE_REALIZATION_NAME_TEMPLATE = "noise.{03d}"
    """The `str.format` template used to generate `SingleImagePlane.name`
    values for noise realization images (`str`).
    """

    NOISE_REALIZATION_NAME_REGEX = r"noise\.(\d+)"
    """A regular expression that can be used to extract (as the first unnamed
    group, which must be coerced to `int`) the index of a noise realization
    image from a `SingleImagePlane.name` value (`str`).
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

    def uncache_image(self) -> None:
        """Remove any cached `image` plane."""
        pass

    @property
    @abstractmethod
    def mask(self) -> Mask:
        """An integer bitmask."""
        raise NotImplementedError()

    def uncache_mask(self) -> None:
        """Remove any cached `mask` plane."""
        pass

    @property
    @abstractmethod
    def variance(self) -> ImageF:
        """Per-pixel variances for the image."""
        raise NotImplementedError()

    def uncache_variance(self) -> None:
        """Remove any cached `variance` plane."""
        pass

    @property
    @abstractmethod
    def mask_fraction_names(self) -> AbstractSet[str]:
        """The names of all mask planes whose fractions are tracked."""
        raise NotImplementedError()

    @abstractmethod
    def get_mask_fraction(self, name: str) -> ImageF:
        """Return the mask fraction with the given name.

        Parameters
        ----------
        name : `str`
            Name of the mask fraction plane.

        Returns
        -------
        mf : `ImageF`
            Mask fraction image.
        """
        raise NotImplementedError()

    @property
    def mask_fractions(self) -> Mapping[str, ImageF]:
        """A mapping from mask plane name to an image of the weighted fraction
        of input pixels with that mask bit set."""
        return MaskFractionMappingView(self)

    def uncache_mask_fraction(self, name: str) -> None:
        """Remove the given cached mask fraction plane.

        Parameters
        ----------
        name : `str`
            Name of the mask fraction plane to uncache.
        """
        pass

    @property
    @abstractmethod
    def n_noise_realizations(self) -> int:
        """The number of noise realizations present."""
        raise NotImplementedError()

    @abstractmethod
    def get_noise_realization(self, index: int) -> ImageF:
        """Return the noise realization with the given index.

        Parameters
        ----------
        index : `int`
            Index of the noise realization.

        Returns
        -------
        noise : `ImageF`
            Noise image.
        """
        raise NotImplementedError()

    @property
    def noise_realizations(self) -> Sequence[ImageF]:
        """A sequence of noise realizations that were coadded with the same
        operations that were appled to the data image.
        """
        return NoiseRealizationsSequence(self)

    def uncache_noise_realization(self, index: int) -> None:
        """Remove the given cached noise realization plane.

        Parameters
        ----------
        index : `int`
            Index of the noise realization to uncache.
        """
        pass

    def asMaskedImage(self) -> MaskedImageF:
        """An `lsst.afw.image.MaskedImage` view of the image, mask, and
        variance planes.
        """
        return MaskedImageF(self.image, self.mask, self.variance)

    def _make_single_plane(
        self,
        name: str,
        description: str,
        image_type: Type[typing_helpers.ImageLikeType],
        get: Callable[[], typing_helpers.ImageLikeType],
        uncache: Callable[[], None],
    ) -> SingleImagePlane[typing_helpers.ImageLikeType]:
        return SingleImagePlane(
            name=name,
            description=description,
            image_type=image_type,
            get=get,
            uncache=uncache,
        )

    def __iter__(self) -> Iterator[SingleImagePlane]:
        yield self._make_single_plane(
            name="image",
            description="main image",
            image_type=ImageF,
            get=self.get_image,
            uncache=self.uncache_image,
        )
        yield self._make_single_plane(
            name="mask",
            description="integer bitmask",
            image_type=Mask,
            get=self.get_mask,
            uncache=self.uncache_mask,
        )
        yield self._make_single_plane(
            name="variance",
            description="per-pixel variance image",
            image_type=ImageF,
            get=self.get_variance,
            uncache=self.uncache_variance,
        )
        for name in self.mask_fraction_names:
            yield self._make_single_plane(
                name=self.MASK_FRACTION_NAME_TEMPLATE.format(name),
                description=f"fraction of input pixels with mask bit {name} set",
                image_type=ImageF,
                get=partial(self.get_mask_fraction, name),
                uncache=partial(self.uncache_mask_fraction, name),
            )
        for index in range(self.n_noise_realizations):
            yield self._make_single_plane(
                name=self.NOISE_REALIZATION_NAME_TEMPLATE.format(index),
                description=f"noise realization {index}",
                image_type=ImageF,
                get=partial(self.get_noise_realization, index),
                uncache=partial(self.uncache_noise_realization, index),
            )

    def get_image(self) -> ImageF:
        """Return the main image.

        Returns
        -------
        image : `ImageF`
            The same as the `image` property.

        Notes
        -----
        This method duplicates the functionality of the `image` property (and
        always delegates to it) for code written in a more functional style.
        Derived classes should always override the `image` property, not this
        method.
        """
        return self.image

    def get_mask(self) -> Mask:
        """Return the integer bitmask.

        Returns
        -------
        mask : `Mask`
            The same as the `mask` property.

        Notes
        -----
        This method duplicates the functionality of the `mask` property (and
        always delegates to it) for code written in a more functional style.
        Derived classes should always override the `mask` property, not this
        method.
        """
        return self.mask

    def get_variance(self) -> ImageF:
        """Return the per-pixel variance image.

        Returns
        -------
        variance : `ImageF`
            The same as the `variance` property.

        Notes
        -----
        This method duplicates the functionality of the `variance` property
        (and always delegates to it) for code written in a more functional
        style.  Derived classes should always override the `variance` property,
        not this method.
        """
        return self.variance

    @staticmethod
    def mask_fraction_getter(name: str) -> Callable[[ImagePlanes], ImageF]:
        """Return a callable that accesses a single mask fraction image.

        Parameters
        ----------
        name : `str`
            Name of the mask fraction image to return.

        Returns
        -------
        getter : `Callable`
            A callable that takes a single `ImagePlanes` argument and returns
            a single mask fraction image.

        Notes
        -----
        This method should never be overridden by derived classes.
        """

        def getter(self: ImagePlanes) -> ImageF:
            return self.get_mask_fraction(name)

        return getter

    @staticmethod
    def noise_realization_getter(index: int) -> Callable[[ImagePlanes], ImageF]:
        """Return a callable that accesses a noise realization image.

        Parameters
        ----------
        index : `index`
            Index of the noise realization image to return.

        Returns
        -------
        getter : `Callable`
            A callable that takes a single `ImagePlanes` argument and returns
            a noise realization image.

        Notes
        -----
        This method should never be overridden by derived classes.
        """

        def getter(self: ImagePlanes) -> ImageF:
            return self.get_noise_realization(index)

        return getter


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
        mask_fractions: Optional[Mapping[str, ImageF]] = None,
        noise_realizations: Sequence[ImageF] = (),
    ):
        if mask_fractions is None:
            mask_fractions = {}
        self._image = image
        self._mask = mask
        self._variance = variance
        self._mask_fractions = dict(mask_fractions)
        self._noise_realizations = list(noise_realizations)

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
    def mask_fraction_names(self) -> AbstractSet[str]:
        # Docstring inherited.
        return self._mask_fractions.keys()

    def get_mask_fraction(self, name: str) -> ImageF:
        # Docstring inherited.
        return self._mask_fractions[name]

    @property
    def mask_fractions(self) -> Dict[str, ImageF]:
        # Docstring inherited.
        return self._mask_fractions

    @property
    def n_noise_realizations(self) -> int:
        # Docstring inherited.
        return len(self._noise_realizations)

    def get_noise_realization(self, index: int) -> ImageF:
        # Docstring inherited.
        return self._noise_realizations[index]

    @property
    def noise_realizations(self) -> List[ImageF]:
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
        self,
        target: ImagePlanes,
        make_view: Callable[[typing_helpers.ImageLike], typing_helpers.ImageLike],
        bbox: Optional[Box2I] = None,
    ):
        self._target = target
        self._bbox = bbox if bbox is not None else self._target.bbox
        self._make_view = make_view

    @property
    def bbox(self) -> Box2I:
        # Docstring inherited.
        return self._bbox

    @property
    def image(self) -> ImageF:
        # Docstring inherited.
        return self._make_view(self._target.image)

    @property
    def mask(self) -> Mask:
        # Docstring inherited.
        return self._make_view(self._target.mask)

    @property
    def variance(self) -> ImageF:
        # Docstring inherited.
        return self._make_view(self._target.variance)

    def mask_fraction_names(self) -> AbstractSet[str]:
        # Docstring inherited.
        return self._target.mask_fraction_names

    def get_mask_fraction(self, name: str) -> ImageF:
        # Docstring inherited.
        return self._make_view(self._target.get_mask_fraction(name))

    @property
    def n_noise_realizations(self) -> int:
        # Docstring inherited.
        return self._target.n_noise_realizations

    def get_noise_realization(self, index: int) -> ImageF:
        # Docstring inherited.
        return self._make_view(self._target.get_noise_realization(index))


class MaskFractionMappingView(Mapping[str, ImageF]):
    def __init__(self, parent: ImagePlanes):
        self._parent = parent

    def __getitem__(self, key: str) -> ImageF:
        return self._parent.get_mask_fraction(key)

    def __iter__(self) -> Iterator[str]:
        return iter(self._parent.mask_fraction_names)

    def __len__(self) -> int:
        return len(self._parent.mask_fraction_names)


class NoiseRealizationsSequence(Sequence[ImageF]):
    def __init__(self, parent: ImagePlanes):
        self._parent = parent

    @overload
    def __getitem__(self, index: int) -> ImageF:
        pass

    @overload
    def __getitem__(self, s: slice) -> Sequence[ImageF]:
        pass

    def __getitem__(self, i: Any) -> Any:
        if isinstance(i, int):
            return self._parent.get_noise_realization(i)
        elif isinstance(i, slice):
            start, stop, step = i.indices(len(self))
            return tuple(self._parent.get_noise_realization(i) for i in range(start, stop, step))
        else:
            raise TypeError(f"Unexpected index for sequence: {i!r}.")

    def __len__(self) -> int:
        return self._parent.n_noise_realizations


@dataclass(frozen=True, eq=False)
class SingleImagePlane(Generic[typing_helpers.ImageLikeType]):
    name: str
    description: str
    image_type: Type[typing_helpers.ImageLikeType]
    get: Callable[[], typing_helpers.ImageLikeType]
    uncache: Callable[[], None]
