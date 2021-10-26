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

__all__ = ("StitchedCellCoadd",)

from typing import TYPE_CHECKING, Callable, Optional, Protocol, Sequence, TypeVar

from lsst.afw.image import ImageF, Mask
from lsst.geom import Box2I

from ._image_planes import ImagePlanes

if TYPE_CHECKING:
    from ._multiple_cell_coadd import MultipleCellCoadd
    from ._single_cell_coadd import SingleCellCoadd


_S = TypeVar("_S")


class _BoxSubset(Protocol):
    """Interface for objects that can be indexed by a `Box2I`, returning a
    view of the same type.
    """

    def __getitem__(self: _S, bbox: Box2I) -> _S:
        pass

    def __setitem__(self: _S, bbox: Box2I, other: _S) -> None:
        pass


_T = TypeVar("_T", bound=_BoxSubset)


class StitchedCellCoadd(ImagePlanes):
    """A lazy-evaluation coadd that stitches together images from adjacent
    cells.

    Parameters
    ----------
    cell_coadd : `MultipleCellCoadd`
        Cell-based coadd to stitch together.
    bbox : `Box2I`, optional
        The region over which a contiguous coadd is desired.  Defaults to
        ``cell_coadd.inner_bbox``.

    Notes
    -----
    This class simply inserts subimages from each cell into the full image,
    doing so when an attribute is first accessed to avoid stitching together
    planes that may never be accessed.
    """

    def __init__(self, cell_coadd: MultipleCellCoadd, *, bbox: Optional[Box2I] = None):
        self._cell_coadd = cell_coadd
        if bbox is None:
            bbox = cell_coadd.inner_bbox
        elif not cell_coadd.inner_bbox.contains(bbox):
            raise ValueError(
                f"Cell coadd inner bounding box {cell_coadd.inner_bbox} does not "
                f"contain stitch target area {bbox}."
            )
        self._bbox = bbox
        self._image: Optional[ImageF] = None
        self._mask: Optional[Mask] = None
        self._variance: Optional[ImageF] = None
        self._missing_fraction: Optional[ImageF] = None
        self._noise_realizations: Optional[Sequence[ImageF]] = None

    @property
    def bbox(self) -> Box2I:
        # Docstring inherited.
        return self._bbox

    @property
    def image(self) -> ImageF:
        # Docstring inherited.
        if self._image is None:
            self._image = self._make_plane(ImageF(self.bbox), lambda planes: planes.image)
        return self._image

    @property
    def mask(self) -> Mask:
        # Docstring inherited.
        if self._mask is None:
            self._mask = self._make_plane(Mask(self.bbox), lambda planes: planes.mask)
        return self._mask

    @property
    def variance(self) -> ImageF:
        # Docstring inherited.
        if self._variance is None:
            self._variance = self._make_plane(ImageF(self.bbox), lambda planes: planes.variance)
        return self._variance

    @property
    def missing_fraction(self) -> ImageF:
        # Docstring inherited.
        if self._missing_fraction is None:
            self._missing_fraction = self._make_plane(
                ImageF(self.bbox), lambda planes: planes.missing_fraction
            )
        return self._missing_fraction

    @property
    def noise_realizations(self) -> Sequence[ImageF]:
        # Docstring inherited.
        if self._noise_realizations is None:
            # Could make this lazier with a custom Sequence class (only stitch
            # a noise plane if that plane is requested), but not clear it's
            # worth the effort.x
            self._noise_realizations = tuple(
                self._make_plane(ImageF(self.bbox), lambda planes: planes.noise_realizations[i])
                for i in range(self._cell_coadd.n_noise_realizations)
            )
        return self._noise_realizations

    def _make_plane(self, result: _T, getter: Callable[[ImagePlanes], _T]) -> _T:
        """Stitch together a single image plane.

        Parameters
        ----------
        result : image-like
            The out `lsst.afw.image.Image` or `lsst.afw.image.Mask` instance
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
        cell: SingleCellCoadd
        for cell in self._cell_coadd.cells.flat:  # type: ignore
            common_bbox = cell.inner.bbox.clippedTo(self.bbox)
            if not common_bbox.isEmpty():
                result[common_bbox] = getter(cell.inner)[common_bbox]
        return result
