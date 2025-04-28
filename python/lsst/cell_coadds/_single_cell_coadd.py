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

__all__ = ("SingleCellCoadd",)

from collections.abc import Iterable, Set
from typing import TYPE_CHECKING

from lsst.afw.image import ImageD, ImageF
from lsst.geom import Box2I

from ._coadd_ap_corr_map import EMPTY_AP_CORR_MAP
from ._common_components import CommonComponents, CommonComponentsProperties
from ._image_planes import ViewImagePlanes
from .typing_helpers import ImageLike, SingleCellCoaddApCorrMap

if TYPE_CHECKING:
    from ._identifiers import CellIdentifiers, ObservationIdentifiers
    from ._image_planes import ImagePlanes, OwnedImagePlanes


class SingleCellCoadd(CommonComponentsProperties):
    """A single coadd cell, built only from input images that completely
    overlap that cell.

    Parameters
    ----------
    outer : `OwnedImagePlanes`
        The actual coadded images.
    psf : `ImageD`
        The coadded PSF image.
    inner_bbox : `Box2I`
        The bounding box of the inner region of this cell; must be disjoint
        with but adjacent to all other cell inner regions.
    inputs : `Iterable` [`ObservationIdentifiers`]
        Identifiers of observations that contributed to this cell.
    common : `CommonComponents`
        Image attributes common to all cells in a patch.
    identifiers : `CellIdentifiers`
        Struct of identifiers for this cell.
    aperture_correction_map : `frozendict` [`str`, `float`], optional
        Mapping of algorithm name to aperture correction value for this cell.

    Notes
    -----
    At present we assume a single PSF image per cell is sufficient to capture
    spatial variability, which seems adequate given the results we have so far
    and the cell sizes we intend to use.
    """

    def __init__(
        self,
        outer: OwnedImagePlanes,
        *,
        psf: ImageD,
        inner_bbox: Box2I,
        inputs: Iterable[ObservationIdentifiers],
        common: CommonComponents,
        identifiers: CellIdentifiers,
        # aperture_correction_map: SingleCellCoaddApCorrMap = EMPTY_AP_CORR_MAP,
    ):
        assert outer.bbox.contains(
            inner_bbox
        ), f"Cell inner bbox {inner_bbox} is not contained by outer bbox {outer.bbox}."
        self._outer = outer
        self._psf = psf
        self._inner_bbox = inner_bbox
        self._inner = ViewImagePlanes(outer, bbox=inner_bbox, make_view=self.make_view)
        self._common = common
        # Remove any duplicate elements in the input, sorted them and make
        # them an immutable sequence.
        # TODO: Remove support for inputs as None when bumping to v1.0 .
        self._inputs = tuple(sorted(set(inputs))) if inputs else ()
        self._identifiers = identifiers
        # self._aperture_correction_map = aperture_correction_map

    @property
    def inner(self) -> ImagePlanes:
        """Image planes within the inner region of this cell that is disjoint
        with all other cell inner regions.
        """
        return self._inner

    @property
    def outer(self) -> ImagePlanes:
        """Image planes within the full outer region of this cell."""
        return self._outer

    @property
    def psf_image(self) -> ImageF:
        """The coadded PSF image."""
        return self._psf

    @property
    # TODO: Remove the option of returning empty tuple in v1.0.
    def inputs(self) -> tuple[ObservationIdentifiers, ...] | tuple[()]:
        """Identifiers for the input images that contributed to this cell,
        sorted by their `visit` attribute first, and then by `detector`.
        """
        return self._inputs

    @property
    def visit_count(self) -> int:
        """Number of visits that contributed to this cell."""
        return len(self.inputs)

    @property
    def identifiers(self) -> CellIdentifiers:
        """Struct of unique identifiers for this cell."""
        # This overrides the patch-level property from
        # CommonComponentsProperties to provide cell-level information.
        return self._identifiers

    @property
    def common(self) -> CommonComponents:
        # Docstring inherited.
        return self._common

    # @property
    # def aperture_correction_map(self) -> SingleCellCoaddApCorrMap:
        """Mapping of algorithm name to aperture correction values.

        Returns
        -------
        aperture_correction_map : `frozendict` [`str`, float]
            Mapping of algorithm name to aperture correction values.
        """
        # return self._aperture_correction_map

    # @property
    # def aperture_corrected_algorithms(self) -> Set[str]:
        """An iterable of algorithm names that have aperture correction values.

        Returns
        -------
        aperture_corrected_algorithms : `tuple` [`str`, ...]
            List of algorithms that have aperture correction values.
        """
        # if self._aperture_correction_map:
            # return self._aperture_correction_map.keys()

        # return set()

    def make_view(self, image: ImageLike, bbox: Box2I | None = None) -> ImageLike:
        """Make a view of an image, optionally within a given bounding box.

        Parameters
        ----------
        image : `ImageLike`
            The image to make a view of.
        bbox : `Box2I`, optional
            The bounding box within which to make the view.

        Returns
        -------
        image_view: `ImageLike`
            The view of the image.
        """
        if bbox is None:
            bbox = self._inner_bbox
        return image[bbox]
