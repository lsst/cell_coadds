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

__all__ = ("ImagePlanes",)

from collections.abc import Sequence
from typing import Self

import lsst.shoefits as shf
import pydantic


class ImagePlanes(pydantic.BaseModel):
    """Struct combining the image-like planes in a cell-based coadd."""

    @classmethod
    def view_of(cls, other: ImagePlanes, bbox: shf.Box) -> Self:
        return cls(
            image=other.image[bbox],
            mask=other.mask[bbox],
            variance=other.variance[bbox],
            mask_fractions=(other.mask_fractions[bbox] if other.mask_fractions is not None else None),
            noise_realizations=tuple([n[bbox] for n in other.noise_realizations]),
        )

    @property
    def bbox(self) -> shf.Box:
        """The bounding box common to all image planes."""
        return self.image.bbox

    image: shf.Image
    """The data image itself."""

    mask: shf.Mask
    """An integer bitmask."""

    variance: shf.Image
    """Per-pixel variances for the image."""

    mask_fractions: shf.Image | None
    """The (weighted) fraction of masked pixels that contribute to each
    pixel.
    """

    noise_realizations: Sequence[shf.Image]
    """A sequence of noise realizations that were coadded with the same
    operations that were appled to the data image.
    """

    # TODO DM-45189: validate that bboxes are consistent.
