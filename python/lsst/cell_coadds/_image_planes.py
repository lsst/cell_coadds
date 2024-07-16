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
from typing import Annotated, Any, Self

import astropy.units
import lsst.shoefits as shf
import numpy as np
import pydantic


class ImagePlanes(pydantic.BaseModel):
    """Struct combining the image-like planes in a cell-based coadd."""

    @classmethod
    def from_bbox(
        cls,
        bbox: shf.Box,
        mask_schema: shf.MaskSchema,
        include_mask_fractions: bool = False,
        n_noise_realizations: int = 0,
        unit: astropy.io.units.Unit | None = None,
        **kwargs: Any,
    ) -> Self:
        return cls(
            image=shf.Image(0.0, bbox=bbox, unit=unit, dtype=np.float32),
            mask=shf.Mask(0, schema=mask_schema, bbox=bbox),
            variance=shf.Image(0.0, bbox=bbox, unit=unit, dtype=np.float32),
            mask_fractions=(shf.Image(0.0, bbox=bbox, dtype=np.float32) if include_mask_fractions else None),
            noise_realizations=tuple(
                [shf.Image(0.0, bbox=bbox, dtype=np.float32, unit=unit) for _ in range(n_noise_realizations)]
            ),
            **kwargs,
        )

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

    # TODO DM-45189: set FITS compression settings, but tiles set dynamically
    # somehow (need to let shoefits.FitsWriteContext override this).

    image: Annotated[shf.Image, shf.FitsOptions(extname="image")]
    """The data image itself."""

    mask: Annotated[shf.Mask, shf.FitsOptions(extname="mask")]
    """An integer bitmask."""

    variance: Annotated[shf.Image, shf.FitsOptions(extname="variance")]
    """Per-pixel variances for the image."""

    mask_fractions: Annotated[shf.Image, shf.FitsOptions(extname="mask_fractions")] | None
    """The (weighted) fraction of masked pixels that contribute to each
    pixel.
    """

    noise_realizations: Annotated[Sequence[shf.Image], shf.FitsOptions(extname="noise_realization")]
    """A sequence of noise realizations that were coadded with the same
    operations that were appled to the data image.
    """

    # TODO DM-45189: validate that bboxes are consistent.

    def copy_from(self, other: ImagePlanes) -> None:
        assert self.bbox == other.bbox
        assert self.mask.schema == other.mask.schema
        assert len(self.noise_realizations) == len(other.noise_realizations)
        self.image.array = other.image.array
        self.mask.array = other.mask.array
        self.variance.array = other.variance.array
        if self.mask_fractions is not None and other.mask_fractions is not None:
            self.mask_fractions = other.mask_fractions
        for a, b in zip(self.noise_realizations, other.noise_realizations, strict=True):
            a.array = b.array
