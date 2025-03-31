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
    "CoaddApCorrMapStacker",
    "EMPTY_AP_CORR_MAP",
)


from collections.abc import Iterable

import numpy as np
from frozendict import frozendict

from lsst.afw.image import ApCorrMap
from lsst.afw.math import BoundedField
from lsst.geom import Point2D

from .typing_helpers import SingleCellCoaddApCorrMap

EMPTY_AP_CORR_MAP: SingleCellCoaddApCorrMap = frozendict()
"""Default empty aperture correction map for a single cell coadd."""


class CoaddApCorrMapStacker:
    """Online aperture correction map for a cell-based coadd.

    This class is responsible for implementing the logic to coadd the
    aperture correction values and their uncertainties.

    Parameters
    ----------
    evaluation_point : `~lsst.geom.Point2D`
        The point at which the input aperture correction is evaluated.
    do_coadd_inverse_ap_corr : `bool`, optional
        If True, the inverse aperture correction is applied to the coadd.

    Notes
    -----
    At least one of class variables are set dynamically the first time the
    ``add`` method is called on any instance of this class. This behavior
    is based on the practical assumption that all ``ApCorrMap`` instances will
    have the same set of field names during the entire processing. A schema
    is therefore not expected at the time of initialization.
    """

    _ap_corr_names: Iterable[str] = ()
    """An iterable of algorithm names that have aperture correction values."""
    # This is set when the first time the add method is called on any instance.

    def __init__(self, evaluation_point: Point2D, do_coadd_inverse_ap_corr: bool = True) -> None:
        # Initialize frozen attributes.
        self._evaluation_point = evaluation_point
        self._do_coadd_inverse_ap_corr = do_coadd_inverse_ap_corr

        # Initialize mutable attributes.
        self._total_weight = 0.0
        self._intermediate_ap_corr_map: dict[str, float] = {}

    @classmethod
    def _setup_ap_corr_names(cls, ap_corr_map: ApCorrMap) -> None:
        """Set up the aperture correction name set.

        Parameters
        ----------
        ap_corr_map : `~lsst.meas.base.ApCorrMap`
            The aperture correction map to add.

        Raises
        ------
        RuntimeError
            Raised if the keys in `ap_corr_map` do not end in "_instFlux" or
            "_instFluxErr".
        """
        ap_corr_name_set = set()
        for field_name in ap_corr_map:
            algorithm_name, suffix = field_name.split("_instFlux")
            if suffix not in ("", "Err"):
                raise RuntimeError(f"Invalid field name {field_name} in aperture correction map.")

            ap_corr_name_set.add(algorithm_name)

        ap_corr_name_set.discard("ext_gaap_GaapFlux_1_15x_0_5_instFlux")
        ap_corr_name_set.discard("ext_gaap_GaapFlux_1_15x_0_5_instFluxErr")

        cls._ap_corr_names = tuple(sorted(ap_corr_name_set))

    @property
    def evaluation_point(self) -> Point2D:
        """The point at which the aperture correction is evaluated."""
        return self._evaluation_point

    @property
    def do_coadd_inverse_ap_corr(self) -> bool:
        """If True, the inverse aperture correction is applied to the coadd."""
        return self._do_coadd_inverse_ap_corr

    @property
    def ap_corr_names(self) -> Iterable[str]:
        """Iterable of algorithm names that have aperture correction values."""
        return self._ap_corr_names

    @property
    def total_weight(self) -> float:
        """The total weight of the aperture correction map."""
        return self._total_weight

    def add(self, ap_corr_map: ApCorrMap, weight: float) -> None:
        """Add an aperture correction map to the coadd.

        Parameters
        ----------
        ap_corr_map : `~lsst.meas.base.ApCorrMap`
            The aperture correction map to add.
        weight : `float`
            The weight to apply to the aperture correction map.

        Raises
        ------
        RuntimeError
            Raised if the keys in `ap_corr_map` do not end in "_instFlux" or
            "_instFluxErr".
        ValueError
            Raised if the aperture correction value or its error is missing.
        """
        if not self.ap_corr_names:
            # Lazily initialize the aperture correction name set.
            self._setup_ap_corr_names(ap_corr_map)

        if not self._intermediate_ap_corr_map:
            self._intermediate_ap_corr_map = dict.fromkeys(
                [f"{algorithm_name}_instFlux" for algorithm_name in self.ap_corr_names]
                + [f"{algorithm_name}_instFluxErr" for algorithm_name in self.ap_corr_names],
                0.0,
            )

        # Accumulate the aperture correction values in a temporary dict.
        # This is so that if we error out in the middle, we don't leave the
        # aperture correction map in an inconsistent state.
        temp_ap_corr_map = dict.fromkeys(self._intermediate_ap_corr_map, 0.0)

        for algorithm_name in self.ap_corr_names:
            # Accumulate the aperture correction values.
            ap_corr_field: BoundedField | None
            if (ap_corr_field := ap_corr_map.get(f"{algorithm_name}_instFlux")) is None:
                raise ValueError(f"Missing {algorithm_name} aperture correction map.")

            ap_corr_value = ap_corr_field.evaluate(self.evaluation_point)

            # Calculate the term to accumulate depending on the boolean.
            if self.do_coadd_inverse_ap_corr:
                if ap_corr_value == 0:
                    raise ValueError("This should not have happened. ap_corr_value is zero.")
                else:
                    term = weight / ap_corr_value
            else:
                term = weight * ap_corr_value

            temp_ap_corr_map[f"{algorithm_name}_instFlux"] = term

            # Accumulate the aperture correction error values.
            ap_corr_err_field: BoundedField | None
            if (ap_corr_err_field := ap_corr_map.get(f"{algorithm_name}_instFluxErr")) is None:
                raise ValueError(f"Missing {algorithm_name} aperture correction error map.")

            ap_corr_err_value = ap_corr_err_field.evaluate(self.evaluation_point)

            # Calculate the term to accumulate depending on the boolean.
            if self.do_coadd_inverse_ap_corr:
                term = (weight * ap_corr_err_value) ** 2 / ap_corr_value**4
            else:
                term = (weight * ap_corr_err_value) ** 2

            temp_ap_corr_map[f"{algorithm_name}_instFluxErr"] = term

        # Update the intermediate aperture correction map.
        for key in self._intermediate_ap_corr_map:
            self._intermediate_ap_corr_map[key] += temp_ap_corr_map[key]

        # Add the weight to the total weight.
        self._total_weight += weight

    @property
    def final_ap_corr_map(self) -> SingleCellCoaddApCorrMap:
        """Final coadded aperture correction map.

        This should be called after all aperture correction maps have been
        added.

        Raises
        ------
        RuntimeError
            Raised if the total weight is zero.
        """
        if self.total_weight == 0 or not self.ap_corr_names:
            raise RuntimeError("Cannot get an empty aperture correction map.")

        final_ap_corr_map = dict.fromkeys(self._intermediate_ap_corr_map, 0.0)

        # The transformation equation is different for the aperture correction
        # values and their uncertainties and it also depends on whether we
        # accumulate the aperture corrections or their inverse.

        if self.do_coadd_inverse_ap_corr:
            for algorithm_name in self.ap_corr_names:
                if (
                    inverse_ap_corr_value := self._intermediate_ap_corr_map[f"{algorithm_name}_instFlux"]
                ) > 0:
                    final_ap_corr_map[f"{algorithm_name}_instFlux"] = (
                        self.total_weight / inverse_ap_corr_value
                    )
                final_ap_corr_map[f"{algorithm_name}_instFluxErr"] = (
                    final_ap_corr_map[f"{algorithm_name}_instFlux"] ** 2
                    * np.sqrt(self._intermediate_ap_corr_map[f"{algorithm_name}_instFluxErr"])
                    / self.total_weight
                )
        else:
            for algorithm_name in self.ap_corr_names:
                final_ap_corr_map[f"{algorithm_name}_instFlux"] = (
                    self._intermediate_ap_corr_map[f"{algorithm_name}_instFlux"] / self.total_weight
                )
                final_ap_corr_map[f"{algorithm_name}_instFluxErr"] = (
                    np.sqrt(self._intermediate_ap_corr_map[f"{algorithm_name}_instFluxErr"])
                    / self.total_weight
                )

        # Return the finalized (immutable) aperture correction map.
        return frozendict(final_ap_corr_map)
