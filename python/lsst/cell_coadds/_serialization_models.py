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

__all__ = ()

from typing import Self

import pydantic

from ._common_components import CommonComponents
from ._exploded_coadd import ExplodedCoadd
from ._identifiers import ObservationIdentifiers
from ._multiple_cell_coadd import MultipleCellCoadd
from ._to_upstream import CellIndex


class CellInputsModel(pydantic.BaseModel):
    index: CellIndex
    observations: list[ObservationIdentifiers] = pydantic.Field(default_factory=list)


class MultipleCellCoaddModel(pydantic.BaseModel):
    metadata: CommonComponents
    data: ExplodedCoadd
    inputs: list[CellInputsModel] = pydantic.Field(default_factory=list)

    @classmethod
    def from_cell_coadd(cls, cell_coadd: MultipleCellCoadd) -> Self:
        result = cls(
            metadata=cell_coadd.common,
            data=ExplodedCoadd.from_cell_coadd(cell_coadd),
        )
        for cell in cell_coadd.cells.values():
            result.inputs.append(CellInputsModel(index=cell.identifiers.cell, observations=list(cell.inputs)))
        return result
