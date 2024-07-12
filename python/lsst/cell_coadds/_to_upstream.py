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
    "Pair",
    "PixelIndex",
    "PixelShape",
    "CellIndex",
    "CellShape",
    "PatchIndex",
)

from collections.abc import Iterable, Sequence
from typing import Self, TypeVar, cast, final, overload

_T = TypeVar("_T")


class Pair(Sequence[_T]):
    def __init__(self, *, y: _T, x: _T):
        self._tuple = (y, x)

    @classmethod
    def from_xy(cls, xy: Iterable[_T]) -> Self:
        x, y = xy
        return cls(x=x, y=y)

    @classmethod
    def from_yx(cls, yx: Iterable[_T]) -> Self:
        y, x = yx
        return cls(y=y, x=x)

    @property
    def y(self) -> _T:
        return self._tuple[0]

    @property
    def x(self) -> _T:
        return self._tuple[1]

    def __str__(self) -> str:
        return f"(y={self.y}, x={self.x})"

    def __repr__(self) -> str:
        return f"{type(self).__name__}{self!s}"

    @overload
    def __getitem__(self, key: int) -> _T: ...

    @overload
    def __getitem__(self, key: slice) -> Sequence[_T]: ...

    def __getitem__(self, key: int | slice) -> _T | Sequence[_T]:
        return self._tuple[key]

    def __len__(self) -> int:
        return 2

    def __eq__(self, other: object) -> bool:
        if type(other) is type(self):
            return self._tuple == cast(Pair, other)._tuple
        return False

    # TODO: pydantic_core serialization hooks


@final
class PixelIndex(Pair[int]):
    pass


@final
class PixelShape(Pair[int]):
    pass


@final
class CellIndex(Pair[int]):
    pass


@final
class CellShape(Pair[int]):
    pass


@final
class PatchIndex(Pair[int]):
    pass
