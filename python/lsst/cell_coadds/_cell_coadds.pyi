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
    "GridContainerBuilder",
    "GridContainer",
    "StitchedPsf",
    "UniformGrid",
)

from typing import Callable, Generic, Iterator, Tuple, TypeVar, overload

from lsst.afw.image import ImageD
from lsst.geom import Box2I, Extent2I, Point2I
from lsst.meas.algorithms import ImagePsf
from lsst.skymap import Index2D

_T = TypeVar("_T")
_T_co = TypeVar("_T_co", covariant=True)
_U = TypeVar("_U")

class GridContainerBuilder(Generic[_T]):
    @overload
    def __init__(self, shape: Tuple[int, int]) -> None: ...
    @overload
    def __init__(self, shape: Tuple[int, int], offset: Tuple[int, int]) -> None: ...
    @property
    def offset(self) -> Index2D: ...
    @property
    def shape(self) -> Index2D: ...
    def __len__(self) -> int: ...
    @overload
    def __setitem__(self, arg0: Tuple[int, int], arg1: _T) -> None: ...
    @overload
    def __setitem__(self, arg0: int, arg1: _T) -> None: ...
    def finish(self) -> GridContainer[_T]: ...

class GridContainer(Generic[_T_co]):
    def __init__(self, builder: GridContainerBuilder[_T_co]) -> None: ...
    @property
    def offset(self) -> Index2D: ...
    @property
    def shape(self) -> Index2D: ...
    @property
    def first(self) -> _T_co: ...
    @property
    def last(self) -> _T_co: ...
    @overload
    def __getitem__(self, arg0: Tuple[int, int]) -> _T_co: ...
    @overload
    def __getitem__(self, arg0: int) -> _T_co: ...
    def __iter__(self) -> Iterator[_T_co]: ...
    def __len__(self) -> int: ...
    def rebuild_empty(self) -> GridContainerBuilder[_U]: ...
    def rebuild_transformed(self, callback: Callable[[_T_co], _U]) -> GridContainerBuilder[_U]: ...

class StitchedPsf(ImagePsf):
    def __init__(self, images: GridContainer[ImageD], grid: UniformGrid) -> None: ...

class UniformGrid:
    @overload
    def __init__(self, bbox: Box2I, cell_size: Extent2I) -> None: ...
    @overload
    def __init__(self, bbox: Box2I, shape: Tuple[int, int]) -> None: ...
    @overload
    def __init__(self, cell_size: Extent2I, shape: Tuple[int, int], min: Point2I = Point2I()) -> None: ...
    def bbox_of(self, position: Tuple[int, int]) -> Box2I: ...
    def index(self, position: Point2I) -> Index2D: ...
    @property
    def bbox(self) -> Box2I: ...
    @property
    def cell_size(self) -> Extent2I: ...
    @property
    def shape(self) -> Index2D: ...
