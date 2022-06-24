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

from itertools import product
from typing import Any, Dict, List

from lsst.skymap import Index2D
from lsst.utils import continueClass

from ._cell_coadds import GridContainer, GridContainerBuilder

__all__: List[str] = []


@continueClass
class GridContainerBuilder:  # noqa: F811
    def __getstate__(self) -> Dict:
        gc: GridContainer = self.finish()
        x_iter = range(self.offset.x, self.offset.x + self.shape.x)
        y_iter = range(self.offset.y, self.offset.y + self.shape.y)
        yx_iter = product(y_iter, x_iter)
        state = {}
        for (y, x) in yx_iter:
            state[Index2D(x, y)] = gc[Index2D(x, y)]
        return state

    def __setstate__(self, state: Dict) -> None:
        x_iter = range(self.offset.x, self.offset.x + self.shape.x)
        y_iter = range(self.offset.y, self.offset.y + self.shape.y)
        yx_iter = product(y_iter, x_iter)
        for (y, x), value in zip(yx_iter, state.values()):
            self.set(x=x, y=y, value=value)
        self.finish()

    def __reduce__(self) -> tuple:
        state = self.__getstate__()
        return (GridContainerBuilder, (self.shape, self.offset), state)


@continueClass
class GridContainer:  # noqa: F811
    def __getstate__(self) -> Dict:
        """Support pickle with pickle.dump(s)."""
        state: Dict[Any, Any] = {"shape": self.shape, "offset": self.offset}
        for y in range(self.offset.y, self.offset.y + self.shape.y):
            for x in range(self.offset.x, self.offset.x + self.shape.x):
                state[Index2D(x, y)] = self[Index2D(x, y)]
        return state

    def __setstate__(self, state: Dict) -> GridContainer:
        """Support unpickle with pickle.load(s)."""
        builder: GridContainerBuilder = GridContainerBuilder(self.shape, self.offset)
        x_iter = range(self.offset.x, self.offset.x + self.shape.x)
        y_iter = range(self.offset.y, self.offset.y + self.shape.y)
        yx_iter = product(y_iter, x_iter)
        for (y, x), value in zip(yx_iter, state):
            builder.set(x=x, y=y, value=value)
        return builder.finish()

    def __reduce__(self) -> tuple:
        builder: GridContainerBuilder = GridContainerBuilder(self.shape, self.offset)
        x_iter = range(self.offset.x, self.offset.x + self.shape.x)
        y_iter = range(self.offset.y, self.offset.y + self.shape.y)
        yx_iter = product(y_iter, x_iter)
        for (y, x), value in zip(yx_iter, self):
            builder.set(x=x, y=y, value=value)
        builder.finish()
        return (builder.finish, (), list(self.__iter__()))
