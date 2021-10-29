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

# MyPy really does not like the continueClass pattern.
# type: ignore

from __future__ import annotations

__all__ = ()

from lsst.utils import continueClass
from ._cell_coadds import SimpleGrid


# Different flake8 versions complain about the same problem on different lines.
@continueClass  # noqa: F811
class SimpleGrid:  # noqa: F811
    def __getitem__(self, slices):
        try:
            y_slice, x_slice = slices
            y_start, y_stop, y_step = y_slice.indices(self.shape[0])
            x_start, x_stop, x_step = y_slice.indices(self.shape[1])
        except (TypeError, AttributeError):
            raise TypeError("Grid index must be a 2-d slice.") from None
        return self.subset((y_start, x_start), (y_stop - 1), (x_stop - 1))
