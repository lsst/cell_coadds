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

import unittest

import lsst.utils.tests
from lsst.cell_coadds import CellCoaddFitsReader
from lsst.cell_coadds._fits import FILE_FORMAT_VERSION


class FitsTestCase(lsst.utils.tests.TestCase):
    """Class for testing FITS specific aspects.

    Some FITS specific aspects are also implemented in test_coadds.py
    """

    def test_read_compatibility(self):
        """Test that the isCompatibleWith method works."""
        self.assertTrue(CellCoaddFitsReader.isCompatibleWith("0.1"))
        self.assertFalse(CellCoaddFitsReader.isCompatibleWith("12345.67"))

        # Check that the reader is compatible with itself.
        self.assertTrue(CellCoaddFitsReader.isCompatibleWith(FILE_FORMAT_VERSION))


class TestMemory(lsst.utils.tests.MemoryTestCase):
    """Check for resource/memory leaks."""


def setup_module(module):  # noqa: D103
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
