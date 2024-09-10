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

import lsst.cell_coadds.test_utils as test_utils
import lsst.utils.tests
from lsst.cell_coadds import CoaddUnits, CommonComponents, PatchIdentifiers


class CommonComponentsTestCase(lsst.utils.tests.TestCase):
    """Test the construction and interfaces of CommonComponents."""

    def setUp(self) -> None:  # noqa: D102
        self.units = CoaddUnits.nJy
        quantum_data_id = test_utils.generate_data_id()
        self.wcs = test_utils.generate_wcs()
        self.band = quantum_data_id.get("band")
        self.identifiers = PatchIdentifiers.from_data_id(quantum_data_id)
        self.common = CommonComponents(
            units=self.units, wcs=self.wcs, band=self.band, identifiers=self.identifiers
        )

    def test_commonComponentsProperties(self):
        """Test the common properties equal those in the common attribute."""
        self.assertEqual(self.common.units, self.units)
        self.assertEqual(self.common.wcs, self.wcs)
        self.assertEqual(self.common.band, self.band)
        self.assertEqual(self.common.identifiers, self.identifiers)


class TestMemory(lsst.utils.tests.MemoryTestCase):
    """Check for resource/memory leaks."""


def setup_module(module):  # noqa: D103
    lsst.utils.tests.init()


if __name__ == "__main__":
    lsst.utils.tests.init()
    unittest.main()
