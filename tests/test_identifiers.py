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

from lsst.cell_coadds import CellIdentifiers, ObservationIdentifiers, PatchIdentifiers
from lsst.cell_coadds.test_utils import generate_data_id
from lsst.skymap import Index2D


class IdentifiersTestCase(unittest.TestCase):
    """Tests for UniformGrid and GridIndex/Index2D's C++/Python
    translation.
    """

    @classmethod
    def setUpClass(cls) -> None:  # noqa: D102
        cls.data_id = generate_data_id()  # type: ignore

    def test_cell_identifiers(self) -> None:
        """Test we can construct a CellIdentifiers from a DataCoordinate."""
        cellIdentifier = CellIdentifiers.from_data_id(
            self.data_id, cell=Index2D(x=4, y=2)  # type: ignore [attr-defined]
        )
        self.assertEqual(cellIdentifier.cell, Index2D(x=4, y=2))
        self.assertEqual(cellIdentifier.tract, 9813)

        # Test that we cannot create a CellIdentifiers with the signature of a
        # PatchIdentifiers factory method by accident.
        with self.assertRaises(TypeError):
            CellIdentifiers.from_data_id(self.data_id)  # type: ignore

    def test_observation_identifiers(self) -> None:
        """Test we can construct an ObservationIdentifiers from a
        DataCoordinate.
        """
        observationIdentifier = ObservationIdentifiers.from_data_id(
            self.data_id  # type: ignore [attr-defined]
        )
        self.assertEqual(observationIdentifier.visit, 1234)
        self.assertEqual(observationIdentifier.detector, 9)

    def test_observation_identifiers_with_backup_detector(self) -> None:
        """Test that the optional detector keyword argument does not override
        the value present in the data_id.
        """
        observationIdentifier = ObservationIdentifiers.from_data_id(
            self.data_id,
            backup_detector=42,
        )
        self.assertEqual(observationIdentifier.detector, 9)

    def test_patch_identifiers(self) -> None:
        """Test we can construct a PatchIdentifiers from a DataCoordinate."""
        patchIdentifier = PatchIdentifiers.from_data_id(self.data_id)  # type: ignore [attr-defined]
        self.assertEqual(patchIdentifier.tract, 9813)
        self.assertEqual(patchIdentifier.patch, Index2D(x=4, y=2))

        # Test that we cannot create a PatchIdentifiers with the signature of a
        # CellIdentifiers factory method by accident.
        with self.assertRaises(TypeError):
            PatchIdentifiers.from_data_id(self.data_id, cell=Index2D(x=4, y=2))  # type: ignore


if __name__ == "__main__":
    unittest.main()
