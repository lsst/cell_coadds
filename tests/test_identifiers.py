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
from lsst.daf.butler import DataCoordinate, DimensionUniverse
from lsst.skymap import Index2D


class IdentifiersTestCase(unittest.TestCase):
    """Tests for UniformGrid and GridIndex/Index2D's C++/Python
    translation."""

    @classmethod
    def setUpClass(cls) -> None:
        universe = DimensionUniverse()

        instrument = universe["instrument"]
        instrument_record = instrument.RecordClass(name="test", detector_max=109, visit_max=10000)

        skymap = universe["skymap"]
        skymap_record = skymap.RecordClass(name="test_skymap")

        band = universe["band"]
        band_record = band.RecordClass(name="r")

        visit = universe["visit"]
        visit_record = visit.RecordClass(id=12345, instrument="test")

        detector = universe["detector"]
        detector_record = detector.RecordClass(id=9, instrument="test")

        physical_filter = universe["physical_filter"]
        physical_filter_record = physical_filter.RecordClass(name="r", instrument="test", band="r")

        patch = universe["patch"]
        patch_record = patch.RecordClass(skymap="test_skymap", tract=9813, patch=42, cell_x=4, cell_y=2)

        record = dict(
            instrument=instrument_record,
            visit=visit_record.id,
            detector=detector_record.id,
            patch=patch_record,
            tract=9813,
            band=band_record,
            skymap=skymap_record,
            physical_filter=physical_filter_record,
        )
        data_id = DataCoordinate.standardize(record, universe=universe)  # type: ignore
        cls.data_id = data_id.expanded(record)  # type: ignore

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
        self.assertEqual(observationIdentifier.visit, 12345)
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
