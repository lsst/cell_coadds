lsst-cell_coadds v30.0.0 (2026-01-18)
=========================================

New Features
------------

- Persist information about the input visits (+detectors) to each cell in separate HDUs. (`DM-40563 <https://rubinobs.atlassian.net/browse/DM-40563>`_)
- Included information about `day_obs` and `physical_filter` in `ObservationIdentifiers` so those metadata are available when reading a file. (`DM-43516 <https://rubinobs.atlassian.net/browse/DM-43516>`_)
- Started including a semantic version number for the file format. (`DM-43586 <https://rubinobs.atlassian.net/browse/DM-43586>`_)
- Add a new APCORR HDU carrying the aperture correction values for each cell. (`DM-48683 <https://rubinobs.atlassian.net/browse/DM-48683>`_)
- Cell-based coadds are no longer assumed to be edge-free. Coadds with edges are fully supported. (`DM-52724 <https://rubinobs.atlassian.net/browse/DM-52724>`_)
- When converting a StitchedCoadd to an afw Exposure, it contains a CoaddInputs.
  However, the records do not contain the individual bounding boxes and WCSs.
  Instead, the validPolygon is already in the coadd coordinates. (`DM-53215 <https://rubinobs.atlassian.net/browse/DM-53215>`_)


API Changes
-----------

- Replace "BAND" in the header with a more standard "FILTER" keyword. (`DM-40563 <https://rubinobs.atlassian.net/browse/DM-40563>`_)
- Allowed `inputs` argument when constructing a `SingleCellCoadd` instance to be any iterable, not just a frozenset. (`DM-43516 <https://rubinobs.atlassian.net/browse/DM-43516>`_)
- Removed `packed` attribute from `ObservationIdentifiers`. (`DM-44233 <https://rubinobs.atlassian.net/browse/DM-44233>`_)
- Add TUNIT4 to the FITS header to denote the units of the variance plane.
  Add BUNIT to the FITS header when converting to an afw Exposure. (`DM-48683 <https://rubinobs.atlassian.net/browse/DM-48683>`_)
- StitchedCoadd objects have a new `set_cell_edges` method to mark the cell boundaries with a mask bit. (`DM-52724 <https://rubinobs.atlassian.net/browse/DM-52724>`_)


Bug Fixes
---------

- Replace TUNI1 with TUNIT2 to denote the units of the image array. (`DM-48683 <https://rubinobs.atlassian.net/browse/DM-48683>`_)
- Save the mask bit definitions in a new MASKDEF binary table. (`DM-52724 <https://rubinobs.atlassian.net/browse/DM-52724>`_)
