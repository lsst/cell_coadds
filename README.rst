###########
cell_coadds
###########


|build|
|coverage|

.. |build| image:: https://github.com/lsst-dm/cell_coadds/actions/workflows/test.yaml/badge.svg?branch=main
   :target: https://github.com/lsst-dm/cell_coadds/actions/workflows/test.yaml

.. |coverage| image:: https://codecov.io/github/lsst-dm/cell_coadds/branch/main/graph/badge.svg
   :target: https://codecov.io/github/lsst-dm/cell_coadds


``cell_coadds`` is a package being developed for eventual inclusion in the `LSST Science Pipelines <https://pipelines.lsst.io>`_.

It will eventually include code for building coadds of astronomical images in small (few arcsecond) cells, in which only input images that fully contain a cell are included.
This helps mitigate problems with PSF discontinuities that are present in traditional coadds.

At present this package is a rapidly-evolving prototype, and it is being developed alongside the `LSST Dark Energy Science Collaboration's <https://lsstdesc.org/>`_ `descwl_coadd <https://github.com/LSSTDESC/descwl_coadd/>`_ package.

Development workflow
====================

This package currently relies exclusively on GitHub Actions for CI, which test it against the `stackvana <https://anaconda.org/conda-forge/stackvana>`_ conda-forge distribution of the latest weekly release of the LSST Science Pipelines.

Development uses a number of tools that are not in broad use in the LSST codebase (and enforces their use in CI); in some respects this package can be considered a pilot program for using them more broadly, but in the meantime the CI checks that utilize these tools require some special care for developers unfamiliar with them.
They include:

- Python formatting via `black <https://pypi.org/project/black/>`_ and `isort <https://pypi.org/project/isort/>`_.
    LSST developers should install these into their conda environment (using the ``conda-forge`` channel, as usual) and run

    .. code-block:: sh

        $ isort python/ tests/
        $ black python/ tests/

    prior to each commit.

- Python static type-checking via `MyPy <http://mypy-lang.org/>`_.
    MyPy is also available from ``conda-forge``.
    Run it here with:

    .. code-block:: sh

        $ mypy python/ tests/

The best way to use all of these tools is via editor integrations, which should be possible for all major editors.
All necessary configuration files for these tools are included in the repository (``pyproject.toml``, ``mypy.ini`` etc.), and these configurations should be allowed to take precedence over any others to ensure the CI checks that use those configurations are satisfied.
Additionally, you may `install a pre-commit <https://pre-commit.com/#installation>`_ hook to ensure that the staged changes are in accordance with these conventions.

.. code-block:: sh

    $ pip install pre-commit
    $ pre-commit install

Because the CI happens solely on Github Actions and because it can take a long time, auto-merge is enabled for this repository.
As an extra safeguard, an approval is required before merging unlike other LSST repositories.
since it is common for reviewers to approve PRs after suggesting minor changes, **the auto-merge option must be used only after the reviewer's comments are addressed**.
