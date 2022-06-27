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

- C++ formatting via `clang-format <https://clang.llvm.org/docs/ClangFormat.html>`_.
    This can also be installed into the usual LSST conda environment (from ``conda-forge``), and then
    run with

    .. code-block:: sh

        $ clang-format include/lsst/cell_coadds/*
        $ clang-format src/*

- Python static type-checking via `MyPy <http://mypy-lang.org/>`_.
    MyPy is also available from ``conda-forge``.
    Run it here with:

    .. code-block:: sh

        $ mypy python/ tests/

    In order to extend type checking to the Python bindings of C++ classes,
    this package includes a typing stub file (``python/lsst/cell_coadds/_cell_coadds.pyi``).
    This should be updated whenever changes are made to the Python bindings.
    Major changes (especially new code) can be mostly handled by the
    ``stubgen`` tool included with ``mypy``, but generally require some hand-editing afterwards (the existing type stubs should certainly not be blindly replaced by ``stubgen`` output).

The best way to use all of these tools is via editor integrations, which should be possible for all major editors.
All necessary configuration files for these tools are included in the repository (``pyproject.toml``, ``mypy.ini``, and ``.clang-format``), and these configurations should be allowed to take precedence over any others to ensure the CI checks that use those configurations are satisfied.
