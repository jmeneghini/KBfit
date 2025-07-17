KBfit Documentation
===================

**KBfit** is a C++ application for fitting K-matrix parameters in lattice QCD calculations, 
specifically for two-hadron systems in finite volume using the Lüscher method.

.. toctree::
   :maxdepth: 2
   :caption: User Guide:

   readme
   installation
   getting_started
   physics_background

.. toctree::
   :maxdepth: 2
   :caption: API Reference:

   api
   
.. toctree::
   :maxdepth: 2
   :caption: Development:

   contributing
   changelog

Key Features
------------

* **Spectrum Fitting**: Advanced root finding in center-of-mass frame for energy shift predictions
* **Determinant Residual Fitting**: Direct minimization of quantization determinants
* **MPI Parallelization**: Efficient parallel processing for large-scale fits
* **Multiple Momentum Classes**: Support for AR, OA, PD, and CD momentum configurations
* **Optimized Performance**: Memory-efficient algorithms with reduced coordinate transformations

Physics Background
------------------

KBfit implements finite-volume scattering theory for lattice QCD using the Lüscher method.
The core physics relationship is the quantization condition:

.. math::

   \det[1 - B^{(P)}(E_{cm}) \tilde{K}] = 0

Where:

* :math:`B^{(P)}(E_{cm})`: Box matrix (geometry-dependent)
* :math:`\tilde{K}`: Well-behaved K-matrix (energy-independent scattering parameters)
* :math:`E_{cm}`: Center-of-mass energy

Recent Optimizations
--------------------

The spectrum fitting methodology has been optimized to work directly in the center-of-mass frame,
eliminating redundant coordinate transformations and improving computational efficiency by ~50%.

Quick Start
-----------

.. code-block:: bash

   # Activate conda environment
   conda activate spectrum
   
   # Build the project
   cmake -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build
   
   # Run a fit
   ./build/KBfit input.xml

For detailed installation and usage instructions, see the :doc:`installation` and :doc:`getting_started` guides.

Indices and Tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
