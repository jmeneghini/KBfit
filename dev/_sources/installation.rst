Installation Guide
==================

KBfit is designed to run in the ``spectrum`` conda environment, which includes all necessary dependencies.

Prerequisites
-------------

* **Operating System**: Linux or macOS
* **Compiler**: C++20 compatible compiler (GCC 10+, Clang 12+)
* **Build System**: CMake 3.16+
* **MPI**: Required for parallel fitting
* **Conda**: For dependency management

Dependencies
------------

Required Dependencies:

* **MINUIT2**: χ² minimization
* **HDF5**: Data storage format
* **MPI**: Parallel processing (OpenMPI or MPICH)
* **LAPACK/BLAS**: Linear algebra
* **muParser**: String expression parsing

Testing Dependencies:

* **doctest**: Unit testing framework

Environment Setup
-----------------

1. **Activate the conda environment:**

   .. code-block:: bash

      conda activate spectrum

2. **Verify MPI installation:**

   .. code-block:: bash

      mpirun --version

Build Instructions
------------------

1. **Clone the repository:**

   .. code-block:: bash

      git clone <repository-url>
      cd KBfit

2. **Create build directory:**

   .. code-block:: bash

      mkdir build
      cd build

3. **Configure with CMake:**

   .. code-block:: bash

      cmake .. -DCMAKE_BUILD_TYPE=Release

   Available build types:
   
   * ``Release``: Optimized build with ``-O3`` (default)
   * ``Debug``: Debug build with symbols

4. **Build the project:**

   .. code-block:: bash

      make -j$(nproc)

   Or using CMake:

   .. code-block:: bash

      cmake --build . --parallel

Configuration Options
---------------------

Key CMake variables can be configured in ``CMakeLists.txt``:

.. code-block:: cmake

   # Maximum orbital angular momentum (≤6)
   set(AR_LMAX 4)
   set(OA_LMAX 4)
   set(PD_LMAX 4)
   set(CD_LMAX 4)
   
   # Maximum intrinsic spin (×2)
   set(AR_SX2MAX 2)
   set(OA_SX2MAX 2)
   set(PD_SX2MAX 2)
   set(CD_SX2MAX 2)

Testing
-------

Run the test suite:

.. code-block:: bash

   ctest

Or run tests directly:

.. code-block:: bash

   ./KBfit_tests

Troubleshooting
---------------

**Common Issues:**

* **MPI not found**: Ensure MPI is installed and ``mpirun`` is in PATH
* **HDF5 errors**: Check that HDF5 is properly installed in conda environment
* **Compilation errors**: Verify C++20 compiler support

**Environment Issues:**

Make sure you're in the correct conda environment:

.. code-block:: bash

   conda info --envs
   conda activate spectrum

**Build Issues:**

Clean build directory if encountering issues:

.. code-block:: bash

   rm -rf build
   mkdir build
   cd build
   cmake ..
   make

Performance Considerations
--------------------------

* **OpenMP**: Single-threaded BLAS is recommended for optimal performance
* **Compiler Flags**: Release build uses ``-O3`` optimization
* **Memory**: Ensure sufficient RAM for large-scale fits
* **MPI**: Use appropriate number of processes for your system

For optimal performance, consider:

.. code-block:: bash

   export OMP_NUM_THREADS=1
   export MKL_NUM_THREADS=1
   mpirun -np 4 ./KBfit input.xml