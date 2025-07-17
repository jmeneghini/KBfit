Changelog
=========

This document tracks changes and improvements to KBfit.

Version 1.0.0 (Development)
----------------------------

**Major Features**

* **Spectrum Fitting Optimization**: Root finding now works directly in center-of-mass frame
  
  - Eliminates redundant Ecm ↔ Elab coordinate transformations
  - ~50% reduction in computational overhead for root finding operations
  - Maintains numerical accuracy while improving performance
  - Uses ``BoxQuantization::getDeltaEnergyPredictionsOptimized()`` for efficient computation

* **Enhanced Documentation System**: Comprehensive Doxygen and Sphinx integration
  
  - Automatic webpage generation on every push
  - GitHub Pages deployment for main and development branches
  - Enhanced API documentation with class diagrams
  - Physics background documentation with mathematical formulas

* **Automated Code Quality**: Integrated clang-format workflow
  
  - Automatic formatting checks on pull requests
  - Auto-formatting for pull request contributions
  - Consistent code style enforcement across the project

**Performance Improvements**

* **Memory Access Optimization**: Pre-allocated temporary vectors in ``evalResidualsAndInvCovCholesky``
  
  - Reduced memory allocations during fitting
  - Cache-friendly memory access patterns
  - Efficient parameter passing using pointer arithmetic

* **Root Finding Optimization**: Improved ``get_DeltaE_predictions`` implementation
  
  - Reduced function call overhead in hot paths
  - Better caching of frequently-used values
  - Optimized energy matching algorithms

**API Changes**

* **BoxQuantization Enhancements**: 
  
  - Added ``getDeltaEnergyPredictionsOptimized()`` method
  - Enhanced energy prediction with flexible output frame selection
  - Improved non-interacting energy calculations

* **SpectrumFit Improvements**:
  
  - Enhanced ``evalResidualsAndInvCovCholesky()`` with CM frame optimization
  - Better documentation of performance-critical methods
  - Improved error handling and validation

**Infrastructure**

* **CI/CD Pipeline**: Enhanced GitHub Actions workflow
  
  - Automated documentation building and deployment
  - Code formatting verification and auto-correction
  - Multi-branch deployment (main and development)

* **Build System**: CMake configuration improvements
  
  - Better dependency detection
  - Enhanced cross-platform support
  - Improved conda environment integration

**Documentation**

* **User Guides**: Comprehensive installation and getting started guides
* **API Reference**: Complete Doxygen-generated API documentation
* **Physics Background**: Detailed theoretical foundation documentation
* **Contributing Guide**: Clear guidelines for contributors

**Bug Fixes**

* **Coordinate Transformation**: Fixed redundant transformations in spectrum fitting
* **Memory Management**: Resolved memory allocation issues in performance-critical paths
* **Error Handling**: Improved error messages and validation

**Technical Debt**

* **Code Organization**: Better separation of concerns in fitting classes
* **Documentation**: Comprehensive inline documentation for complex algorithms
* **Testing**: Enhanced test coverage for optimization features

Previous Versions
-----------------

**v0.9.x (Pre-release)**

* **MPI Implementation**: Full MPI support for parallel fitting
* **Concurrency Fixes**: Resolved concurrency issues in parallel operations
* **Progress Indicators**: Added progress bars for long-running operations
* **Sample Aggregation**: Fixed MPI sample aggregation functionality

**v0.8.x (Pre-release)**

* **Core Functionality**: Basic spectrum and determinant residual fitting
* **Box Quantization**: Implementation of quantization condition calculations
* **K-matrix Calculations**: Support for K-matrix and inverse K-matrix modes
* **Multiple Momentum Classes**: AR, OA, PD, and CD momentum support

Future Roadmap
--------------

**Planned Features**

* **GPU Acceleration**: CUDA/OpenCL support for box matrix calculations
* **Advanced Algorithms**: Improved root finding algorithms
* **Visualization Tools**: Built-in plotting and analysis tools
* **Extended Physics**: Support for additional particle types and interactions

**Performance Goals**

* **Vectorization**: SIMD optimizations for critical loops
* **Parallelization**: Enhanced MPI and OpenMP integration
* **Memory Efficiency**: Further memory usage optimizations
* **Scalability**: Support for larger lattice sizes and more channels

**Quality Improvements**

* **Testing**: Expanded test suite with more comprehensive coverage
* **Documentation**: Interactive tutorials and examples
* **Error Handling**: More robust error detection and recovery
* **Validation**: Enhanced input validation and physics consistency checks

**Breaking Changes**

This section will document any breaking changes in future versions.

Migration Guide
---------------

**From v0.9.x to v1.0.0**

* **API Changes**: No breaking changes in public API
* **Performance**: Automatic performance improvements, no code changes needed
* **Configuration**: Existing XML configuration files remain compatible
* **Build System**: No changes required for existing build setups

**Deprecation Notices**

Currently no deprecated features.

Contributing
------------

To contribute to KBfit:

1. Check the issue tracker for open issues
2. Follow the contribution guidelines in ``contributing.rst``
3. Submit pull requests with clear descriptions
4. Add tests for new features
5. Update documentation as needed

See the :doc:`contributing` guide for detailed instructions.

Acknowledgments
---------------

* **KBfit Development Team**: Core implementation and optimization
* **Lattice QCD Community**: Physics guidance and validation
* **Open Source Contributors**: Bug reports, feature requests, and improvements
* **GitHub Actions**: Automated CI/CD infrastructure
* **Sphinx/Doxygen**: Documentation generation tools