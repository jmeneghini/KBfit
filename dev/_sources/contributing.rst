Contributing to KBfit
====================

We welcome contributions to KBfit! This guide will help you get started with contributing to the project.

Development Setup
-----------------

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:

   .. code-block:: bash

      git clone https://github.com/yourusername/KBfit.git
      cd KBfit

3. **Set up the development environment**:

   .. code-block:: bash

      conda activate spectrum
      cmake -B build -DCMAKE_BUILD_TYPE=Debug
      cmake --build build

Code Style
-----------

KBfit follows specific coding conventions:

* **Style**: Enforced by ``.clang-format`` file in the root directory
* **Indentation**: 2 spaces, no tabs
* **Line Length**: 80 characters maximum
* **Naming**: Clear variable names, consistent with existing code
* **Comments**: Doxygen-style documentation for classes and methods

### Automatic Formatting

The project uses clang-format for code formatting:

.. code-block:: bash

   # Format all source files
   find source -name "*.h" -o -name "*.cc" | xargs clang-format -i
   
   # Check formatting
   clang-format --dry-run --Werror source/path/to/file.cc

The CI system will automatically check code formatting and can auto-format pull requests.

Documentation Style
-------------------

* **Header Comments**: Include file description, purpose, and author information
* **Class Documentation**: Use Doxygen format with ``@brief``, ``@param``, ``@return``
* **Method Documentation**: Document performance-critical methods thoroughly
* **Inline Comments**: Explain complex algorithms and optimization decisions

Example:

.. code-block:: cpp

   /**
    * @brief Performance-critical method for evaluating residuals
    * @param fitparams Current fit parameter values
    * 
    * This method is called thousands of times during fitting.
    * Any optimization here has significant impact on performance.
    */
   void evalResidualsAndInvCovCholesky(const std::vector<double>& fitparams);

Testing
-------

* **Unit Tests**: Use doctest framework in ``source_testing/``
* **Test Data**: Include test data in ``source_testing/test_quant_conds/``
* **Integration Tests**: Test with realistic input examples
* **Performance Tests**: Benchmark critical sections

Run tests:

.. code-block:: bash

   ctest
   # or
   ./build/KBfit_tests

Writing Tests
^^^^^^^^^^^^^

.. code-block:: cpp

   #include "doctest.h"
   #include "your_class.h"

   TEST_CASE("Test case description") {
       // Test implementation
       CHECK(result == expected);
   }

Pull Request Process
--------------------

1. **Create a feature branch**:

   .. code-block:: bash

      git checkout -b feature/your-feature-name

2. **Make your changes** following the code style guidelines

3. **Add tests** for new functionality

4. **Update documentation** if needed

5. **Commit your changes**:

   .. code-block:: bash

      git commit -m "Brief description of changes
      
      Detailed explanation of what was changed and why.
      
      Co-Authored-By: Your Name <your.email@example.com>"

6. **Push to your fork**:

   .. code-block:: bash

      git push origin feature/your-feature-name

7. **Create a pull request** on GitHub

### Pull Request Guidelines

* **Title**: Clear, descriptive title
* **Description**: Explain what changes were made and why
* **Tests**: Include relevant test cases
* **Documentation**: Update docs if needed
* **Breaking Changes**: Clearly mark any breaking changes

Areas for Contribution
----------------------

### High Priority

* **Performance Optimization**: Improve computational efficiency
* **Memory Management**: Reduce memory usage in large fits
* **MPI Improvements**: Better parallel processing
* **Error Handling**: More robust error messages

### Medium Priority

* **Documentation**: Improve user guides and examples
* **Testing**: Expand test coverage
* **Build System**: CMake improvements
* **Code Quality**: Refactoring and cleanup

### Low Priority

* **Utilities**: Helper scripts and tools
* **Examples**: More comprehensive examples
* **Visualization**: Tools for result analysis

Performance Considerations
--------------------------

When contributing performance improvements:

* **Profile first**: Use profiling tools to identify bottlenecks
* **Benchmark changes**: Measure performance impact
* **Memory efficiency**: Consider memory access patterns
* **Algorithmic improvements**: Focus on algorithmic optimizations
* **Maintain correctness**: Ensure numerical accuracy is preserved

Common Contribution Types
-------------------------

### Bug Fixes

* **Reproduce the issue**: Create a minimal test case
* **Fix the root cause**: Don't just treat symptoms
* **Add regression tests**: Prevent future occurrences
* **Document the fix**: Explain the issue and solution

### New Features

* **Discuss design**: Open an issue to discuss the feature
* **Follow existing patterns**: Use consistent design patterns
* **Add comprehensive tests**: Test edge cases and error conditions
* **Update documentation**: Include user-facing documentation

### Documentation

* **User guides**: Help users understand features
* **API documentation**: Document all public interfaces
* **Physics background**: Explain theoretical foundations
* **Examples**: Provide working examples

### Optimization

* **Measure impact**: Quantify performance improvements
* **Preserve behavior**: Maintain identical numerical results
* **Document changes**: Explain optimization techniques
* **Consider trade-offs**: Balance performance vs. readability

Code Review Process
-------------------

All contributions go through code review:

1. **Automated checks**: CI runs formatting, tests, and builds
2. **Peer review**: Other developers review the code
3. **Feedback incorporation**: Address review comments
4. **Final approval**: Maintainer approves the changes

### Review Criteria

* **Correctness**: Does the code work as intended?
* **Performance**: Are there performance implications?
* **Maintainability**: Is the code easy to understand and modify?
* **Testing**: Are there adequate tests?
* **Documentation**: Is the code properly documented?

Getting Help
------------

* **GitHub Issues**: Report bugs and request features
* **Discussions**: Ask questions and get help
* **Documentation**: Check existing documentation first
* **Code Examples**: Look at existing code for patterns

Contact
-------

* **GitHub**: Open issues and pull requests
* **Email**: Contact maintainers directly for sensitive issues
* **Community**: Join development discussions

Thank you for contributing to KBfit!