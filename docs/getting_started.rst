Getting Started
===============

This guide will help you get started with KBfit for K-matrix parameter fitting in lattice QCD.

Basic Usage
-----------

KBfit is driven by XML configuration files. The basic usage pattern is:

.. code-block:: bash

   # Single process
   ./KBfit input.xml
   
   # Parallel execution
   mpirun -np 4 ./KBfit input.xml

Input File Structure
--------------------

KBfit uses XML-driven task execution with the following basic structure:

.. code-block:: xml

   <KBFit>
     <Initialize>
       <!-- Initialization parameters -->
     </Initialize>
     
     <TaskSequence>
       <!-- Task definitions -->
     </TaskSequence>
   </KBFit>

Example Input Files
-------------------

The best starting point for understanding KBfit inputs is the example files located at:

.. code-block:: bash

   /pi-mnt/latticeQCD/spectrum_analysis/channels/phirho/luscher

These XML files demonstrate:

* Spectrum fitting configuration
* Determinant residual fitting
* Multiple ensemble handling
* K-matrix parameterization

Data Exploration
----------------

Before running fits, explore your data using ``sigmond_query``:

.. code-block:: bash

   # For binary files
   sigmond_query -h
   
   # For HDF5 files
   sigmond_query -i test.hdf5[/samples]

Example data files are available at:

.. code-block:: bash

   /pi-mnt/latticeQCD/spectrum_analysis/channels/phirho/levels

Fitting Methodologies
---------------------

KBfit supports two main fitting approaches:

Spectrum Fitting
^^^^^^^^^^^^^^^^

**Best for**: Direct energy level analysis

.. code-block:: xml

   <SpectrumFit>
     <OmegaMu>0.8</OmegaMu>
     <QuantizationCondition>StildeCB</QuantizationCondition>
     
     <RootFinder>
       <MaxIterations>100</MaxIterations>
       <Tolerance>1e-9</Tolerance>
     </RootFinder>
     
     <KtildeMatrix>
       <!-- K-matrix parameters -->
     </KtildeMatrix>
     
     <MCEnsembleParameters>
       <!-- Ensemble configuration -->
     </MCEnsembleParameters>
     
     <KBBlock>
       <!-- Momentum block definition -->
       <LabFrameEnergyShift>
         <MCObs>...</MCObs>
         <NonInteractingPair>pi(1)pi(0)</NonInteractingPair>
       </LabFrameEnergyShift>
       <CMFrameEnergyMin>2.0</CMFrameEnergyMin>
       <CMFrameEnergyMax>4.0</CMFrameEnergyMax>
       <!-- OR -->
       <CMFrameEnergyAutoRangeMargin>0.1<CMFrameEnergyAutoRangeMargin>
     </KBBlock>
   </SpectrumFit>

Determinant Residual Fitting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Best for**: Lower computational cost

.. code-block:: xml

   <DetresFit>
     <OmegaMu>8.0</OmegaMu>
     <QuantizationCondition>StildeCB</QuantizationCondition>
     
     <KtildeMatrix>
       <!-- K-matrix parameters -->
     </KtildeMatrix>
     
     <!-- Similar structure to SpectrumFit -->
   </DetresFit>

Momentum Classes
----------------

KBfit supports different momentum configurations:

* **AR**: At-rest (P = 0)
* **OA**: On-axis (P along coordinate axis)
* **PD**: Planar diagonal (P in coordinate plane)
* **CD**: Cubic diagonal (P along space diagonal)

Example momentum block:

.. code-block:: xml

   <KBBlock>
     <TotalMomentumRay>oa</TotalMomentumRay>
     <TotalMomentumIntSquared>2</TotalMomentumIntSquared>
     <LGIrrep>T1u</LGIrrep>
     <LmaxValues>4 4</LmaxValues>
     
     <!-- Energy data -->
   </KBBlock>

Common Workflows
----------------

1. **Data Preparation**:
   
   * Prepare energy level data in HDF5 format
   * Use ``sigmond_query`` to verify data structure
   * Identify momentum classes and irreps

2. **Initial Fit Setup**:
   
   * Choose fitting methodology (Spectrum vs Determinant)
   * Configure K-matrix parameterization
   * Set energy bounds and tolerances

3. **Parameter Optimization**:
   
   * Start with broad energy ranges
   * Adjust ``OmegaMu`` and root finder tolerances
   * Optimize for computational efficiency

4. **Result Analysis**:
   
   * Examine fit quality and residuals
   * Check parameter correlations
   * Validate physics consistency

Performance Tips
----------------

* **Memory**: Use pre-allocated vectors for large fits
* **MPI**: Match process count to data structure
* **Root Finding**: Adjust tolerances for balance of speed/accuracy
* **Energy Bounds**: Use appropriate CM frame ranges

Common Pitfalls
---------------

* **Coordinate Frames**: Ensure consistent energy frame usage
* **Non-Interacting Pairs**: Verify correct particle assignments
* **Momentum Classes**: Check little group irrep assignments
* **Tolerance Settings**: Balance numerical precision with performance

Next Steps
----------

* Review the :doc:`physics_background` for theoretical context
* Explore the :doc:`api` for detailed class documentation
* Check example inputs for your specific use case
* Join the development discussion for advanced features

For specific questions or issues, consult the troubleshooting section or contact the development team.