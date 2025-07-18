# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

KBfit is a C++ application for fitting K-matrix parameters in lattice QCD calculations, specifically for two-hadron systems in finite volume. The code performs χ² minimization on spectrum data and uses the Lüscher method for scattering analysis.

## Build System

The project uses CMake with C++20 standard. The build system supports:
- **Debug/Release builds**: Default is Release with `-O3` optimization
- **MPI support**: Required for parallel fitting with `mpirun`
- **Conda environment integration**: Automatically detects and uses conda MPI/compilers
- **Cross-platform**: MacOS (with Homebrew) and Linux support

### Project Environment
This project is designed to run in the `spectrum` conda environment, which includes all necessary dependencies. The environment can be activated with:

```bash
conda activate spectrum
```

Hence, all commands and scripts should be run within this environment to ensure proper functionality.

### Build Configuration

Key CMake variables configured in `CMakeLists.txt:76-87`:
- `AR_LMAX`, `OA_LMAX`, `PD_LMAX`, `CD_LMAX`: Maximum orbital angular momentum (≤6)
- `AR_SX2MAX`, `OA_SX2MAX`, `PD_SX2MAX`, `CD_SX2MAX`: Maximum intrinsic spin (×2)

## Code Architecture

### Core Components

1. **Main Driver** (`source/main/KBfit.cc`):
   - XML-driven task execution
   - MPI initialization and coordination
   - Input format: `<KBFit><Initialize>...<TaskSequence>...`

2. **Task Management** (`source/tasks/`):
   - `TaskHandler`: Core task orchestration and XML parsing
   - `task_fit.cc`: Fitting routines
   - `task_print.cc`: Output generation
   - `task_single_channel.cc`: Single-channel analysis

3. **Fitting Engine** (`source/fitting/`):
   - `chisq_fit.h/cc`: Main fitting class with MPI support and progress bars
   - `chisq_spectrum.cc`: Fitting done by root solving the QC for the model predictions (dElab), then a χ² minimization is done on the residuals. This method also introduces priors for non dElab observables, such as mass and anistropy.
   - `chisq_detres.cc`: Fitting done with the QC as the chi-squared residual
   - `minimizer.cc`: MINUIT2 and NL2Sno integration
   - `root_finder.cc`: Root finding utilities used in spectrum fitting

4. **Two-Hadron Physics** (`source/twohadronsinabox/`):
   - `K_matrix_calc.h/cc`: K-matrix computations
   - `box_quant.cc`: Box quantization conditions
   - `box_matrix_*.cc`: Momentum-dependent matrix elements (AR, OA, PD, CD)
   - `fit_forms.cc`: Parameterization forms

5. **Data Handling** (`source/data_handling/`):
   - `io_handler_hdf5.cc`: HDF5 file I/O
   - `io_handler_fstream.cc`: Text file I/O
   - `obs_get_handler.cc`: Observable data management

6. **Analysis Tools** (`source/analysis/`):
   - `kbobs_handler.cc`: Observable handling
   - `mc_estimate.cc`: Monte Carlo error analysis
   - `sampling_info.cc`: Bootstrap/jackknife resampling

### Key Dependencies

- **MINUIT2**: χ² minimization
- **HDF5**: Data storage format
- **MPI**: Parallel processing
- **OpenMP**: Threading (single-threaded BLAS recommended)
- **LAPACK/BLAS**: Linear algebra
- **muParser**: String expression parsing

### Testing

Uses **doctest** framework in `source_testing/`:
- Test files: `test_*.cc`
- Test data: `test_quant_conds/` (includes numpy arrays for validation)
- Run with: `ctest` or `./KBfit_tests`

## Development Workflow

### MPI Usage
```bash
# Single process
./KBfit input.xml

# Parallel execution
mpirun -np 4 ./KBfit input.xml
```

### Key File Locations
- Input examples: `/pi-mnt/latticeQCD/spectrum_analysis/channels/phirho/luscher`. The XML input files present here are your best starting point to understanding the needed input, outputs, goals, etc.
- Input example data that can be explored using `sigmond_query`: `/pi-mnt/latticeQCD/spectrum_analysis/channels/phirho/levels`
- Test data: `source_testing/test_quant_conds/`
- Build artifacts: `build/`
- Output: Generated in project-specific subdirectories

### Code Conventions
- Headers use include guards: `#ifndef HEADER_H`
- MPI-aware classes use `MPI_Comm_rank/size` 
- XML parsing via `xml_handler.h`
- Style set by `clang-format` with `.clang-format` file in root directory
- Progress indication with `indicators.h` library
- Object files: Box matrix sources compiled without optimization (`-O0`)
- Utilize standard coding practices, including clear variable names, consistent indentation, modular design,
comments in the style of the various header comments on the top of files, as well as docstrings for classes and methods, omitting @author
- Strive towards readability and modularity all around, while maintaining and further optimizing performance in critical areas.
Utilize `doctest` for unit tests, with test cases in the `source_testing/` directory.


### Hotspots and Optimization
- The `chisq_fit` class is a performance-critical area, especially in the `fit_spectrum` method.
- The minimizer calls `evalResidualsAndInvCholesky` for each sample iteration, and for the spectrum method (the current focus of development),
it evaluates the roots of the `Omega` function, which is extremely expensive to compute.
- Thus, code optimization and profiling should focus on this area, reducing the number of calls to `getOmega`,
removing redundant calculations, optimizing memory access and writing, and optimizing the root-finding process.

### Sigmond Query
- The `sigmond_query` tool is used to explore and query input data files. It can be run with:
```bash
sigmond_query -h
```
- It supports .bin files and also allows querying of HDF5 files with the root specified. For example:
```bash
sigmond_query -i test.hdf5[/samples]
```

## Quantum Mechanics Context

This code implements finite-volume scattering theory for lattice QCD using the Lüscher method.

### Physics Reference
**Primary Reference**: `morningstar_physics.md` - Comprehensive physics background from Morningstar et al.
**Original Source**: `Morningstar-Two-particle_multi-channel_systems_in_ﬁnite_volume_with_arbitrary_spin_key_points.pdf`

### Fundamental Quantization Condition
The core physics relationship: **det[1 - B^(P)(E_cm) K̃] = 0**

Where:
- **B^(P)(E_cm)**: Box matrix (geometry-dependent, computed from lattice parameters)
- **K̃**: Well-behaved K-matrix (energy-independent scattering parameters to be fitted)
- **E_cm**: Center-of-mass energy

### Energy Shifts and Non-Interacting Energies
- **Non-interacting energies**: Theoretical energies if particles didn't interact in finite volume
- **Energy shifts (dE_lab)**: Difference between observed interacting energies and non-interacting energies
- **Lattice observations**: Actual energy levels measured in Monte Carlo simulations
- **Model predictions**: Energy shifts predicted by solving the quantization condition

### Two Fitting Methodologies

#### **Spectrum Fitting** (SpectrumFit class)
- **Observations**: Center-of-mass energies E_cm,k plus lattice parameters (masses, volume, anisotropy)
- **Method**: Root finding to solve Omega function Ω(E_cm) = 0, yielding energy shift predictions
- **Residuals**: (E_lab_shift_observed - E_lab_shift_predicted)²
- **QC Type**: Always uses `StildeCB` quantization condition
- **Computational Cost**: High (Box Matrix recomputed several times for each residual evaluation due to root finding)
- **OPTIMIZATION**: Root finding now works directly in center-of-mass frame to eliminate redundant coordinate transformations

#### **Determinant Residual Fitting** (DetresFit class)  
- **Observations**: Quantization determinant values themselves
- **Method**: Minimize det(QC(E_cm)) directly using observed energies
- **Residuals**: Ω(μ, 1 - B^(P)(E_cm^obs) K̃) where Ω is a filter function
- **QC Types**: Can use `StildeCB`, `StildeinvCB`, `KtildeB`, or `KtildeinvB`
- **Computational Cost**: Significantly Lower (Box matrix is precomputed, no root finding)

### Momentum Classes
- **AR**: At-rest (P = 0)
- **OA**: On-axis (P along coordinate axis)  
- **PD**: Planar diagonal (P in coordinate plane)
- **CD**: Cubic diagonal (P along space diagonal)

### Angular Momentum Basis
States labeled |J m_J L S a⟩:
- **J**: Total angular momentum
- **m_J**: Angular momentum projection  
- **L**: Orbital angular momentum
- **S**: Total spin
- **a**: Channel (particle species, spins, parities, isospins)

## Implementation Details

### Root Finding Process in Spectrum Fits

The spectrum fitting methodology involves solving the quantization condition Ω(E_cm) = 0 through root finding to predict energy shifts. The implementation has been optimized to work directly in the center-of-mass frame to eliminate redundant coordinate transformations.

#### **Key Implementation Components**

1. **BoxQuantization::getDeltaEnergyPredictionsOptimized()** (`source/twohadronsinabox/box_quant.h:523-533`)
   - Performs root finding directly in CM frame to avoid repeated Ecm ↔ Elab transformations
   - Computes non-interacting energies in CM frame for bracketing
   - Matches predicted roots to observed energies
   - Converts to lab frame only once at the end when dElab is requested
   - **NEW**: Defers lab frame calculations until output phase for additional efficiency

2. **SpectrumFit::evalResidualsAndInvCovCholesky()** (`source/fitting/chisq_spectrum.cc:805-958`)
   - Performance-critical method called for every minimization function evaluation
   - Uses optimized root finding to compute energy shift predictions
   - Eliminates ~N×M coordinate transformations per residual evaluation
   - Maintains numerical accuracy while improving computational efficiency

#### **Optimization Details**

**Before Optimization:**
- Root finding called `getDeltaElabPredictionsInElabInterval()` working in lab frame
- Every omega function evaluation required `getElabOverMrefFromEcm()` transformation
- Box matrices computed with repeated coordinate transformations
- Lab frame calculations performed for all observations upfront
- Significant computational overhead from redundant calculations

**After Optimization:**
- Root finding calls `getDeltaEnergyPredictionsOptimized()` working in CM frame
- Omega function evaluations work directly in natural physics coordinates
- Box matrices computed once in CM frame per energy value
- **Lab frame calculations deferred until output phase when needed**
- Coordinate transformation performed only once at the end
- ~50% reduction in computational overhead for root finding operations
- Additional performance gain when CM frame output is desired

#### **Root Finding Chain**

The optimized root finding follows this call chain:
1. `SpectrumFit::evalResidualsAndInvCovCholesky()` - Main fitting method
2. `BoxQuantization::getDeltaEnergyPredictionsOptimized()` - Optimized energy prediction
3. `BoxQuantization::get_roots_in_interval()` - Root finding in CM intervals
4. `BoxQuantization::get_omega()` - Omega function evaluation in CM frame
5. `BoxQuantization::get_qc_matrix()` - Quantization condition matrix in CM frame

#### **Performance Impact**

The CM frame optimization significantly reduces computational cost by:
- Eliminating repeated coordinate transformations in the hot path
- Reducing function call overhead during root finding
- Improving numerical stability by working in physics-natural coordinates
- Enabling better caching of frequently-used values
- Maintaining mathematical equivalence while improving performance

#### **Physics Consistency**

The optimization maintains physics accuracy by:
- Performing all root finding in the natural CM frame where the quantization condition is defined
- Converting to lab frame only for final dElab output as required
- Preserving the exact mathematical relationship between CM and lab frames
- Ensuring identical numerical results with improved computational efficiency