# Two-Particle Multi-Channel Systems in Finite Volume with Arbitrary Spin - Key Points

## Physics Context: Finite Volume Scattering Theory

### Quantization Condition
The fundamental equation relating finite-volume energies to infinite-volume scattering amplitudes:

**det[1 - B^(P)(E_cm) K̃] = 0**

Where:
- **B^(P)(E_cm)**: Box matrix (depends on momentum P and energy E_cm)
- **K̃**: Well-behaved K-matrix (energy-independent scattering parameters)
- **E_cm**: Center-of-mass energy

### Energy Shifts and Non-Interacting Energies
- **Non-interacting energies**: Theoretical energies if particles didn't interact
- **Energy shifts**: Difference between observed interacting energies and non-interacting energies
- **Lattice observations**: Actual energy levels measured in Monte Carlo simulations

## Two Fitting Methods

### 15.1 Spectrum Method
**Used in SpectrumFit class**

**Observations**: Center-of-mass energies E_cm,k plus lattice parameters (masses, volume, anisotropy)

**Procedure**:
1. For each momentum P and irrep Λ, obtain lab-frame two-particle energies E_k
2. Convert to center-of-mass energies E_cm,k  
3. **Root finding**: Solve Omega function Ω(E_cm) = 0 to find energy shift predictions
4. These roots are the **model predictions** for the energy shifts
5. Compute χ² residuals: (E_shift^obs - E_shift^pred)²

**Key challenge**: Model predictions depend on observations (need masses, lattice size for box matrix)
**Solution**: Include masses, lattice parameters as both observations and model parameters

### 15.2 Determinant Residual Method  
**Used in DetresFit class**

**Observations**: Only the quantization determinant values themselves

**Procedure**:
1. Use observed energies and lattice parameters directly in box matrix B^(P)(E_cm^obs)
2. **Residuals**: Ω(μ, 1 - B^(P)(E_cm^obs) K̃) where Ω is a filter function
3. **No root finding**: Minimize the determinant values directly
4. Much simpler computationally than spectrum method

**Filter Function**: Ω(μ, A) = det(A)/det[(μ² + AA†)^(1/2)]
- Bounded between -1 and 1
- Handles large matrix eigenvalues gracefully
- Zero when quantization condition is satisfied

## Key Physical Quantities

### Center-of-Mass Energy
E_cm = √(E² - P²)

### Momentum and Channel Parameters
For each decay channel a:
- q_cm,a²: Center-of-mass momentum squared
- u_a²: Scaled momentum parameter
- s_a: Kinematic factor

### Angular Momentum Basis
States labeled |J m_J L S a⟩:
- J: Total angular momentum
- m_J: Angular momentum projection
- L: Orbital angular momentum  
- S: Total spin
- a: Channel (particle species, spins, parities)

### Momentum Classes
- **AR**: At-rest (P = 0)
- **OA**: On-axis (P along coordinate axis)
- **PD**: Planar diagonal (P in coordinate plane)
- **CD**: Cubic diagonal (P along space diagonal)