# Changelog

## [v1.3.0] - Code Management - PFF/hydrogen chemo-mechanical coupling (Unreleased)

Major architectural refactor of pre-/post-processing and PFF/hydrogen chemo-mechanical coupling.

### Backward Incompatibilities
- Deprecate support for `.inp` mesh on favor of Gmsh `.msh` format.
- Dedicated `MeshManager` class to manage the mesh. 
- Separated mesh and RVE data from the input file:
  - `.mesh.hdf5`: Stores global nodal coordinates and connectivity.
  - `.rve.hdf5`: Dedicated file for phase-field RVE microstructures.
- A physics-based extension naming convention (`mech`, `diff` and `pff`) for all `.in` and `.out` files to support modular staggered multiphysics coupling. 
- `CaseStudies` update: All examples are modified for compatibility with v1.3.0.
  
### Added
- `Dockerfile` providing a pre-configured, high-performance environment with optimized `PETSc`, `Eigen`, and `HDF5`.
- `FEM_utils/MeshManager` a dedicated class to manage the mesh. 
- `INSTALL.md` new installation instructions guide.
- `WSL_QuickStart.md` WSL environment setup guide for Windows users.
- `configure_env.sh` to automate setting up environment variables for **Manual Installation**, eliminating the need for manual `.bashrc` edits.

### Changed
- `PreProcessing`
  - Deprecate geometry-based functions (`TensileDisp2D`, `TensileDisp3D`, `PermeationX`, `WriteConBCs`, `WriteDispBCs`)
  - Implemented `AssignDirichletBC` and `WriteBCVTK` as general-purpose functions leveraging Gmsh's `Physical groups`.

### Build System & Infrastructure
- Shared library `.so` build support compared to the previous static library `.a` for reduced driver code binary size.
- Shifted to a modular `CMake` architecture for the core library, enabling cleaner dependency management and faster compilation.
- Introduced **Docker-based** workflows for consistent, platform-independent simulation environments.

## [v1.2.0] - Phase-field fracture (01-10-2025)

Implementing the geometric phase-field fracture for modeling damage evolution.  

### Backward Incompatibilities
- The `MechModel` constructor now takes an additional optional argument `const string kspType` to specify the KSP solver type used internally (e.g., `KSPPREONLY`, `KSPGMRES`). 

### Changed
- `PreProcessing` 
  - Takes new `SimulType`: `PFF`, to manage phase-field fracture simulations.
  - Plane stress is no longer accepted for plasticity models (not yet supported).
- `Models/MechModel` monitors the convergence of `SNES` solver via `Logger` class and saves it to the simulation log file.
- `Materials/Mechanics/IsoHard` 
  - Power-law hardening can calculate the passive dislocation density evolution based on KME model. As a result, KME parameters are now required in preprocessing. If dislocation evolution is not needed, set the KME parameters to zero.
  - Include new function `ReturnMapping2D_PFF` for return mapping algorithm including damage effects from PFF. 

### Added
- `Materials/PFF/PFF` phase-field fracture material. 
- `FiniteElements/PFF/Quad4PFF` quad4 element to manage phase-field fracture models. 
- `Models/PFFModel` model class to manage phase-field fracture simulations. 

## [v1.1.1] - Hydrogen-GB interactions (11-05-2025)

Implementing hydrogen interactions with high- and low-angle grain boundaries. Note that this is only for `Tri3TH` elements. 

### Changed
- `FiniteElements/Trapping/Tri3TH` implementing HLGB hydrogen interaction.

### Added
- `Materials/Trapping/TrapHLGB` material model for hydrogen interactions with high- and low-angle grain boundaries. 
- `CaseStudies/FluxGB` case study for hydrogen accumulation and flux at grain boundaries. 
- `CaseStudies/HLGB` case study for the effect of GB misorientation on hydrogen flux. 

## [v1.1.0] - Hydrogen-mechanics interactions (01-05-2025)

Including the hydrogen-mechanics interactions, i.e. hydrostatic stresses and dislocations.

### Backward Incompatibilities
- The `SimulationData` preprocessing dict now should include a constant concentration value `conB` for trapping simulations. 

### Added
- `CaseStudies/CrackTip` directory for hydrogen interactions using boundary layer model. 
- `Materials/Trapping/MechTrap` for modeling hydrogen interactions with hydrostatic stresses and dislocations. 

### Changed
- `Materials/Mechanics/IsoHard` added the `Kocks-Mecking-Estrin` hardening model based on dislocation density.
- `FiniteElements/Trapping/Quad4TH` hydrogen-mechanics interactions.
- `Models/TrappingModel` equilibrium boundary concentrations for mechanics interactions. 

## [v1.0.0] - First Release 🚀 (16-02-2025)

This marks the **first stable** release of PHIMATS, featuring case studies and core functionality.

### Added
- `CaseStudies` directory.

## [v0.5.0]
### Added
- `WriteOutputFile` in PreProcessing for centralized creation of _out.hdf5 output files.
- `environment.yml` file for managing `phimats` conda environment dependencies. 

### Changed
- Error handling in `OpenFileHDF5` and `CloseFileHDF5` in PreProcessing.
- `mpi` handles COMM in `Logger` instead of `PETSC`.
- `Logger` now controls all command-line interface messages.

## [v0.4.0]
### Added
- `BoundaryConditions.py` dedicated module for handling boundary conditions.
- `ReadMesh` function for parsing INP files and reading mesh.
- `WriteXDMF` class for centralized handing of XDMF files.

### Changed
- `preprocessing.py` to `PreProcessing.py`.
- `postprocessing.py` to `PostProcessing.py`.
- Restructuring of `PreProcessing`

## [v0.3.0]
### Added
- `Logger` for consistent messaging, integrated with `PETSc`.

### Changed
- `makefile` now uses machine-independent path variables: `PHIMATSINCLUDES`, `EIGEN`, `H5ID`, and `H5LD`.
- Compile PhiMATS as a library `-libphimats`.
- `LinearTransport`
  - `GMRES` option as solver. 
  - Include `logger` for handling user interface.
- `TrappingModel`
  - Updated matrix allocation strategy.
  - Optimized `Assemble` by removing redundant loops.
  - Include `logger` for handling user interface.

---

## [v0.2.0] - 2025-01-10
### Added
- `IsoHard` material for isotropic hardening in J2 plasticity.

### Removed
- `PlaneStrain.cxx` and `PlaneStress.cxx` : Integrated functionality into `LinearElasti.cxx`.

### Changed
- `MechModel`:
  - Implemented `SNES` non-linear solver using incremental iterative analysis.
  - Optimized `Assemble` by removing redundant loops.
- Added J2 plasticity support to:
  - `Hex8` elements.
  - `Quad4` elements.
- `H5IO`:
  - Added `ReadString` for reading `string` HDF5 input file.
  - Templates in `ReadField1D` and `ReadField2D` to support `int` and `double` types. 
  - Added `try-catch` for improved error handling.

### Fixed
- `MechModel`:
  - Updated matrix allocation strategy to use `MAT_NEW_NONZERO_ALLOCATION_ERR` instead of `MAT_NEW_NONZERO_LOCATION_ERR` for improved memory handling.

---

## [v0.1.0] - 2024-09-01
### Added
- `Tri3H` and `Quad4H` elements for hydrogen diffusion and trapping.
- `TrappingModel` for grain-boundaries and second-phase materials.
- Support for Thermal Desorption Spectroscopy (TDS) in second-phase analysis.

---

## [v0.0.1] - 2024-08-01
### Added
- Initial experimental version with:
  - Integration with the `Eigen` library.
  - `HDF5` file output support.
  - Integration with `PETSc` `KSP` solver.
  - Finite elements:
    - `Tri3`, `Quad4`, and `Hex8` elements for displacement DOFs.
    - `Tri3T` and `Quad4T` elements for temperature/concentration DOFs.
  - Materials:
    - Linear elastic materials (3D, plane-strain, plane-stress).
    - Heat/Mass transfer materials.
  - Models for solving finite elements and handling IO:
    - `MechModel` for mechanics problems.
    - `TransportModel` for heat and mass transfer.
