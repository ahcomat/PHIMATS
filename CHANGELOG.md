# Changelog

## [Unreleased]
### Added
- `Logger` for consistent messaging, integrated with `PETSc`.

### Changed
- `makefile` now uses machine-independent path variables: `PHIMATSINCLUDES`, `EIGEN`, `H5ID`, and `H5LD`.

---

## [v0.2.0] - 2025-01-10
### Added
- `IsoHard` material for isotropic hardening in J2 plasticity.

### Changed
- `MechModel`:
  - Implemented `PETSc` `SNES` non-linear solver using incremental iterative analysis.
  - Optimized `Assemble` by removing redundant loops.
- Added J2 plasticity support to:
  - `Hex8` elements.
  - `Quad4` elements.

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
