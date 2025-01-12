# `PhiMATS`

The `phase-FIEld MATerials Simulator` (pronounced Phi-MATS) is a finite element package for solid mechanics and phase transformations, primarily based on the phase-field method.

## Applications

- Mechanics (linear elasticity, J2 plasticity and crystal plasticity).
- Transport (mainly for hydrogen diffusion and trapping).
- Full-field mesoscale models (phase-field based RVEs).
- Fracture and fatigue (phase-field based).

## Features

- Pre-processing using `Python` interface.
- Modular design using object-oriented `C++` with `inheritance` and `polymorphism`.
- Data storage in `HDF5` format.
- `XDMF` for visualization in `Paraview`.

## Summary

- `H5IO` Wrapper for HDF5, for file IO.
- `Nodes` Handles node coordinates.
- `FiniteElements` For element development.
- `Model` Assembles the global stiffness matrix and applies boundary conditions (nterfaces with `PETSc`). Handles output.
- `Solver` Wrapper for `PETSc` solvers. 

## Documentation

- Theory.
- Publications in `EXAMPLES` folder.
- More details about the functions by running

```bash
doxygen Doxyfile
```

## Installation

The following third party libraries are required

- [Eigen](http://eigen.tuxfamily.org)
- [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
- [PETSc](https://www.mcs.anl.gov/petsc/)

### For installation

- Instructions for PC
- Instructions for HPC

## Usage

In order to combine simulation flexibility and high performance, `PhiMATS` is compiled as a library to be linked with a `driver code`. Examples for such work flow can be found in the folder `EXAMPLES`, which also acts as a tutorial.

## Element naming convention

The elements has naming convention as abbreviation of the geometrical shape, followed by number of nodes, i.e. `Quad4`. Elements formulated for displacement DOFs has no further specifications. for example

- `Quad4`: 4-node quadrilateral.
- `Tri3`: 3-node triangular.
- `Hex8`: 8-node hexahedron.

For heat/mass transfer it will have the letter `T`.

- `Quad4T`: 4-node quadrilateral.
- `Tri3T`: 3-node triangular.
- `Hex8T`: 8-node hexahedron.

For heat/mass transfer coupled with mechanics it will have the letter `MT`. Note that the DOF numbering convention is to have all displacement DOFs first, then the transfer DOFs after. *(make numerical example?)*

- `Quad4MT`: 4-node quadrilateral.
- `Tri3MT`: 3-node triangular.
- `Hex8MT`: 8-node hexahedron.

For general phase-field (Allen-Cahn type) it will have the letter `P`. *(Shall all scalar elements have the same convention?)*

- `Quad4P`: 4-node quadrilateral.
- `Tri3P`: 3-node triangular.
- `Hex8P`: 8-node hexahedron.

For phase-field fracture, i.e. phase-field coupled with mechanics, it will have the letter `MP`.

- `Quad4MP`: 4-node quadrilateral.
- `Tri3MP`: 3-node triangular.
- `Hex8MP`: 8-node hexahedron.

For coupled deformation-fracture-diffusion, i.e. phase-field coupled with mechanics and mass transfer, it will have the letter `MTP`.

- `Quad4MTP`: 4-node quadrilateral.
- `Tri3MTP`: 3-node triangular.
- `Hex8MTP`: 8-node hexahedron.

## License

PhiMATS is available under the [GNU General Public License v3.0 or later](https://www.gnu.org/licenses/gpl-3.0.html).

## Releases

The releases of this project are are available [here](https://github.com/AbduKT/PhiMATSFEM).
