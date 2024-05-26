# Element foumulations

The elements has naming convention as abbreviation of the geometrical shape, followed by number of nodes, i.e. `Quad4`. Elements formulated for displacement DOFs has no further specifications.  for ther

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

# Classes

## `BaseMaterial`

Handles material types. Reads data from HDF5_in file. Includes:

**BaseMechanics**: For mechanical behavior. Includes:

1. **`PlaneStrain`**
2. **`PlaneStress`**
3. **`Elastic3D`**

**BaseTransport**: For heat/mass transport. Includes:

1. **`Diffusion`**
2. **`Thermal`**

### TODO
- [ ] Make materials with name and sets in the HDF5_in file.

## `Nodes`

For reading and storing nodal coordinated from HDF5_in file.

## `Elements`

A container to read and store element data. A general constructor

```cpp
Elements(Material mat, Nodes nodes, H5IO H5File)
```
### virtual functions
- Assemble global stiffness matrix.

### Element specific functions
- Reads data.
- Builds local shape, derivative and stiffness matrix.
- Writes the integration point data (like stresses) to HDF5_out file.

## `Model`

Basically, assembles the global stiffness matrices. Handles the interface with the solver, which is in general PETSc. Also handles boundary conditions. 

```cpp
Model(Elements element, Nodes nodes, H5IO H5File)
```