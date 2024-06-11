# To do list for `PhiMATS`

This is a list for planned features and developments to be added to `PhiMATS`. Smaller modifications are added per class/function `@todo`.

## Models

- [ ] Incremental loading for linear elasticity
- [ ] Transport model
- [ ] Diffusion and trapping model
- [ ] Stress gradients and stress-driven diffusion
- [ ] Phase-field fracture
- [ ] Neumann BCs for mechanics
- [ ] Small strain plasticity
- [ ] Large strain and non linear geometry

## FEM

- [ ] Modify the HDF IO format.
- [ ] DM and DMPlex in PETSc
- [ ] Parallel
- [ ] Adaptive meshing
- [x] 2D and 3D nodes
- [x] Multi-material/multi-element sets

## Software

- [ ] Make as library
- [ ] Make `PhiMATS` `namespace`
- [ ] Inheritance for Models
- [ ] Versions for the code
- [ ] Update `Make` for HPC (does `CMake` work for that?)
- [ ] If not `CMake`, automate `Make` for release using `shell` commands 
- [ ] `XDMF` output
- [x] Put `Eigen` types in separate header file
- [x] Hierarchy and OOP redesign (elements/solvers/materials)
- [x] `Doxygen` and proper commenting
- [x] `git` and branches thing
