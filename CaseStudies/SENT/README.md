### Summary

This case study demonstrates how to run a phase-field fracture simulation in **PHIMATS**. The formulation implemented in **PHIMATS** is the geometric phase-field fracture method introduced by Miehe et al. 

📄 For more background on the model, see:

- Christian Miehe et al. **Phase field modeling of fracture in multi-physics problems. Part I. Balance of crack surface and failure criteria for brittle crack propagation in thermo-elastic solids**. *Computer Methods in Applied Mechanics and Engineering* 294 (2015), pp. 449–485. [DOI](https://doi.org/10.1016/j.cma.2014.11.016)

- Christian Miehe et al. **Phase field modeling of fracture in multi-physics problems. Part II. Coupled brittle-to-ductile failure criteria and crack propagation in thermo-elastic–plastic solids**. *Computer Methods in Applied Mechanics and Engineering* 294 (2015), pp. 486–522. [DOI](https://doi.org/10.1016/j.cma.2014.11.016)

The model is applied to an SENT sample in plane strain condition. The mechanics material model is elastoplastic with power-law hardeing. Three options for the crack driving force are available in the `SENT.cxx` simulation driver code :
- Brittle ```CalcDrivForcB``` with the formula $\tilde{\mathcal{D}}^\mathrm{b} = \frac{\tilde{\psi}^\mathrm{e+}}{w_\mathrm{c}}$.
- Elastoplastic ```CalcDrivForcEP``` with the formula $\tilde{\mathcal{D}}^\mathrm{ep} = \frac{(\tilde{\psi}^\mathrm{e+} + \tilde{w}^\mathrm{p})}{w_\mathrm{c}}$.
- Elastoplastic with a threshold ```CalcDrivForcEP_TH``` with the formula $\langle \tilde{\mathcal{D}}^\mathrm{ep} \rangle = \zeta \left\langle \frac{\tilde{\psi}^\mathrm{e+} + \tilde{w}^\mathrm{p}}{\tilde{w}_\mathrm{c}} - 1 \right\rangle$.

Choose **ONLY** one of them and comment out the other two. 

---

### Instructions

#### 1️⃣ Pre-processing

Run the **`Preprocessing`** notebook:

* Sets up material parameters and generates the **`.in.hdf5`** input file.
* Creates the **`.out.hdf5`** file for writing the results.
* Generates **`.xdmf`** file for visualization in **ParaView**.

---

#### 2️⃣ Compilation & Execution

**Step 1: Compile the driver code**

```sh
make
```

* This compiles the provided **`.cxx` driver** and links it with `libphimats.so`.

**Step 2: Run the simulation**

**Note** Use `PETSc` command line options for line-search to improve convergence.

```sh
./SENT -snes_linesearch_type bt -snes_linesearch_damping 0.8 -snes_linesearch_max_it 50 -snes_linesearch_monitor -snes_max_it 100
```

---

#### 3️⃣ Post-processing

Run the **`PostProcessing`** notebook to obtain the total force and crack mouth opening displacement from the _out.hdf5 file and plot them.

---

The results shown in this study were generated using **PHIMATS v1.3.0**. To activate this specific version, run:

```bash
git checkout v1.3.0-beta.1
```
And recompile.