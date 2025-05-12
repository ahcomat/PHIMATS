### Summary

This case study demonstrates how to model **hydrogen-mechanics** interactions at a **blunting crack-tip** using the **fully kinetic** formulation implemented in PHIMATS. In this formulation, the spatial gradient of the normalized dislocation density, $\nabla \bar{\rho}$, is the diffusion driving force for hydrogen accumulation at plastically deformed regions, as opposed to the local-equilibrium-based models. $\bar{\rho}$ is calculated from a Taylor hardening model with the Kocks-Mecking-Estrin dislocation density evolution equation. These simulations employ one-way coupling, therefore a staggered solution procedure is adopted for numerical efficiency. 

---

###  Simulation Overview

There are three simulations folders:
1️⃣ **`Loading`** → Boundary layer model for localized plastic deformation at the crack-tip.
2️⃣ **`Diffusion`** → The effect of the mechanical fields $\bar{\rho}$ and $\sigma_\mathrm{h}$ on hydrogen transport for precharged and insulated boundary conditions.
3️⃣ **`EquilBC`** → The effect of the mechanical fields $\bar{\rho}$ and $\sigma_\mathrm{h}$ on hydrogen transport for uncharged and equilibrium concentration boundary conditions.

---

### Instructions

#### 1️⃣ Pre-processing

Each case study directory contains a **`PreProcessing` notebook** that:

* Sets up **input parameters** and generates the **`_in.hdf5`** input file.
* Defines the simulation outputs and creates the **`_out.hdf5`** results file.
* Exports an **`.xdmf`** file for visualization in **ParaView** (recommended: v5.9.1).

#### 2️⃣ Compilation & Execution

**Step 1: Compile the driver code**

```bash
make
```

* Run this in each case study directory.
* It compiles the corresponding **`.cxx` driver** and links with `libphimats.a`.

**Step 2: Run the simulation**

```bash
./Loading
```

* This step solves the **mechanical loading** and take longer due to slow convergence in regions of localized plastic deformation.
* 🔍 Tip: use `./Loading -snes_monitor` to print the residual norm to the terminal during execution.
* 🕒 **Wait for this step to finish** before running the transport simulations.

Then run:

```bash
./Diffusion     # for precharged diffusion with insulated boundary conditions
./EquilBC       # for simulations with equilibrium concentration boundary conditions
```

#### 3️⃣ Post-processing

Each case includes a **`PostProcessing` notebook** that provides Python scripts to:

* Extract nodal results ahead of the crack tip from the output files.
* Plot and analyze **Hydrostatic stress, dislocation density and hydrogen concentration** fields.

---

