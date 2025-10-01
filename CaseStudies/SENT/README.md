### Summary

This case study demonstrates how to run a phase-field fracture simulation in **PHIMATS**. The formulation implemented in **PHIMATS** is the geometric phase-field fracture method introduced by Miehe et al. 

📄 For more background on the model, see:

Christian Miehe et al. **Phase field modeling of fracture in multi-physics problems. Part I. Balance of crack surface and failure criteria for brittle crack propagation in thermo-elastic solids**. *Computer Methods in Applied Mechanics and Engineering* 294 (2015), pp. 449–485. [DOI](https://doi.org/10.1016/j.cma.2014.11.016)

Christian Miehe et al. **Phase field modeling of fracture in multi-physics problems. Part II. Coupled brittle-to-ductile failure criteria and crack propagation in thermo-elastic–plastic solids**. *Computer Methods in Applied Mechanics and Engineering* 294 (2015), pp. 486–522. [DOI](https://doi.org/10.1016/j.cma.2014.11.016)

The model is applied to an SENT sample in plane strain condition. The mechanics material model is elastoplastic with power-law hardeing. Three options for the crack driving force are available in the `SENT.cxx` simulation driver code :
- Brittle ```CalcDrivForcB``` with the formula $\tilde{\mathcal{D}}^\mathrm{b} = \frac{\tilde{\psi}^\mathrm{e+}}{w_\mathrm{c}}$.
- Elastoplastic ```CalcDrivForcEP``` with the formula $\tilde{\mathcal{D}}^\mathrm{ep} = \frac{(\tilde{\psi}^\mathrm{e+} + \tilde{w}^\mathrm{p})}{w_\mathrm{c}}$.
- Elastoplastic with a threshold ```CalcDrivForcEP_TH``` with the formula $\langle \tilde{\mathcal{D}}^\mathrm{ep} \rangle = \zeta \left\langle \frac{\tilde{\psi}^\mathrm{e+} + \tilde{w}^\mathrm{p}}{\tilde{w}_\mathrm{c}} - 1 \right\rangle$.

Choose **ONLY** one of them and comment out the other two. 

---

### Instructions

#### 1️⃣ Pre-processing

Run the **`Preprocessing`** notebook:

* Sets up material parameters and generates the **`_in.hdf5`** input file.
* Creates the **`_out.hdf5`** file for writing the results.
* Generates **`.xdmf`** file for visualization in **ParaView (v5.9.1 recommended)**.

---

#### 2️⃣ Compilation & Execution

**Step 1: Compile the driver code**

```sh
make
```

* This compiles the provided **`.cxx` driver** and links it with `libphimats.a`.

**Step 2: Run the simulation**

```sh
./SENT
```

---

#### 3️⃣ Post-processing

Run the **`PostProcessing`** notebook to obtain the total force and crack mouth opening displacement from the _out.hdf5 file and plot them.

---

The results shown in this study were generated using **PHIMATS v1.2.0**. To activate this specific version, run:

```bash
git checkout v1.2.0
```
