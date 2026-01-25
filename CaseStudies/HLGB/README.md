### Summary

This case study demonstrates how to model the effect of grain boundary character on hydrogen permeation using **PHIMATS**. The grain boundaries are represented using the double-obstacle potential from the order parameters of the phase-field method, allowing hydrogen interactions to be driven by $\nabla g(\phi)_\mathrm{HAGB}$ and $\nabla g(\phi)_\mathrm{LAGB}$. 

📄 For more background on the model, see:
**"The Effect of grain Boundary misorientation on hydrogen flux using a phase‐field‐based diffusion and trapping model"**
*Advanced Engineering Materials*, **26(22)** (2024) 2401561
[DOI](https://doi.org/10.1002/adem.202401561) | [arXiv PDF](https://arxiv.org/pdf/2412.19129)

While the simulations in that paper were conducted using **finite difference method**, this case study uses the **finite element implementation** in PHIMATS, based on the formulation developed in:

**"A full-field model for hydrogen diffusion and trapping in two-phase microstructures: Application to thermal desorption spectroscopy of duplex stainless steel"**
*Acta Materialia* (2025) 121042
[Open Access](https://doi.org/10.1016/j.actamat.2025.121042)

---

### Instructions

#### 1️⃣ Pre-processing

Run the **`Preprocessing`** notebook:

* Loads the mesh file.
  📌 *Note:* Due to its large size, the mesh file is hosted on **Zenodo**. The notebook provides a direct download link and instructions.
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

```sh
./HLGB
```

---

#### 3️⃣ Post-processing

Run the **`PostProcessing`** notebook to:

* Plot **averaged** quantities, such as hydrogen concentration profiles and exit flux.
* Visualize **full-field** results, including the hydrogen flux vector field.
  
---

PHIMATS is under continuous development and newer versions might function differently. The results shown in this study were generated using **PHIMATS v1.3.0**. To activate this specific version, run:

```bash
git checkout v1.3.0
```

And recompile.