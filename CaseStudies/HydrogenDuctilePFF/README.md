### Summary

This case study demonstrates the chemo-mechanics coupling for hydrogen embrittlement in **PHIMATS**. The fully kinetic hydrogen transport model is coupled with a modular geometric phase-field fracture. A novel crack driving force is proposed that uses a weighted combination of elastic strain energy and plastic work densities. A hypertangent scaling function of triaxiality is used such that plastic dissipation contribute to damage only in tensile state. Several case studies were presented to reproduce experimentally observed damage patters in ductile metals with and without hydrogen. More importantly, modeling hydrogen accumulation at dislocations was essential for modeling surface cracking at the necking region of smooth round bar and ductile tearing to embrittled crack in compact tension specimen.  

📄 For more background on the model, see:
- **"A coupled fully kinetic hydrogen transport and ductile phase-field
fracture framework for modeling hydrogen embrittlement"**
(2026)

---

Selected case studies are provided in this directory. Each contain pre-post processing notebooks, driver code and .geo/.msh files by Gmsh. The analysis is either axisymmetirc or plane-strain formulation. In each `PreProcessing` notebook, the optimized `PETSc` cmd options are stated to overcome convergence difficulties.  

- `NRB` notched round bar (cup-and-cone).
- `DoubleNotch` double notched sample.
- `SRB` smooth round bar.
  - `SRB_Air` in air.
  - `SRB_Boundary` in hydrogen with different pressures.
  - `SRB_Precharge` 30 min precharging in 30 MPa hydrogen.
  - `SRB_Rate` effect of loading rate in 30 MPa hydrogen. **NOTE** `SRB_Precharge` has to be completed before running this simulation.   
- `CT`: Compact tension specimen.
  - `CT_Air`: in air.
  - `CT_Hydrogen` in hydrogen.

---

### Instructions

#### 1️⃣ Pre-processing

Run the **`Preprocessing`** notebook:

* Sets up material parameters and generates the **`.in.hdf5`** input files.
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

**Note** Use `PETSc` cmd options specified in the `PreProcessing` notebook to improve convergence.

---

#### 3️⃣ Post-processing

Run the **`PostProcessing`** notebook to plt the load-displacement curves.

---

The results shown in this study were generated using **PHIMATS v1.3.0**. To activate this specific version, run:

```bash
git checkout v1.3.0
```
And recompile.