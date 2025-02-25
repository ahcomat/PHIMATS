### Summary
This case study demonstrates how to model hydrogen charging and thermal desorption spectroscopy (TDS) in duplex stainless steels using **PHIMATS**.  

<<<<<<< HEAD
For theoretical details, please refer to:  
**A. Hussein et al.**, *A full-field model for hydrogen diffusion and trapping in two-phase microstructures: Application to thermal desorption spectroscopy of duplex stainless steel*. **Acta Mater.**  

=======
>>>>>>> dev-PlasInter
To make these simulations feasible for a standard PC, this setup reduces the domain size to one-third of the original RVE width.  
- **PC version**: Requires ~3.2 GB of RAM.
- **Full RVE version** (as in the paper): **Required ~15 GB RAM** and was run on the **Flemish Supercomputer Center (VSC)**. 

###  Simulation Overview
There are **two main simulations** in this case study:  
1️⃣ **`Charging`** → Simulates **0.2-day hydrogen charging**  
2️⃣ **`TDS`** → Simulates **thermal desorption** at **900 K/h**  

⚠️ **Important:** `Charging` must be run **first**, as it generates the **initial state** required for `TDS`.

---

### Instructions
####  1. Pre-processing
Each directory contains a **`PreProcessing` notebook** that:
- Sets up **input parameters** and generates the **`_in.hdf5` input file**.
- Outputs the **`_out.hdf5` simulation results**.
- Generates **`.xdmf` visualization files** (compatible with **ParaView v5.9.1**).

####  2. Compilation & Execution
1️⃣ **Compile the driver code**  
   - Open a terminal in the case study directory.  
   - Run:  
     ```sh
     make
     ```
   - This compiles the **`.cxx` driver code** and links with `libphimats.a`.  

2️⃣ **Run the simulations**  
   - Execute:  
     ```sh
     ./Charging
     ```
     (Wait for it to finish before running `TDS`.)  
   - Then execute:  
     ```sh
     ./TDS
     ```
     (This step **takes longer** due to stiffness matrix updates.)  

####  3. Expected Run Times
💻 **Charging**: ~20 minutes  
💻 **TDS**: ~2.5 hours (requires stiffness matrix updates at each time step)  

####  4. Post-processing
- The **`PostProcessing` notebook** includes Python scripts for:
  - **Plotting homogenized data**.
  - **Analyzing results** from the output files.