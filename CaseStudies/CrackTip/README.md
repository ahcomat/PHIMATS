### Summary

In this case study, the fully kinetic model for hydrogen–dislocation interactions is employed. A boundary-layer model (BLM) is used to simulate a blunting crack tip with localized plastic deformation. Several aspects are demonstrated:

- Isotropic hardening based on dislocation density evolution from the Kocks-Mecking-Estrin equation.
- How the hydrogen chemical potential is affected by the interaction of hydrostatic stress and dislocation density.
- The effect of loading rate on hydrogen accumulation.
- Dislocation pipe diffusion.
- The effect of equilibrium on concentration Dirichlet boundary conditions. 

More technical details could be found in the article: A fully kinetic model for hydrogen transport near a blunting crack tip. Int. J. Plast. (2025): 104406. https://doi.org/10.1016/j.ijplas.2025.104406 (Open access).

###  Simulation Overview

Because of the one-way coupling in this case study, i.e. mechanical loading influences hydrogen diffusion, a staggered coupling scheme is employed for numerical efficiency. The mechanical loading is solved first to compute the hydrostatic stress and normalized dislocation density fields, which are then used by the hydrogen transport simulations. Since rate-independent plasticity is used with quasi-static time increments, the physical time increments are determined from the mass transport simulations. There are **three simulations** in this case study:

1️⃣ **`Loading`** → Simulates the mechanical loading of the BLM to $K_I = 89.2 \, \text{MPa}\sqrt{\text{m}}$  in 500 time increments
2️⃣ **`Diffusion`** → Simulates hydrogen redistribution in an insulated model with initial uniform concentration of 3.453e-3 mol/m³
3️⃣ **`EquilBC`** → Simulates hydrogen diffusion in initially uncharged sample with prescribed concentration at the crack face. 

⚠️ **Important:** `Loading` must be run first before the other two simulations. 


The results shown in this study were generated using **PHIMATS v1.1.1**. To activate this specific version, run:

```bash
git checkout v1.1.1
```
