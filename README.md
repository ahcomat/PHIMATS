<img src="PHIMATS.png" width="200"> 

**PHIMATS** (Phase-field Multiphysics Materials Simulator) is a finite element library in C++, originally designed for simulating mesoscale **hydrogen-material interactions** using phase-field based representative volume element (RVE). However, its modular design allows it to be applied to a **wide range of multiphysics problems**.

---

## Applications
- **Phase-field based RVE models** for **hydrogen-material interactions**  
- **Solid mechanics**   
- **Heat and mass transfer** modeling  

---

## Features
- **Pre/Post processing via Python interface** for flexible input/output handling  
- **Object-oriented C++** with **inheritance & polymorphism** for modularity  
- **HDF5-based data storage** for efficient data handling  
- **XDMF format support** for visualization in **ParaView**  

---

## Documentation
- **Companion theory manual:** *[Finite Element Theory for PHIMATS](https://arxiv.org/abs/2502.16283)*  
- **Publications & examples:** Available in the `CaseStudies` directory  
- **Doxygen documentation:** [Online API Reference](https://ahcomat.github.io/doxygen-phimats/index.html)  
 
---

### Usage

To balance **simulation flexibility** with **high performance**, `PHIMATS` is compiled as a **library (`libphimats.so`)**, which must be linked with a **custom driver code** in a simulation folder.  

**Example driver codes** and **hands-on tutorials** can be found in the `CaseStudies` directory. The `makefile` in the subdirectories handles the compilation and linking automatically. These examples demonstrate how to use **PHIMATS** for various **hydrogen-material interaction simulations** and related applications.  

---

### Citation  
For citing **PHIMATS**, please use:  

- [![DOI](https://zenodo.org/badge/755482016.svg)](https://doi.org/10.5281/zenodo.18376381)

- **A. Hussein**, *Finite Element Theory for PHIMATS*, 2025. doi: 10.48550/arXiv.2502.16283.

Additionally, consider citing relevant references from the `CaseStudies` directory if applicable.  

---

### License
**PHIMATS** is released under the **[GNU General Public License v3.0](https://www.gnu.org/licenses/gpl-3.0.html)**, or later.  

