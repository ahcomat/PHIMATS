import numpy as np
import h5py
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Any

@dataclass
class PhysicsConfig:
    """Schema for Simulation Parameters."""
    SimulName: str
    PhysicsType: str      # Mechanical, Transport, PFF
    PhysicsCategory: str  # Elastic, Plastic, MechTrappingPFF, 2PhaseTrapping, etc.
    nSteps: int
    dt: Optional[float] = None
    Temperature: Optional[float] = 293.15   # Default Room Temp
    R: Optional[float] =8.31446261815324    # Default to SI units
    conB: Optional[float] = 0               # Boundary Concentration
    presBCs: List[Any] = field(default_factory=list)
    exitNodes: List[int] = field(default_factory=list)

@dataclass
class MeshConfig:
    """Metadata passed from MeshManager."""
    nTotNodes: int
    nTotElements: int
    nDim: int
    materialNames: List[str]

class PreProcessing:
    def __init__(self, config: PhysicsConfig, mesh: MeshConfig, materials: Dict):
        self.config = config
        self.mesh = mesh
        self.materials = materials
        
        # Configuration 
        self.ext_map = {"Mechanical": "mech", "Transport": "diff", "PFF": "pff"}
        self.physicsRegistry = {
            "Mechanical": ["Elastic", "Plastic"],
            "Transport": ["Transport", "2PhaseTrapping", "GBTrapping", "HLGBTrapping", "MechTrapping", "MechTrappingPFF"],
            "PFF": ["PFF", "ChemoMech"]
        }
        self._allowed_Analysis = ["3D", "PlaneStrain", "PlaneStress", "PlaneStrainPFF"]
        self._allowed_Isotropy = ["Isotropic", "Cubic"]
        
        # Mapping requirements for dictionary validation
        self._req_map = {
            "Elastic": {"Elastic": ["AnalysisType", "Isotropy"]},
            "Plastic": {"Elastic": ["AnalysisType", "Isotropy"], "Plastic": ["Plasticity", "HardeningLaw", "sig_y0"]},
            "Transport": ["Dx", "Dy", "s"],
            "MechTrappingPFF": ["D0x", "D0y", "DQx", "DQy", "m", "s", "Vh", "zeta_rho", "Zd"],
            "2PhaseTrapping": ["D0x1", "D0y1", "DQx1", "DQy1", "D0x2", "D0y2", "zeta_j", "zeta_jj"]
        }

        # Logic initialization
        self.tag = self.ext_map[self.config.PhysicsType]
        self.nElementSets = len(self.mesh.materialNames)
        self.is_insulated = False
        
        self._validate_inputs()
        self._setup_physics_stats()
        self._print_summary()

    def _validate_inputs(self):
        """Unified validation for Physics, Categories, and Materials."""
        # Check Physics & Category
        if self.config.PhysicsType not in self.physicsRegistry:
            raise ValueError(f"Unknown PhysicsType: {self.config.PhysicsType}")
        
        if self.config.PhysicsCategory not in self.physicsRegistry[self.config.PhysicsType]:
            raise ValueError(f"Category {self.config.PhysicsCategory} invalid for {self.config.PhysicsType}")

        # Check Materials against req_map
        required_data = self._req_map.get(self.config.PhysicsCategory)
        if not required_data:
            return # Skip if no requirement defined for this category yet

        for mat_name in self.mesh.materialNames:
            props = self.materials.get(mat_name)
            if not props:
                raise KeyError(f"Material {mat_name} not found in input dictionary.")
            
            # Mechanical Nested Check
            if isinstance(required_data, dict):
                for group, keys in required_data.items():
                    if group not in props:
                        raise KeyError(f"Material '{mat_name}' missing sub-group: '{group}'")
                    for k in keys:
                        if k not in props[group]:
                            raise KeyError(f"Material '{mat_name}' -> '{group}' missing field: '{k}'")
                
                # Nested Logical Checks for Elasticity
                el = props["Elastic"]
                if el["AnalysisType"] not in self._allowed_Analysis:
                    raise ValueError(f"Invalid AnalysisType in {mat_name}")
                if el["Isotropy"] == "Isotropic" and not all(k in el for k in ["Emod", "nu"]):
                    raise KeyError(f"Isotropic material '{mat_name}' missing Emod/nu")
                if el["Isotropy"] == "Cubic" and not all(k in el for k in ["C11", "C12", "C44"]):
                    raise KeyError(f"Cubic material '{mat_name}' missing C11/C12/C44")
            else:
                # Flat check for Transport/PFF
                for key in required_data:
                    if key not in props:
                        raise KeyError(f"Material '{mat_name}' missing: {key}")

    def _setup_physics_stats(self):
        """Calculate DOFs and check simulation-wide constraints."""
        dof_map = {'Mechanical': self.mesh.nDim, 'Transport': 1, 'PFF': 1}
        self.nDofsPerNode = dof_map[self.config.PhysicsType]
        self.nTotDofs = self.mesh.nTotNodes * self.nDofsPerNode
        
        if self.config.PhysicsType == "Transport":
            if not self.config.dt or self.config.dt <= 0:
                raise ValueError("Transport requires a positive 'dt'.")
            self.is_insulated = len(self.config.presBCs) == 0

        if self.config.PhysicsType == "Mechanical" and not self.config.presBCs:
            raise ValueError("Mechanical simulation requires prescribed BCs (Supports).")

    def WriteInputFile(self, overwrite=True):
        filename = f"{self.config.SimulName}.{self.tag}.in.hdf5"
        mode = "w" if overwrite else "x"
        
        with h5py.File(filename, mode) as f:
            params = f.create_group("SimulationParameters")
            params.create_dataset("PhysicsType", data=np.bytes_(self.config.PhysicsType))
            params.create_dataset("Category", data=np.bytes_(self.config.PhysicsCategory))
            params.create_dataset("nDim", data=self.mesh.nDim, dtype=np.int64)
            params.create_dataset("nTotNodes", data=self.mesh.nTotNodes, dtype=np.int64)
            params.create_dataset("nTotDofs", data=self.nTotDofs, dtype=np.int64)
            params.create_dataset("nTotElements", data=self.mesh.nTotElements, dtype=np.int64)
            params.create_dataset("nPresDofs", data=len(self.config.presBCs), dtype=np.int64)
            params.create_dataset("nElementSets", data=self.nElementSets, dtype=np.int64)
            params.create_dataset("nSteps", data=self.config.nSteps, dtype=np.int64)
            params.create_dataset("R", data=self.config.R, dtype=np.float64)
            params.create_dataset("T", data=self.config.Temperature, dtype=np.float64)
            if self.config.dt: params.create_dataset("dt", data=self.config.dt, dtype=np.float64)
            if self.config.PhysicsCategory in ["2PhaseTrapping", "GBTrapping", "HLGBTrapping", "MechTrapping", "MechTrappingPFF"]:
                params.create_dataset("conB", data=self.config.conB, dtype=np.float64)
                params.create_dataset("nExitNodes", data=len(self.config.exitNodes), dtype=np.int64)
                f.create_dataset("ExitNodes", data=np.array(self.config.exitNodes, dtype=np.int64))
                
            # --- BC Logic ---
            if self.config.presBCs:
                paramsBC = f.create_group("PrescribedDOFs")
                for i, bc_val in enumerate(self.config.presBCs):
                    paramsBC.create_dataset(f"Prescribed_{i}", data=np.array(bc_val, dtype=np.float64))

            # Materials
            mat_grp = f.create_group("Materials")
            for i, name in enumerate(self.mesh.materialNames, start=1):
                m_sub = mat_grp.create_group(f"Material_{i}")
                self._write_dict_to_hdf5(m_sub, self.materials[name])
                
        print(f"  Input file initialized: {filename}")

    def _write_dict_to_hdf5(self, h5_group, data_dict):
        for key, val in data_dict.items():
            if isinstance(val, dict):
                self._write_dict_to_hdf5(h5_group.create_group(key), val)
            else:
                data = np.bytes_(val) if isinstance(val, str) else val
                h5_group.create_dataset(key, data=data)

    def WriteOutputFile(self, overwrite=True, **kwargs):
        """Initializes groups for the .out.hdf5 file."""
        filename = f"{self.config.SimulName}.{self.tag}.out.hdf5"
        mode = "w" if overwrite else "x"
        
        with h5py.File(filename, mode) as f:
            
            if self.config.PhysicsType == "Mechanical":
                for grp in ['Disp', 'Force', 'Stress', 'Strain']: f.create_group(grp)
                if self.config.PhysicsCategory == "Plastic":
                    for grp in ['Strain_e', 'Strain_p', 'Strain_eq', 'Stress_eq', 'Stress_h', 'Rho']: 
                        f.create_group(grp)
            
            elif self.config.PhysicsType == "PFF":
                f.create_group("Phi")
            
            elif self.config.PhysicsType == "Transport":
                f.create_group("Con")
                f.create_group("Time")
                if kwargs.get('AVCON'): f.create_group('AvCon')
                if kwargs.get('FLUX'): f.create_group('Flux')
                if kwargs.get('AVFLUX'):   f.create_group('AvFlux')
                if kwargs.get('TDS'):   f.create_group('Temp')
        
        print(f"  Output file initialized: {filename}")

    def _print_summary(self):
        print("\n" + "="*50)
        print(f"PHIMATS: {self.config.SimulName} | {self.config.PhysicsType} ({self.config.PhysicsCategory})")
        print("="*50)
        if self.is_insulated: print("  !! WARNING: INSULATED (No concentration BCs) !!")
        print(f"  Nodes/Elems : {self.mesh.nTotNodes} / {self.mesh.nTotElements}")
        print(f"  Total DOFs  : {self.nTotDofs}")
        print(f"  BCs Count   : {len(self.config.presBCs)}")
        print("="*50 + "\n")