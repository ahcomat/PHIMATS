import h5py
import numpy as np
import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom

class WriteXDMF:
    def __init__(self, fileName, elementName, nSteps, components, nDim, 
                 mechModel="Elastic", START=0, tOut=1, FLUX=False, TDS=False):
        """
        Args:
            fileName (str): Base name of the simulation.
            elementName (str): 'triangle', 'quad', etc.
            nSteps (int): Number of time steps.
            components (list): Active physics, e.g., ["mech", "diff", "pff"].
            nDim (int): 2 or 3.
            mechModel (str): "Elastic", "Plastic", "J2", or "CP". Defaults to "Elastic".
            START (int): Initial step number. Defaults to 0.
            tOut (int): Number of steps to tOut when writing outputs. Defaults to 1.
            FLUX (bool): Flag for flux field. Defaults to False.
            TDS (bool): Flag for writing Temp in TDS. 
        """
        self.FName = fileName
        self.nDim = nDim
        self.components = components
        self.mechModel = mechModel
        
        # Standardized file map 
        self.mesh_file = f"{fileName}.mesh.hdf5"
        self.file_map = {
            "mech": f"{fileName}.mech.out.hdf5",
            "diff": f"{fileName}.diff.out.hdf5",
            "pff":  f"{fileName}.pff.out.hdf5"
        }

        # Element configuration
        self.element_config = {
            "triangle":   {"ElTopology": "Triangle", "nElNodes": 3, "nElStres": 3},
            "quad":       {"ElTopology": "Quadrilateral", "nElNodes": 4, "nElStres": 3},
            "hexahedron": {"ElTopology": "Hexahedron", "nElNodes": 8, "nElStres": 6},
        }
        
        conf = self.element_config[elementName]
        self.ElTopology = conf["ElTopology"]
        self.nElNodes   = conf["nElNodes"]
        self.nElStres   = conf["nElStres"]

        # Read mesh metadata from the dedicated mesh file
        try:
            with h5py.File(self.mesh_file, "r") as f:
                self.nTotNodes = len(f["NodeCoordinates"])
                self.nTotElements = len(f["NodeConnectivity"])
        except OSError:
            raise OSError(f"Could not find mesh file: {self.mesh_file}")

        # Initialize XDMF
        root = ET.Element("Xdmf", Version="3.0")
        domain = ET.SubElement(root, "Domain")
        self.time_series_grid = ET.SubElement(domain, "Grid", Name="TimeSeries", GridType="Collection", CollectionType="Temporal")

        # Build time steps
        for t in range(START, nSteps, tOut):
            grid = self._create_step_header(t)
            
            # --- TRANSPORT COMPONENTS ---
            if "diff" in self.components:
                f_out = self.file_map["diff"]
                self._add_attr(grid, t, "Con", "Scalar", f_out)
                if FLUX: self._add_attr(grid, t, "Flux", "Vector", f_out)
                if TDS:  self._add_attr(grid, t, "Temp", "Scalar", f_out)
            
            # --- PHASE FIELD FRACTURE COMPONENTS ---
            if "pff" in self.components:
                self._add_attr(grid, t, "Phi", "Scalar", self.file_map["pff"])
            
            # --- MECHANICAL COMPONENTS ---
            if "mech" in self.components:
                f_out = self.file_map["mech"]
                self._add_attr(grid, t, "Disp", "Vector", f_out)
                self._add_attr(grid, t, "Force", "Vector", f_out)
                self._add_attr(grid, t, "Stress", "Tensor", f_out)
                self._add_attr(grid, t, "Strain", "Tensor", f_out)
                
                # Add plasticity specifics
                if self.mechModel in ["Plastic", "J2", "CP"]:
                    self._add_attr(grid, t, "Strain_e", "Tensor", f_out)
                    self._add_attr(grid, t, "Strain_p", "Tensor", f_out)
                    self._add_attr(grid, t, "Strain_eq", "Scalar", f_out)
                    self._add_attr(grid, t, "Stress_eq", "Scalar", f_out)
                    self._add_attr(grid, t, "Stress_h", "Scalar", f_out)

        # Write output
        with open(f"{fileName}.xdmf", "w") as f:
            f.write(self.prettify(root))

    def _create_step_header(self, t):
        step_grid = ET.SubElement(self.time_series_grid, "Grid", Name=f"Step_{t}")
        ET.SubElement(step_grid, "Time", Value=str(float(t)))

        # Topology (from .mesh.hdf5)
        topo = ET.SubElement(step_grid, "Topology", TopologyType=self.ElTopology, NumberOfElements=str(self.nTotElements))
        ET.SubElement(topo, "DataItem", Format="HDF", Dimensions=f"{self.nTotElements} {self.nElNodes}").text = f"{self.mesh_file}:/NodeConnectivity"

        # Geometry (from .mesh.hdf5)
        geo_type = "XYZ" if self.nDim == 3 else "XY"
        geom = ET.SubElement(step_grid, "Geometry", GeometryType=geo_type)
        ET.SubElement(geom, "DataItem", Format="HDF", Dimensions=f"{self.nTotNodes} {self.nDim}").text = f"{self.mesh_file}:/NodeCoordinates"
        
        return step_grid

    def _add_attr(self, grid, t, name, attr_type, h5_path):
        dim_str = str(self.nTotNodes)
        if attr_type == "Vector": dim_str += f" {self.nDim}"
        if attr_type == "Tensor": dim_str += f" {self.nElStres}"

        attr = ET.SubElement(grid, "Attribute", Name=name, AttributeType=attr_type, Center="Node")
        h5_internal = f"/{name}/Step_{t}"
        ET.SubElement(attr, "DataItem", Format="HDF", Dimensions=dim_str).text = f"{h5_path}:{h5_internal}"

    def prettify(self, elem):
        rough_string = ET.tostring(elem, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        return reparsed.toprettyxml(indent="  ")