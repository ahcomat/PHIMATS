import h5py
import numpy as np
import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom


class WriteXDMF:
    
    def __init__(self, fileName, elementName, nSteps, simulationType, START=0, FLUX=False, skip=1):
        """
        A class to write XDMF files for ParaView visualization of results stored in _out.hdf5.

        Args:
            fileName (str): Base name for the XDMF file.
            elementName (str): Mesh element type (meshio convention).
            nSteps (int): Number of time steps.
            simulationType (str): Type of simulation. Allowed: ["Transport2D", "Elastic2D", "Elastic3D", "Plastic2D", "Plastic3D"].
            START (int, optional): Initial step number. Defaults to 1.
            FLUX (bool, optional): Flag for flux field. Defaults to False.
            skip (int, optional): Number of steps to skip when writing outputs. Defaults to 1.
        """
        self.FName = fileName  # Base file name

        # Read mesh parameters from HDF5
        try:
            with h5py.File(f"{fileName}_in.hdf5", "r") as fh5:
                self.nTotNodes = fh5["SimulationParameters"]["nTotNodes"][()]
                self.nTotDofs = fh5["SimulationParameters"]["nTotDofs"][()]
                self.nDim = fh5["SimulationParameters"]["nDim"][()]
                self.nTotElements = fh5["SimulationParameters"]["nTotElements"][()]
        except OSError as e:
            raise OSError(f"Error opening file '{fileName}_in.hdf5': {e}")
        
        self.nSteps = nSteps

        # Validate element type
        self.element_config = {
            "triangle": {"ElTopology": "Triangle", "nElNodes": 3},
            "triangle6": {"ElTopology": "Triangle", "nElNodes": 6},
            "quad": {"ElTopology": "Quadrilateral", "nElNodes": 4, "nElStres": 3},
            "hexahedron": {"ElTopology": "Hexahedron", "nElNodes": 8, "nElStres": 6},
        }
        
        if elementName not in self.element_config:
            raise ValueError(f"Unknown element name '{elementName}'. Allowed: {', '.join(self.element_config.keys())}")
        
        self.ElTopology = self.element_config[elementName]["ElTopology"]
        self.nElNodes = self.element_config[elementName]["nElNodes"]
        self.nElStres = self.element_config[elementName].get("nElStres")

        # Validate simulation type
        allowedSimTypes = ["Transport2D", "Elastic2D", "Elastic3D", "Plastic2D", "Plastic3D"]
        if simulationType not in allowedSimTypes:
            raise ValueError(f"Unknown simulation type '{simulationType}'. Allowed: {', '.join(allowedSimTypes)}")
        
        self.simType = simulationType

        # Initialize XDMF structure
        root = ET.Element("Xdmf", Version="3.0", xmlns="http://www.w3.org/2001/XInclude")
        domain = ET.SubElement(root, "Domain")
        self.time_series_grid = ET.SubElement(domain, "Grid", Name="TimeSeries", GridType="Collection", CollectionType="Temporal")

        # Write time steps
        for t in range(START, self.nSteps, skip):
            if self.simType == "Transport2D":
                self.WriteCon2D(t, FLUX)
            elif self.simType == "Elastic2D":
                self.WriteElastic2D(t)
            elif self.simType == "Plastic2D":
                self.WritePlastic2D(t)
            elif self.simType == "Plastic3D":
                self.WritePlastic2D(t)
                
        # ---------------------------------------- #   
                
        # Write the pretty-printed XML to a file
        with open(self.FName+".xdmf", "w") as f:
            f.write(self.prettify(root))
            
#-----------------------------------------------------------------------------#

    def prettify(self, elem):
        """Returns a pretty-printed XML string for the given element."""
        
        rough_string = ET.tostring(elem, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        
        return reparsed.toprettyxml(indent="  ")
    
#-----------------------------------------------------------------------------#
    
    def WriteCon2D(self, t, FLUX):
        
        timestep_grid = ET.SubElement(self.time_series_grid, "Grid", Name=f"TimeStep_{t:.1f}")
        ET.SubElement(timestep_grid, "Time", Value=f"{t:.1f}")

        # Add topology element
        topology = ET.SubElement(timestep_grid, "Topology", TopologyType=self.ElTopology, NumberOfElements=str(self.nTotElements))
        ET.SubElement(topology, "DataItem", Format="HDF", DataType="Int", Dimensions=str(self.nTotElements)+" "+str(self.nElNodes)).text = self.FName+"_in.hdf5:/NodeConnectivity"

        # Add geometry element
        geometry = ET.SubElement(timestep_grid, "Geometry", GeometryType="XY")
        ET.SubElement(geometry, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(self.nDim)).text = self.FName+"_in.hdf5:/NodeCoordinates"

        # Add attribute element for concentration
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="con", AttributeType="Scalar", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)).text = self.FName+"_out.hdf5:/Con/Step_"+str(t)
        
        if FLUX:
            # Add attribute element for flux
            attribute = ET.SubElement(timestep_grid, "Attribute", Name="flux", AttributeType="Vector", Center="Node")
            ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(self.nDim)).text = self.FName+"_out.hdf5:/Flux/Step_"+str(t)
            
#-----------------------------------------------------------------------------#

    def WriteElastic2D(self, t):
        
        timestep_grid = ET.SubElement(self.time_series_grid, "Grid", Name=f"TimeStep_{t:.1f}")
        ET.SubElement(timestep_grid, "Time", Value=f"{t:.1f}")

        # Add topology element
        topology = ET.SubElement(timestep_grid, "Topology", TopologyType=self.ElTopology, NumberOfElements=str(self.nTotElements))
        ET.SubElement(topology, "DataItem", Format="HDF", DataType="Int", Dimensions=str(self.nTotElements)+" "+str(self.nElNodes)).text = self.FName+"_in.hdf5:/NodeConnectivity"

        # Add geometry element
        geometry = ET.SubElement(timestep_grid, "Geometry", GeometryType="XY")
        ET.SubElement(geometry, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(self.nDim)).text = self.FName+"_in.hdf5:/NodeCoordinates"

        # Add attribute element for displacement
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="disp", AttributeType="Vector", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(self.nDim)).text = self.FName+"_out.hdf5:/Disp/Step_"+str(t)
        
        # Add attribute element for force
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="force", AttributeType="Vector", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(self.nDim)).text = self.FName+"_out.hdf5:/Force/Step_"+str(t)
        
        # Add attribute element for strain
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="strain", AttributeType="Tensor", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(self.nElStres)).text = self.FName+"_out.hdf5:/Strain/Step_"+str(t)
        
        # Add attribute element for stress
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="stress", AttributeType="Tensor", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(self.nElStres)).text = self.FName+"_out.hdf5:/Stress/Step_"+str(t)
                
#-----------------------------------------------------------------------------#

    def WritePlastic2D(self, t):
        
        timestep_grid = ET.SubElement(self.time_series_grid, "Grid", Name=f"TimeStep_{t:.1f}")
        ET.SubElement(timestep_grid, "Time", Value=f"{t:.1f}")

        # Add topology element
        topology = ET.SubElement(timestep_grid, "Topology", TopologyType=self.ElTopology, NumberOfElements=str(self.nTotElements))
        ET.SubElement(topology, "DataItem", Format="HDF", DataType="Int", Dimensions=str(self.nTotElements)+" "+str(self.nElNodes)).text = self.FName+"_in.hdf5:/NodeConnectivity"

        # Add geometry element
        geometry = ET.SubElement(timestep_grid, "Geometry", GeometryType="XY")
        ET.SubElement(geometry, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(self.nDim)).text = self.FName+"_in.hdf5:/NodeCoordinates"

        # Add attribute element for displacement
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="disp", AttributeType="Vector", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(self.nDim)).text = self.FName+"_out.hdf5:/Disp/Step_"+str(t)
        
        # Add attribute element for force
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="force", AttributeType="Vector", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(self.nDim)).text = self.FName+"_out.hdf5:/Force/Step_"+str(t)
        
        # Add attribute element for strain
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="strain", AttributeType="Tensor", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(self.nElStres)).text = self.FName+"_out.hdf5:/Strain/Step_"+str(t)
        
        # Add attribute element for stress
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="stress", AttributeType="Tensor", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(self.nElStres)).text = self.FName+"_out.hdf5:/Stress/Step_"+str(t)
        
        # Add attribute element for Strain_p
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="strain_p", AttributeType="Tensor", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(self.nElStres)).text = self.FName+"_out.hdf5:/Strain_p/Step_"+str(t)
        
        # Add attribute element for Strain_e
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="strain_e", AttributeType="Tensor", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(self.nElStres)).text = self.FName+"_out.hdf5:/Strain_e/Step_"+str(t)
        
        # Add attribute element for Strain_eq
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="strain_eq", AttributeType="Scalar", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(1)).text = self.FName+"_out.hdf5:/Strain_eq/Step_"+str(t)
        
        # Add attribute element for Stress_eq
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="stress_eq", AttributeType="Scalar", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(1)).text = self.FName+"_out.hdf5:/Stress_eq/Step_"+str(t)
        
        # Add attribute element for hydrostatic stress
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="stress_h", AttributeType="Scalar", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(1)).text = self.FName+"_out.hdf5:/Stress_h/Step_"+str(t)
        
    #-----------------------------------------------------------------------------#
    
    def WritePlastic3D(self, t):
        
        timestep_grid = ET.SubElement(self.time_series_grid, "Grid", Name=f"TimeStep_{t:.1f}")
        ET.SubElement(timestep_grid, "Time", Value=f"{t:.1f}")

        # Add topology element
        topology = ET.SubElement(timestep_grid, "Topology", TopologyType=self.ElTopology, NumberOfElements=str(self.nTotElements))
        ET.SubElement(topology, "DataItem", Format="HDF", DataType="Int", Dimensions=str(self.nTotElements)+" "+str(self.nElNodes)).text = self.FName+"_in.hdf5:/NodeConnectivity"

        # Add geometry element
        geometry = ET.SubElement(timestep_grid, "Geometry", GeometryType="XYZ")
        ET.SubElement(geometry, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(self.nDim)).text = self.FName+"_in.hdf5:/NodeCoordinates"

        # Add attribute element for displacement
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="disp", AttributeType="Vector", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(self.nDim)).text = self.FName+"_out.hdf5:/Disp/Step_"+str(t)
        
        # Add attribute element for force
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="force", AttributeType="Vector", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(self.nDim)).text = self.FName+"_out.hdf5:/Force/Step_"+str(t)
        
        # Add attribute element for strain
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="strain", AttributeType="Tensor", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(self.nElStres)).text = self.FName+"_out.hdf5:/Strain/Step_"+str(t)
        
        # Add attribute element for strain_e
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="strain_e", AttributeType="Tensor", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(self.nElStres)).text = self.FName+"_out.hdf5:/Strain_e/Step_"+str(t)
        
        # Add attribute element for strain_p
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="strain_p", AttributeType="Tensor", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(self.nElStres)).text = self.FName+"_out.hdf5:/Strain_p/Step_"+str(t)
        
        # Add attribute element for stress
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="stress", AttributeType="Tensor", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(self.nElStres)).text = self.FName+"_out.hdf5:/Stress/Step_"+str(t)
        
        # Add attribute element for equivalent stress
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="stress_eq", AttributeType="Scalar", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)).text = self.FName+"_out.hdf5:/Stress_eq/Step_"+str(t)
        
        # Add attribute element for equivalent strain
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="strain_eq", AttributeType="Scalar", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)).text = self.FName+"_out.hdf5:/Strain_eq/Step_"+str(t)
        
        # Add attribute element for hydostatic stress
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="stress_h", AttributeType="Scalar", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(1)).text = self.FName+"_out.hdf5:/Stress_h/Step_"+str(t)
        
#-----------------------------------------------------------------------------#

