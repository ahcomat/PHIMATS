import h5py
import numpy as np
import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom


class WriteXDMF:
    
    def __init__(self, fileName, elementName, nSteps, simulationType, skip=1):
        """A class to write xdmf file for paraview visualization of the results in _out.hdf5 file.

        Args:
            inputData (_type_): _description_
        """
        
        # Name of file
        self.FName = fileName     
        
        # Read mesh parameters from _in.hdf5        
        try:
            fh5 = h5py.File(fileName+"_in.hdf5", "r")

        except OSError as e:
            print(f"Error opening file {fileName}: {e}")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

        self.nTotNodes = fh5["SimulationParameters"]["nTotNodes"][()]
        self.nTotDofs = fh5["SimulationParameters"]["nTotDofs"][()]
        self.nDim = fh5["SimulationParameters"]["nDim"][()]
        self.nTotElements = fh5["SimulationParameters"]["nTotElements"][()]
        self.nSteps = fh5["SimulationParameters"]["nSteps"][()]

        fh5.close()  
        
        # ---------------------------------------- #   
        
        # Read element type
        # Allowed elements (naming should match with meshio)
        allowedElements = ["quad", "quad8", "triangle", "triangle6", "hexahedron"]
        
        if not elementName in allowedElements:
            ErrString = "ERROR! Unknown element name < " + elementName + " >\n"
            ErrString += "Allowed elements are: \n"
            for elem in allowedElements:
                ErrString += elem + "\n"
            raise ValueError(ErrString)
        
        if elementName=="triangle":
            self.ElTopology = "Triangle"
            self.nElNodes = 3
        elif elementName=="quad":
            self.ElTopology = "Quadrilateral"
            self.nElNodes = 4
            self.nElStres = 3
        elif elementName=="":
            self.ElTopology = "Hexahedron"
            self.nElNodes = 8
            self.nElStres = 6
            
        # ---------------------------------------- #  
        
        allowedSimTypes = ["Transport2D", "Elastic2D", "Elastic3D", "Plastic2D", "Plastic3D"] 
        
        self.simType = simulationType
        
        if not self.simType in allowedSimTypes:
            ErrString = "ERROR! Unknown simulation type < " + self.simType + " >\n"
            ErrString += "Allowed simulation types are: \n"
            for elem in allowedSimTypes:
                ErrString += elem + "\n"
            raise ValueError(ErrString)
        
        # ---------------------------------------- #   
        
        root = ET.Element("Xdmf", Version="3.0", xmlns="http://www.w3.org/2001/XInclude")
        domain = ET.SubElement(root, "Domain")
        self.time_series_grid = ET.SubElement(domain, "Grid", Name="TimeSeries", GridType="Collection", CollectionType="Temporal")  
        
        for t in range(0, nSteps+1, skip):
            
            if self.simType=="Transport2D":
                self.WriteCon2D(t)
            elif self.simType=="Plastic2D":
                self.WritePlastic2D(t)
                
        # ---------------------------------------- #   
                
        # Write the pretty-printed XML to a file
        with open(self.FName+".xdmf", "w") as f:
            f.write(self.prettify(root))
            
    pass

    def prettify(self, elem):
        
        rough_string = ET.tostring(elem, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        
        return reparsed.toprettyxml(indent="  ")
    
    def WriteCon2D(self, t):
        
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
        
        pass
    
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
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(self.nElStres)).text = self.FName+"_out.hdf5:/Strain_e/Step_"+str(t)# Add attribute element for stress
        
        # Add attribute element for Strain_eq
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="strain_eq", AttributeType="Scalar", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(1)).text = self.FName+"_out.hdf5:/Strain_eq/Step_"+str(t)
        
        # Add attribute element for Stress_eq
        attribute = ET.SubElement(timestep_grid, "Attribute", Name="stress_eq", AttributeType="Scalar", Center="Node")
        ET.SubElement(attribute, "DataItem", Format="HDF", DataType="Float", Dimensions=str(self.nTotNodes)+" "+str(1)).text = self.FName+"_out.hdf5:/Stress_eq/Step_"+str(t)
        
        pass
    
    def WritePlastic2D(self, t):
        
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
        
        pass
