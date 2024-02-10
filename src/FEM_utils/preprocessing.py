
import numpy as numpy
from pathlib import Path
import h5py
import numpy as np

class PreProcessing:
    
    def __init__(self, inputData):
        
        self.Simul = inputData["Simul"]
        self.nNodes = inputData["nNodes"]
        self.nElements = inputData["nElements"]
        self.elementName = inputData["elementName"]
        self.nElementSets = inputData["nElementSets"]
        self.nPresDofs = inputData["nPresDofs"]
        self.presDOFs = inputData["presDOFs"]

        allowedElements = ["quad4", "quad8"]
        elements2D = ["quad4"]
        elements3D = ["quad8"]
        if not self.elementName in allowedElements:
            ErrString = "ERROR! Unknown element type: " + self.elementType
            raise Exception(ErrString)
        
        if self.elementName == "quad4":
            self.elementNodes = 4
            self.elementType = 1
            
        if self.elementName in elements2D:
            self.nDim = 2
        
        self.nodeConnectivity = inputData["nodeConnectivity"]
        self.nodeCoord = inputData["nodeCoord"]
        
        self.Emod = inputData["Emod"]
        self.nu = inputData["nu"]

        
        pass
    
#-----------------------------------------------------------------------------#
    
    def CreateFileHDF5(self):
        
        # path = Path(self.Simul+'.hdf5')
        
        # if path.is_file():
        #     ErrString = 'File '+ self.Simul+ '.hdf5 exists. Can not override!'
        #     raise Exception(ErrString)
        
        self.fh5 = h5py.File(self.Simul+"_in.hdf5", "w")
        
        # SimulationParameters -------
        self.fh5.attrs["Simulation"] =  self.Simul
        self.grp_Sim_Params = self.fh5.create_group('SimulationParameters')
        self.grp_Sim_Params.create_dataset("nNodes", data=self.nNodes, dtype = np.int64)
        self.grp_Sim_Params.create_dataset("nDim", data=self.nDim, dtype = np.int64)
        self.grp_Sim_Params.create_dataset("nElements", data=self.nElements, dtype = np.int64)
        self.grp_Sim_Params.create_dataset("nPresDofs", data=self.nPresDofs, dtype = np.int64)
        self.grp_Sim_Params.create_dataset("elementType", data=self.elementType, dtype = np.int64) 
        self.grp_Sim_Params.create_dataset("nElementSets", data=self.nElementSets, dtype = np.int64) 
        self.grp_Sim_Params.create_dataset("elementNodes", data=self.elementNodes, dtype = np.int64) 

        self.grp_Sim_Params.create_dataset("Emod", data=self.Emod) 
        self.grp_Sim_Params.create_dataset("nu", data=self.nu) 
        
        # Write node coordinates
        self.grp_nNodes = self.fh5.create_group('NodeCoordinates')
        for iset in range(self.nNodes):
            self.grp_nNodes.create_dataset("Node_"+str(iset), data=self.nodeCoord[iset], dtype = np.float64)   
        
        # Write node connectivity
        self.grp_nodeConnectivity = self.fh5.create_group('NodeConnectivity')
        for iset in range(self.nElements):
            self.grp_nodeConnectivity.create_dataset("Element_"+str(iset), data=self.nodeConnectivity[iset], dtype = np.int64)      
            
        # Write prescribed dofs
        self.grp_prescribedDOFs = self.fh5.create_group('PrescribedDOFs')
        for iPreDof in range(self.nPresDofs):
            self.grp_prescribedDOFs.create_dataset("Prescribed_"+str(iPreDof), data=self.presDOFs[iPreDof]) 
       
       # Close hdf5 file
        self.CloseFileHDF5()
        
        pass

#-----------------------------------------------------------------------------#
    
    def CloseFileHDF5(self):
        
        self.fh5.close()
        
        pass
    