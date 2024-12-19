import numpy as numpy
import meshio
from pathlib import Path
import h5py
import numpy as np


class PreProcessing:
    
    def __init__(self, inputData):
                
        #----------------------------------------------------------------------
        # Read input data  
        #----------------------------------------------------------------------
        
        self.Simul = inputData["Simul"]
        self.SimulType = inputData["SimulType"]
        self.mesh = inputData["mesh"]
        self.elementName = inputData["elementName"]
        self.nElementSets = inputData["nElementSets"]
        self.presBCs = inputData["presBCs"]
        self.nPresDofs = len(self.presBCs)
        self.nSteps = inputData["nSteps"]
        
        # Allowed simulation types
        allowedSimulTypes = ["Mechanical", "Transport", "PhaseTrapping", "GBTrapping", "MechTrapping"]
        
        self.TransportSimulTypes = ["Transport", "PhaseTrapping", "GBTrapping", "MechTrapping"]
                
        if not self.SimulType in allowedSimulTypes:
            ErrString = "ERROR! Unknown simulation type < " + self.SimulType + " >\n"
            ErrString += "Allowed elements are: \n"
            for simulType in allowedSimulTypes:
                ErrString += simulType + "\n"
            raise ValueError(ErrString)
        
        if self.SimulType=="Transport":
            self.exitNods = inputData["exitNods"]
            self.dt = inputData["dt"]
        elif self.SimulType=="GBTrapping":
            self.exitNods = inputData["exitNods"]
            self.dt = inputData["dt"]
            self.T = inputData["T"]
            self.gPhi = inputData["gPhi"]
        elif self.SimulType=="PhaseTrapping":
            self.exitNods = inputData["exitNods"]
            self.dt = inputData["dt"]
            self.T = inputData["T"]
            self.martensite = inputData["martensite"]
            self.gPhi_MM = inputData["gPhi_MM"]
            self.gPhi_ff = inputData["gPhi_ff"]
            self.gPhi_fM = inputData["gPhi_fM"]
        elif self.SimulType=="MechTrapping":
            self.exitNods = inputData["exitNods"]
            self.dt = inputData["dt"]
            self.R = inputData["R"]
            self.T = inputData["T"]
            self.sigmaH = inputData["sigmaH"]
        
        #----------------------------------------------------------------------
        # Check for allowed elements and assign element data (number of nodes,
        # dimension and order) 
        #----------------------------------------------------------------------

        # Allowed elements (naming should match with meshio)
        allowedElements = ["quad", "quad8", "triangle", "triangle6", "hexahedron"]
        # 2D elements
        elements2D = ["quad", "quad8", "triangle", "triangle6"]
        # 3D elements
        elements3D = ["hexahedron"]
        # First order
        elementsOrder1 = ["quad", "triangle", "hexahedron"]
        # Second order
        elementsOrder2 = ["quad8", "triangle6"]
        
        if not self.elementName in allowedElements:
            ErrString = "ERROR! Unknown element name < " + self.elementName + " >\n"
            ErrString += "Allowed elements are: \n"
            for elem in allowedElements:
                ErrString += elem + "\n"
            raise ValueError(ErrString)
        
        # 2D or 3D
        if self.elementName in elements2D:
            self.nDim = 2
        if self.elementName in elements3D:
            self.nDim = 3
        
        # First or second order
        if self.elementName in elementsOrder1:
            self.nOrder = 1
        if self.elementName in elementsOrder2:
            self.nOrder = 2
            
        #----------------------------------------------------------------------
        # Read nodes  
        #----------------------------------------------------------------------
        
        # Total number of nodes
        self.nTotNodes = self.mesh.points.shape[0]  
        # Node connectivity
        self.nodeConnectivity = self.mesh.cells_dict[self.elementName] 
        
        # Node coordinates
        if self.nDim == 2:
            self.nodeCoord = self.mesh.points[:,0:2]
        elif self.nDim == 3:    
            self.nodeCoord = self.mesh.points
                    
        # Total number of Dofs
        if self.SimulType == "Mechanical":
            self.nTotDofs = self.nTotNodes*self.nDim
        elif self.SimulType == "Transport":
            self.nTotDofs = self.nTotNodes
        elif self.SimulType == "GBTrapping":
            self.nTotDofs = self.nTotNodes
        elif self.SimulType == "PhaseTrapping":
            self.nTotDofs = self.nTotNodes
        elif self.SimulType == "MechTrapping":
            self.nTotDofs = self.nTotNodes
        
        #----------------------------------------------------------------------
        # Read number of elements  
        #----------------------------------------------------------------------
        
        # Total number of elements
        self.nTotElements = self.mesh.cells_dict[self.elementName].data.shape[0]  

        #----------------------------------------------------------------------
        # Read material data
        #----------------------------------------------------------------------
        
        # Read materials dict
        self.Materials = inputData["Materials"]
        
        # Check isotropy
        Isotropy = ["Isotropic", "Cubic"]
        
        # Check plasticity type
        HardeningLaws = ["Linear", "PowerLaw"]
                
        for mat in self.Materials:
            if "Elastic" in self.Materials[mat]:
                if not self.Materials[mat]["Elastic"]["Isotropy"] in Isotropy:
                    ErrString = "ERROR! undefined hardening law < " + self.Materials[mat]["Elastic"]["Isotropy"] + " >\n"
                    ErrString += "Allowed hardening laws are: \n"
                    for i in HardeningLaws:
                        ErrString += i + "\n"
                    raise ValueError(ErrString)
                
            if "Plastic" in self.Materials[mat]:
                if not self.Materials[mat]["Plastic"]["HardeningLaw"] in HardeningLaws:
                    ErrString = "ERROR! undefined hardening law < " + self.Materials[mat]["Plastic"]["HardeningLaw"] + " >\n"
                    ErrString += "Allowed hardening laws are: \n"
                    for i in HardeningLaws:
                        ErrString += i + "\n"
                    raise ValueError(ErrString)

        pass
    
#-----------------------------------------------------------------------------#
    
    def WriteFileHDF5(self):
        """
        Creates and writes data to hdf5 input file
        """
        
        # Open hdf5 file 
        self.OpenFileHDF5()
        
        # try :
        
        #----------------------------------------------------------------------
        # Write simulation parameters to hdf5  
        #----------------------------------------------------------------------
            
        self.fh5.attrs["Simulation"] =  self.Simul
                    
        #----------------------------------------------------------------------
        # Nodes/Elements data 
        #----------------------------------------------------------------------
        
        self.grp_Sim_Params = self.fh5.create_group('SimulationParameters')
        
        if self.SimulType == "GBTrapping":
            self.grp_Sim_Params.create_dataset("Trapping", data=int(1), dtype = np.int64)
        elif self.SimulType == "PhaseTrapping":
            self.grp_Sim_Params.create_dataset("Trapping", data=int(2), dtype = np.int64)
        
        
        self.grp_Sim_Params.create_dataset("nDim", data=self.nDim, dtype = np.int64)
        self.grp_Sim_Params.create_dataset("nTotNodes", data=self.nTotNodes, dtype = np.int64)
        self.grp_Sim_Params.create_dataset("nTotDofs", data=self.nTotDofs, dtype = np.int64)
        self.grp_Sim_Params.create_dataset("nTotElements", data=self.nTotElements, dtype = np.int64)
        self.grp_Sim_Params.create_dataset("nPresDofs", data=self.nPresDofs, dtype = np.int64)
        self.grp_Sim_Params.create_dataset("nElementSets", data=self.nElementSets, dtype = np.int64)
        self.grp_Sim_Params.create_dataset("nSteps", data=self.nSteps, dtype = np.int64)
        
        if self.SimulType == "Transport":
            self.grp_Sim_Params.create_dataset("dt", data=self.dt)
            
        if self.SimulType in self.TransportSimulTypes:
            self.grp_Sim_Params.create_dataset("nExitNodes", data=len(self.exitNods))
        
        # Case Trapping 
        if self.SimulType == "GBTrapping":
            self.grp_Sim_Params.create_dataset("dt", data=self.dt)
            self.grp_Sim_Params.create_dataset("T", data=self.T)
            self.fh5.create_dataset('gPhi', data=self.gPhi, dtype = np.float64)
            
        # Case Trapping 
        if self.SimulType == "PhaseTrapping":
            self.grp_Sim_Params.create_dataset("dt", data=self.dt)
            self.grp_Sim_Params.create_dataset("T", data=self.T)
            self.fh5.create_dataset('gPhi_MM', data=self.gPhi_MM, dtype = np.float64) 
            self.fh5.create_dataset('gPhi_fM', data=self.gPhi_fM, dtype = np.float64) 
            self.fh5.create_dataset('gPhi_ff', data=self.gPhi_ff, dtype = np.float64) 
            self.fh5.create_dataset('martensite', data=self.martensite, dtype = np.float64) 
            
        # Case MechTrapping 
        if self.SimulType == "MechTrapping":
            self.grp_Sim_Params.create_dataset("dt", data=self.dt)
            self.grp_Sim_Params.create_dataset("R", data=self.R)
            self.grp_Sim_Params.create_dataset("T", data=self.T)
            self.fh5.create_dataset('sigmaH', data=self.sigmaH, dtype = np.float64)  
        
        #----------------------------------------------------------------------
        # Material data
        #----------------------------------------------------------------------
        self.grp_Materials = self.fh5.create_group('Materials')
        
        if self.SimulType == "Mechanical":
        
            counter = 0
            for mat in self.Materials:
                counter+=1
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/Elastic/AnalysisType", data=np.bytes_(self.Materials[mat]['Elastic']["AnalysisType"]))
                if self.Materials[mat]['Elastic']['Isotropy'] == "Isotropic":
                    self.grp_Materials.create_dataset("Material_"+str(counter)+"/Elastic/Isotropy", data=np.bytes_(self.Materials[mat]['Elastic']["Isotropy"]))
                    self.grp_Materials.create_dataset("Material_"+str(counter)+"/Elastic/Emod", data=self.Materials[mat]['Elastic']["Emod"])
                    self.grp_Materials.create_dataset("Material_"+str(counter)+"/Elastic/nu", data=self.Materials[mat]['Elastic']["nu"])
                elif self.Materials[mat]['Elastic']['Isotropy'] == "Cubic":
                    self.grp_Materials.create_dataset("Material_"+str(counter)+"/Elastic/Isotropy", data=np.bytes_(self.Materials[mat]['Elastic']["Isotropy"]))
                    self.grp_Materials.create_dataset("Material_"+str(counter)+"/Elastic/C11", data=self.Materials[mat]['Elastic']["C11"])
                    self.grp_Materials.create_dataset("Material_"+str(counter)+"/Elastic/C12", data=self.Materials[mat]['Elastic']["C12"])
                    self.grp_Materials.create_dataset("Material_"+str(counter)+"/Elastic/C44", data=self.Materials[mat]['Elastic']["C44"])

                # Check plasticity
                if "Plastic" in self.Materials[mat]:
                    self.grp_Materials.create_dataset("Material_"+str(counter)+"/Plastic/Plasticity", data=np.bytes_(self.Materials[mat]['Plastic']["Plasticity"]))
                    if self.Materials[mat]['Plastic']["HardeningLaw"] == "Linear":
                        self.grp_Materials.create_dataset("Material_"+str(counter)+"/Plastic/HardeningLaw", data=np.bytes_(self.Materials[mat]['Plastic']["HardeningLaw"]))
                        self.grp_Materials.create_dataset("Material_"+str(counter)+"/Plastic/sig_y0", data=self.Materials[mat]['Plastic']["sig_y0"])
                        self.grp_Materials.create_dataset("Material_"+str(counter)+"/Plastic/K_hard", data=self.Materials[mat]['Plastic']["K_hard"])
                    elif self.Materials[mat]['Plastic']["HardeningLaw"] == "PowerLaw":
                        self.grp_Materials.create_dataset("Material_"+str(counter)+"/Plastic/HardeningLaw", data=np.bytes_(self.Materials[mat]['Plastic']["HardeningLaw"]))
                        self.grp_Materials.create_dataset("Material_"+str(counter)+"/Plastic/sig_y0", data=self.Materials[mat]['Plastic']["sig_y0"])
                        self.grp_Materials.create_dataset("Material_"+str(counter)+"/Plastic/K_hard", data=self.Materials[mat]['Plastic']["K_hard"])
                        self.grp_Materials.create_dataset("Material_"+str(counter)+"/Plastic/n_pow", data=self.Materials[mat]['Plastic']["n_pow"])                   

        elif self.SimulType == "Transport":
            
            counter = 0
            for mat in self.Materials:
                counter+=1
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/Dx", data=self.Materials[mat]["Dx"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/Dy", data=self.Materials[mat]["Dy"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/s", data=self.Materials[mat]["s"])

                if self.nDim == 3:
                    self.grp_Materials.create_dataset("Material_"+str(counter)+"/Dy", data=self.Materials[mat]["Dz"])
                    
        elif self.SimulType == "MechTrapping":
            for mat in self.Materials:
                counter = 1
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/D0x1", data=self.Materials[mat]["D0x1"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/D0y1", data=self.Materials[mat]["D0y1"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/DQx1", data=self.Materials[mat]["DQx1"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/DQy1", data=self.Materials[mat]["DQy1"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/Vh", data=self.Materials[mat]["Vh"])

                if self.nDim == 3:
                    self.grp_Materials.create_dataset("Material_"+str(counter)+"/Dz", data=self.Materials[mat]["Dz"])
                    
        elif self.SimulType == "GBTrapping":
            for mat in self.Materials:
                counter = 1
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/D0x1", data=self.Materials[mat]["D0x1"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/D0y1", data=self.Materials[mat]["D0y1"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/DQx1", data=self.Materials[mat]["DQx1"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/DQy1", data=self.Materials[mat]["DQy1"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/D0x2", data=self.Materials[mat]["D0x2"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/D0y2", data=self.Materials[mat]["D0y2"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/DQx2", data=self.Materials[mat]["DQx2"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/DQy2", data=self.Materials[mat]["DQy2"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/kappa_GB", data=self.Materials[mat]["kappa_GB"])

                if self.nDim == 3:
                    self.grp_Materials.create_dataset("Material_"+str(counter)+"/Dz", data=self.Materials[mat]["Dz"])
                    
        elif self.SimulType == "PhaseTrapping":
            for mat in self.Materials:
                counter = 1
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/D0x1", data=self.Materials[mat]["D0x1"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/D0y1", data=self.Materials[mat]["D0y1"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/DQx1", data=self.Materials[mat]["DQx1"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/DQy1", data=self.Materials[mat]["DQy1"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/D0x2", data=self.Materials[mat]["D0x2"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/D0y2", data=self.Materials[mat]["D0y2"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/DQx2", data=self.Materials[mat]["DQx2"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/DQy2", data=self.Materials[mat]["DQy2"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/zeta_MM", data=self.Materials[mat]["zeta_MM"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/zeta_M", data=self.Materials[mat]["zeta_M"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/zeta_fM", data=self.Materials[mat]["zeta_fM"])
                self.grp_Materials.create_dataset("Material_"+str(counter)+"/zeta_ff", data=self.Materials[mat]["zeta_ff"])

                if self.nDim == 3:
                    self.grp_Materials.create_dataset("Material_"+str(counter)+"/Dz", data=self.Materials[mat]["Dz"])
                
        #----------------------------------------------------------------------
        # Write node coordinates 
        #----------------------------------------------------------------------
            
        self.fh5.create_dataset("NodeCoordinates", data=self.nodeCoord, dtype = np.float64)
                    
        #----------------------------------------------------------------------
        # Write element node connectivity
        #----------------------------------------------------------------------
        
        # All elements
        self.fh5.create_dataset("NodeConnectivity", data=self.nodeConnectivity, dtype = np.int64)

        # Data per element set
        # element connectivity
        counter = 0
        for elSet in self.mesh.cell_sets_dict.values():
            counter+=1
            self.grp_elemSet = self.fh5.create_group('Elements/ElementSet_'+str(counter))
            elemIDs = list(elSet.values())[0]
            nElems = elemIDs.shape[0]
            print(nElems)
            self.grp_elemSet.create_dataset("nElements", data=nElems, dtype = np.int64)
            self.grp_elemSet.create_dataset("ElementSetIDs", data=elemIDs, dtype = np.int64)  
        # Nodes
        counter = 0
        for pSet in self.mesh.point_sets:
            counter+=1
            elNodes = self.mesh.point_sets[pSet]
            nElTotNodes = len(elNodes)
            print(len(self.mesh.point_sets[pSet]))
            self.fh5.create_dataset('Elements/ElementSet_'+str(counter)+'/nNodes', data=nElTotNodes, dtype = np.int64)
            self.fh5.create_dataset('Elements/ElementSet_'+str(counter)+'/NodeSetIDs', data=elNodes, dtype = np.int64)  
        
        #----------------------------------------------------------------------
        # Write prescribed dofs
        #----------------------------------------------------------------------
        
        self.grp_prescribedDOFs = self.fh5.create_group('PrescribedDOFs')
        
        for iPreDof in range(self.nPresDofs):
            self.grp_prescribedDOFs.create_dataset("Prescribed_"+str(iPreDof), data=self.presBCs[iPreDof]) 
            
        if self.SimulType in self.TransportSimulTypes:
            self.fh5.create_dataset("ExitNodes", data=self.exitNods, dtype = np.int64)

        #----------------------------------------------------------------------
        # Close hdf5 file
        #----------------------------------------------------------------------
        
        self.CloseFileHDF5()
    
        # except:
        #     # Close hdf5 file
        #     self.CloseFileHDF5()
        #     print("Some error occurred with hdf5 file handle. The file is correctly closed.")
            
        pass

#-----------------------------------------------------------------------------#

    def OpenFileHDF5(self, mode="w"):
        """
        Opens HDF5 file
        """
        
        # TODO: Do we need this with separate input/output files?
        
        # #----------------------------------------------------------------------
        # # Check if hdf5 file exists to avoid override
        # #----------------------------------------------------------------------
        
        # path = Path(self.Simul+'.hdf5')

        # if path.is_file():
        #     ErrString = 'File '+ self.Simul+ '.hdf5 exists. Can not override!'
        #     raise Exception(ErrString)
        
        self.fh5 = h5py.File(self.Simul+"_in.hdf5", mode)
        
        pass

#-----------------------------------------------------------------------------#
    
    def CloseFileHDF5(self):
        """
        Closes HDF5 file
        """
        
        self.fh5.close()
        
        pass

#-----------------------------------------------------------------------------#

    @staticmethod
    def TensileDisp2D(ly, yDisp, mesh):
        """
        Applies tensile displacement boundary conditions to a regular 
        quadrilateral in the y direction. The origin point must be (0,0)
        """
        
        #----------------------------------------------------------------------
        # Prepare data  
        #----------------------------------------------------------------------
        
        # List of prescribed degrees of freedom. Order of list [node id, dof, value]
        presBCs = []   
        # Number of nodes
        nNodes = mesh.points.shape[0]  
        # Node coordinates
        nodeCoord = mesh.points[:,0:2]
        
        #----------------------------------------------------------------------
        # Loop through nodes  
        #----------------------------------------------------------------------
    
        for iNod in range(nNodes):
            # Bottom nodes
            if nodeCoord[iNod][1]==0:
                # Find bottom left corner node.
                if nodeCoord[iNod][0] == 0: 
                    # Apply fixed BC
                    presBCs.append([iNod, 0, 0])
                    presBCs.append([iNod, 1, 0]) 
                # Other bottom nodes
                else :  
                    # Y-fixed   
                    presBCs.append([iNod, 1, 0])

            # Top nodes
            elif nodeCoord[iNod][1]==ly:
                # Top left corner node
                if nodeCoord[iNod][0] == 0:
                    # Fix-x 
                    presBCs.append([iNod, 0, 0])
                    presBCs.append([iNod, 1, yDisp])
                # Other top nodes
                else :
                    # Prescribed y
                    presBCs.append([iNod, 1, yDisp])
        
        return presBCs

#-----------------------------------------------------------------------------#

    # TODO: Merge with TensileDisp2D
    @staticmethod
    def TensileDisp3D(lz, zDisp, mesh): 
        """
        Applies tensile displacement boundary conditions to a regular hexahedron
        in the z direction. The origin point must be (0,0,0)
        """
        
        #----------------------------------------------------------------------
        # Prepare data  
        #----------------------------------------------------------------------
        
        # List of prescribed degrees of freedom. Order of list [node id, dof, value]
        presBCs = []   
        # Number of nodes
        nNodes = mesh.points.shape[0]  
        # Node coordinates
        nodeCoord = mesh.points
        
        #----------------------------------------------------------------------
        # Loop through nodes  
        #----------------------------------------------------------------------
    
        for iNod in range(nNodes):
            
            # x = 0
            if nodeCoord[iNod][0]==0:
                presBCs.append([iNod, 0, 0])
                # x = 0 & y = 0
                if nodeCoord[iNod][1]==0:
                    presBCs.append([iNod, 1, 0])
                    # x = 0 & y = 0 & z = 0
                    if nodeCoord[iNod][2]==0:
                        presBCs.append([iNod, 2, 0])
                    # x = 0 & y = 0 & z = lz
                    elif nodeCoord[iNod][2]==lz:
                        presBCs.append([iNod, 2, zDisp])
                # x = 0 & z = 0
                elif nodeCoord[iNod][2]==0:
                    presBCs.append([iNod, 2, 0])
                # x = 0 & z = lz
                elif nodeCoord[iNod][2]==lz:
                    presBCs.append([iNod, 2, zDisp])
            # y = 0
            elif nodeCoord[iNod][1]==0:
                presBCs.append([iNod, 1, 0])
                # z = 0
                if nodeCoord[iNod][2]==0:
                    presBCs.append([iNod, 2, 0])
                # z = lz
                elif nodeCoord[iNod][2]==lz:
                    presBCs.append([iNod, 2, zDisp])
                    
            # z = 0
            elif nodeCoord[iNod][2]==0:                    
                # z-fixed   
                presBCs.append([iNod, 2, 0])
            # z = lz
            elif nodeCoord[iNod][2]==lz:
                presBCs.append([iNod, 2, zDisp])
                
        return presBCs

#-----------------------------------------------------------------------------#

    @staticmethod
    def PermeationX(lx, Con_b, mesh, Charging=False):
        """
        Applies boundary conditions for permeation simulation to a regular 
        quadrilateral in the x direction. The origin point must be (0,0)
        """
        
        #----------------------------------------------------------------------
        # Prepare data  
        #----------------------------------------------------------------------
        
        # List of prescribed degrees of freedom. Order of list [node id (dof) value]
        conBCs = []   
        exitNods = []     # for flux calculation
        # Number of nodes
        nNodes = mesh.points.shape[0]  
        # Node coordinates
        nodCoord = mesh.points[:,0:2]
        
        #----------------------------------------------------------------------
        # Loop through nodes  
        #----------------------------------------------------------------------
    
        for iNod in range(nNodes):
            # Left nodes
            if nodCoord[iNod][0]==0.0:
                conBCs.append([iNod, Con_b])                    
            # Right nodes
            elif nodCoord[iNod][0] == lx:
                if Charging:
                    conBCs.append([iNod, Con_b])
                else:
                    conBCs.append([iNod, 0])
                exitNods.append([iNod])
                
        return conBCs, exitNods

#-----------------------------------------------------------------------------#

    @staticmethod
    def Uptake(lx, ly, Con_b, mesh):
        """
        Applies boundary conditions for uptake simulation to a regular 
        quadrilateral in the x direction. The origin point must be (0,0)
        """
        
        #----------------------------------------------------------------------
        # Prepare data  
        #----------------------------------------------------------------------
        
        # List of prescribed degrees of freedom. Order of list [node id (dof) value]
        conBCs = []   
        # Number of nodes
        nNodes = mesh.points.shape[0]  
        # Node coordinates
        nodCoord = mesh.points[:,0:2]
        
        #----------------------------------------------------------------------
        # Loop through nodes  
        #----------------------------------------------------------------------
    
        for iNod in range(nNodes):
            # Left nodes
            if nodCoord[iNod][0]==0.0:
                conBCs.append([iNod, Con_b])                    
            # Right nodes
            elif nodCoord[iNod][0] == lx:
                conBCs.append([iNod, Con_b])
            # Bottom nodes
            elif nodCoord[iNod][1] == 0:
                conBCs.append([iNod, Con_b])
            # Top nodes
            elif nodCoord[iNod][1] == ly:
                conBCs.append([iNod, Con_b])
                
        return conBCs
    
#-----------------------------------------------------------------------------#

    @staticmethod
    def WriteDispBCs(Simul, elementName, mesh, presBCs, dispDofs=2):
        
        #----------------------------------------------------------------------
        # Prepare data  
        #----------------------------------------------------------------------
        
        # Number of prescribed dofs
        nPresDofs = len(presBCs)
        # Number of nodes
        nNodes = mesh.points.shape[0] 
        # Node coordinates
        nodeCoord = mesh.points
        # Node connectivity
        nodeConnectivity = mesh.cells_dict[elementName] 
        
        # Vector of displacement dofs values (like solution vector)
        uDisp = np.zeros(dispDofs*nNodes)
        # Vector of displacement dofs flag 
        flagsDisp = np.zeros(dispDofs*nNodes, dtype=int)
        
        # Vector of displacement dofs values (compatible with meshio)
        sdisp = np.zeros((nNodes, dispDofs))
        # Vector of displacement dofs flag
        fdisp = np.zeros((nNodes, dispDofs), dtype=int)

        #----------------------------------------------------------------------
        # Create new mesh object for writing BCs  
        #----------------------------------------------------------------------
        
        cells = [
            (elementName, nodeConnectivity),
        ]
        BCmesh = meshio.Mesh(
            nodeCoord,
            cells,
        )
        
        #----------------------------------------------------------------------
        # Arrange BCs as solution vector  
        #----------------------------------------------------------------------

        for i in range(nPresDofs):
            pDOF = presBCs[i][0]*dispDofs + presBCs[i][1]
            uDisp[pDOF] = presBCs[i][2]
            flagsDisp[pDOF] = 1
            
        #----------------------------------------------------------------------
        # Rearrange BCs array of `nNodes x dispDofs` for meshio output
        #----------------------------------------------------------------------
            
        n = nNodes*dispDofs

        sdisp[:,0] = uDisp[0:n:dispDofs]
        sdisp[:,1] = uDisp[1:n:dispDofs]         

        fdisp[:,0] = flagsDisp[0:n:dispDofs]
        fdisp[:,1] = flagsDisp[1:n:dispDofs]
        
        if dispDofs == 3:
            sdisp[:,2] = uDisp[2:n:dispDofs]
            fdisp[:,2] = flagsDisp[2:n:dispDofs]
        
        #----------------------------------------------------------------------
        # Append and write to vtk
        #----------------------------------------------------------------------
        
        BCmesh.point_data.update({"disp": sdisp})
        BCmesh.point_data.update({"FlagBC": fdisp})
        BCmesh.write(Simul+"_BC.vtk")
        
        return BCmesh
    
#-----------------------------------------------------------------------------#

    @staticmethod
    def WriteConBCs(Simul, elementName, mesh, presBCs, dims=2):
        
        #----------------------------------------------------------------------
        # Prepare data  
        #----------------------------------------------------------------------
        
        # Number of prescribed dofs
        nPresDofs = len(presBCs)
        # Number of nodes
        nNodes = mesh.points.shape[0] 
        # Node coordinates
        nodeCoord = mesh.points
        # Node connectivity
        nodeConnectivity = mesh.cells_dict[elementName] 

        #----------------------------------------------------------------------
        # Create new mesh object for writing BCs  
        #----------------------------------------------------------------------
        
        cells = [
            (elementName, nodeConnectivity),
        ]
        BCmesh = meshio.Mesh(
            nodeCoord,
            cells,
        )
        
        #----------------------------------------------------------------------
        # Arrange BCs as solution vector  
        #----------------------------------------------------------------------

        uCon = np.zeros(nNodes)
        flagsCon = np.zeros(nNodes, dtype=int)

        for i in range(nPresDofs):
            pDOF = presBCs[i][0]
            uCon[pDOF] = presBCs[i][1]
            
            flagsCon[pDOF] = 1
        
        #----------------------------------------------------------------------
        # Append and write to vtk
        #----------------------------------------------------------------------
        
        BCmesh.point_data.update({"con": uCon})
        BCmesh.point_data.update({"FlagBC": flagsCon})
        BCmesh.write(Simul+"_BC.vtk")
        
        return BCmesh
    
#-----------------------------------------------------------------------------#