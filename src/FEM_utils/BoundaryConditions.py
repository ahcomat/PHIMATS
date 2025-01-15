"""
@file BoundaryConditions.py
@brief A module containing some functions for assigning boundary conditions for 'PIHMATS'.
"""

import numpy as np
import meshio
from pathlib import Path

def TensileDisp2D(ly, yDisp, mesh):
    """Applies tensile displacement boundary conditions to a regular 
    quadrilateral in the y direction. The origin point must be (0,0)

    Args:
        ly (float): y-length
        yDisp (float): y-displacement
        mesh (meshio): Mesh object
    Returns:
        list: List containing [[nod id, dof, value]]
    """
    
    # List of prescribed degrees of freedom. Order of list [node id, dof, value]
    presBCs = []   
    # Number of nodes
    nNodes = mesh.points.shape[0]  
    # Node coordinates
    nodeCoord = mesh.points[:,0:2]
    
    # Loop through nodes  
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

def TensileDisp3D(lz, zDisp, mesh): 
    """Applies tensile displacement boundary conditions to a regular 
    quadrilateral in the z direction. The origin point must be (0,0,)

    Args:
        lz (float): z-length
        zDisp (float): z-displacement
        mesh (meshio): Mesh object
    Returns:
        list: List containing [[nod id, dof, value]]
    """
       
    # List of prescribed degrees of freedom. Order of list [node id, dof, value]
    presBCs = []   
    # Number of nodes
    nNodes = mesh.points.shape[0]  
    # Node coordinates
    nodeCoord = mesh.points
    
    # Loop through nodes  
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


def PermeationX(lx, Con_b, mesh, Charging=False):
    """Applies boundary conditions for charging/discharging BCs to a regular 
    quadrilateral in the x direction. The origin point must be (0,0)
    
    Args:
        lx (float): Edge length 
        Con_b (float): Boundary concentration
        mesh (meshio): Mesh object
        Charging (bool, optional): Flag for charging/discharging. Defaults to False.

    Returns:
        list: List containing [[nod id, value]]
        list: Exit nodes ids
    """
    
    # List of prescribed degrees of freedom. Order of list [node id (dof) value]
    conBCs = []   
    exitNods = []     # for flux calculation
    # Number of nodes
    nNodes = mesh.points.shape[0]  
    # Node coordinates
    nodCoord = mesh.points[:,0:2]
    
    # Loop through nodes  
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

def WriteDispBCs(Simul, elementName, mesh, presBCs, dispDofs=2):
    """Function to help visualize the displacement BCs in Paraview. Writes mesh with prescibed BCs as vtk.

    Args:
        Simul (str): Simulation name
        elementName (str): meshio compatible element name.
        mesh (meshio): Mesh object.
        presBCs (list): List containing [[nod id, dof, value]].
        dispDofs (int, optional): Displacement dofs. Defaults to 2.

    Returns:
        meshio: Mesh object with BCs.
    """
       
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

    # Create new mesh object for writing BCs  
    cells = [
        (elementName, nodeConnectivity),
    ]
    BCmesh = meshio.Mesh(
        nodeCoord,
        cells,
    )
    
    # Arrange BCs as solution vector  
    for i in range(nPresDofs):
        pDOF = presBCs[i][0]*dispDofs + presBCs[i][1]
        uDisp[pDOF] = presBCs[i][2]
        flagsDisp[pDOF] = 1
        
    # Rearrange BCs array of `nNodes x dispDofs` for meshio output
    n = nNodes*dispDofs

    sdisp[:,0] = uDisp[0:n:dispDofs]
    sdisp[:,1] = uDisp[1:n:dispDofs]         

    fdisp[:,0] = flagsDisp[0:n:dispDofs]
    fdisp[:,1] = flagsDisp[1:n:dispDofs]
    
    if dispDofs == 3:
        sdisp[:,2] = uDisp[2:n:dispDofs]
        fdisp[:,2] = flagsDisp[2:n:dispDofs]
    
    
    # Append and write to vtk
    BCmesh.point_data.update({"disp": sdisp})
    BCmesh.point_data.update({"FlagBC": fdisp})
    BCmesh.write(Simul+"_BC.vtk")
    
    return BCmesh

#-----------------------------------------------------------------------------#

@staticmethod
def WriteConBCs(Simul, elementName, mesh, presBCs, dims=2):
    """Function to help visualize the concentration BCs in Paraview. Writes mesh with prescibed BCs as vtk.

    Args:
        Simul (str): Simulation name
        elementName (str): meshio compatible element name.
        mesh (meshio): Mesh object.
        presBCs (list): List containing [[nod id, value]].
        dispDofs (int, optional): Displacement dofs. Defaults to 2.

    Returns:
        meshio: Mesh object with BCs.
    """
    
    # Number of prescribed dofs
    nPresDofs = len(presBCs)
    # Number of nodes
    nNodes = mesh.points.shape[0] 
    # Node coordinates
    nodeCoord = mesh.points
    # Node connectivity
    nodeConnectivity = mesh.cells_dict[elementName] 

    # Create new mesh object for writing BCs      
    cells = [
        (elementName, nodeConnectivity),
    ]
    BCmesh = meshio.Mesh(
        nodeCoord,
        cells,
    )
    
    # Arrange BCs as solution vector  
    uCon = np.zeros(nNodes)
    flagsCon = np.zeros(nNodes, dtype=int)

    for i in range(nPresDofs):
        pDOF = presBCs[i][0]
        uCon[pDOF] = presBCs[i][1]
        
        flagsCon[pDOF] = 1
    
    # Append and write to vtk    
    BCmesh.point_data.update({"con": uCon})
    BCmesh.point_data.update({"FlagBC": flagsCon})
    BCmesh.write(Simul+"_BC.vtk")
    
    return BCmesh