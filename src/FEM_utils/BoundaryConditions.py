"""
@file BoundaryConditions.py
@brief A module containing some functions for assigning boundary conditions for 'PIHMATS'.
"""

import numpy as np
import meshio
from pathlib import Path

def AssignDirichletBC(nodeIDs, dof, value):
    """
    Creates a list of BCs for a given set of nodes.
    
    Args:
        nodeIDs (ndarray): Node indices from MeshManager.
        dof (int): DOF index (0:x/scalar, 1:y, 2:z)
        value (float): Prescribed value.
    """
            
    return [[nodeID, dof, value] for nodeID in nodeIDs]

def WriteBCVTK(Simul, meshManager, presBCs, dofNames=["ux", "uy"]):
    """
    Visualizes boundary conditions in Paraview.
    
    Args:
        Simul (str): Base name for the file.
        meshManager (MeshManager): Mesh manager object.
        presBCs (list): The combined list of all [node, dof, value].
        dofNames (list): Name(s) for the DOFs (e.g. ["ux", "uy", "phi", "con"]).
    """
    nNodes = meshManager.get_nTotNodes()
    nDofs = len(dofNames)
    
    # Initialize data arrays
    bcValues = np.zeros((nNodes, nDofs))
    bcFlags = np.zeros((nNodes, nDofs), dtype=int)

    # Fill arrays from the presBCs list
    for entry in presBCs:
        nodeID, dof, val = entry
        if dof < nDofs:
            bcValues[nodeID, dof] = val
            bcFlags[nodeID, dof] = 1

    # Create mesh for VTK
    cells = [(meshManager.elementName, meshManager.mesh.cells_dict[meshManager.elementName])]
    outMesh = meshio.Mesh(meshManager.mesh.points, cells)

    # Add data to points
    for i, name in enumerate(dofNames):
        outMesh.point_data[f"val_{name}"] = bcValues[:, i]
        outMesh.point_data[f"flag_{name}"] = bcFlags[:, i]

    outMesh.write(f"{Simul}_BC.vtu")
    print(f"BC visualization written to {Simul}_BC.vtu")