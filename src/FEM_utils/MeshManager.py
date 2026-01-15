"""
@file MeshManager.py
@author Abdelrahman Hussein (a.h.a.hussein@outlook.com)
@brief Class for managing Gmsh 4.1 meshes, extracting sets, and handling HDF5 separation.
@date 2026-01-14

@copyright Copyright (C) 2026 Abdelrahman Hussein

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.
"""

import meshio
import numpy as np
import h5py
from pathlib import Path
import sys

class MeshManager:
	def __init__(self, mesh_file: str, elementName: str):
		"""
		Initializes the manager and loads the mesh.
		Validates file existence and format.
		"""
  
		self.path = Path(mesh_file)
		
		# Check if file exists
		if not self.path.exists():
			print(f"[ERROR] Mesh file '{mesh_file}' not found.")
			raise FileNotFoundError(f"Could not locate: {self.path.absolute()}")
		
		# Check extension, must be .msh
		if self.path.suffix.lower() != ".msh":
			print(f"[ERROR] Invalid file format: {self.path.suffix}")
			print(f"PHIMATS strictly requires Gmsh 4.1 (.msh) files.")
			raise ValueError(f"Unsupported mesh file extension: {self.path.suffix}")

		self.elementName = elementName

		# Check for allowed elements and assign element data (number of nodes,
		# dimension and order) 
		# Format: "name": (dimension, order)
		element_registry = {
			"quad":      (2, 1),
			"quad8":     (2, 2),
			"triangle":  (2, 1),
			"triangle6": (2, 2),
			"hexahedron": (3, 1),
		}

		# Check
		if self.elementName not in element_registry:
			allowed = ", ".join(element_registry.keys())
			raise ValueError(
				f"\n[ERROR] Unknown element name <{self.elementName}>\n"
				f"Allowed elements are: {allowed}"
			)

		# Assign
		self.nDim, self.nOrder = element_registry[self.elementName]

		try:
			print(f"Reading mesh: {self.path.name}...")
			self.mesh = meshio.read(str(self.path))
   
			# Consistency Check: Ensure no mixed-topology in the domain
			# We ignore 'line' and 'vertex' because those are for boundaries/points
			domain_element_types = [
				etype for etype in self.mesh.cells_dict.keys() 
				if etype not in ["line", "vertex"]
			]

			for etype in domain_element_types:
				if etype != self.elementName:
					raise TypeError(
						f"\n[ERROR] Mixed mesh detected! The simulation is set for <{self.elementName}>, "
						f"but the mesh contains <{etype}> elements.\n"
						f"Check your Gmsh recombination settings."
					)
			
			print(f"Mesh consistency verified: All domain elements are <{self.elementName}>.")
			
			# Field data maps Name -> [Tag, Dimension]
			self.field_data = self.mesh.field_data
			# Invert to map Tag -> Name
			self.tag_to_name = {v[0]: k for k, v in self.field_data.items()}

			self.InitMaterialSets()
			self.PrintMeshSummary()
			
		except Exception as e:
			print(f"[ERROR] Failed to parse mesh file. {e}")
			sys.exit(1)

	def InitMaterialSets(self):
		"""
		Automatically detects Material Sets from Physical Groups and 
		establishes the index-based order for PHIMATS.
		"""
		detectedGroups = [
			(name, data[0]) for name, data in self.field_data.items() 
			if data[1] == self.nDim
		]
		# Sort by Tag ID
		detectedGroups.sort(key=lambda x: x[1])
		self.materialNames = [group[0] for group in detectedGroups]

		# Inform the user of the order (CRITICAL for Material_i coupling)
		print("\n" + "="*50)
		print("PHIMATS MESH INITIATION")
		print("="*50)
		print(f"Detected {len(self.materialNames)} Material Sets. Follow this order in your Input Data:")
		for i, name in enumerate(self.materialNames, start=1):
			print(f"  [Index {i}] -> Physical Group: '{name}'")
		print("="*50 + "\n")
  
	def PrintMeshSummary(self):
		"""
		Prints a concise summary of the mesh metadata to the console.
		"""
		print("-" * 30)
		print("MESH SUMMARY")
		print("-" * 30)
		print(f"File Name:      {self.path.name}")
		print(f"Element Type:   {self.elementName}")
		print(f"Dimension:      {self.nDim}D")
		print(f"Node Order:     {self.nOrder}")
		print(f"Total Nodes:    {self.get_nTotNodes()}")
		print(f"Total Elements: {self.get_nTotElements()}")
		
		# Optional: Print the specific cell blocks found by meshio
		print("Cell Blocks:")
		for cellType, data in self.mesh.cells_dict.items():
			print(f"  - {cellType}: {len(data)} items")
		print("-" * 30 + "\n")
			
	def getNodesByGroup(self, name: str):
		"""
		Extracts unique node IDs from any Physical Group (Point, Curve, or Surface).
		"""
		if name not in self.field_data:
			raise KeyError(f"'{name}' not found in Physical Groups.")
			
		targetTag, dim = self.field_data[name]
		nodes = []

		# Point (Dim 0)
		if dim == 0:
			# Point data in meshio is often in 'vertex' cells
			for cell_block, tags in zip(self.mesh.cells, self.mesh.cell_data["gmsh:physical"]):
				if cell_block.type == "vertex":
					mask = (tags == targetTag)
					nodes.extend(cell_block.data[mask].flatten())

		# Curve (Dim 1)
		elif dim == 1:
			for cell_block, tags in zip(self.mesh.cells, self.mesh.cell_data["gmsh:physical"]):
				if cell_block.type == "line":
					mask = (tags == targetTag)
					nodes.extend(cell_block.data[mask].flatten())

		# Surface/Volume (Dim 2 or 3)
		else:
			for cell_block, tags in zip(self.mesh.cells, self.mesh.cell_data["gmsh:physical"]):
				if cell_block.type == self.elementName:
					mask = (tags == targetTag)
					nodes.extend(cell_block.data[mask].flatten())

		if not nodes:
			print(f"WARNING: No nodes found for Physical Group '{name}' (Dim {dim}).")
			
		return np.unique(nodes).astype(np.int64)

	def getNodesByPoint(self, targetPoint: list, tolerance: float = 1e-6):
			"""
			Finds the ID of the node closest to a specific [x, y, z] coordinate.
			
			Args:
				targetPoint (list): [x, y] or [x, y, z] coordinates.
				tolerance (float): Maximum distance to consider a match.
			Returns:
				int: The node ID of the closest match.
			"""
			# Ensure targetPoint has 3 components for distance calculation
			point = np.zeros(3)
			point[:len(targetPoint)] = targetPoint
			
			# Calculate Euclidean distances from the target to all mesh points
			distances = np.linalg.norm(self.mesh.points - point, axis=1)
			
			# Find the index of the minimum distance
			minIdx = np.argmin(distances)
			
			if distances[minIdx] > tolerance:
				print(f"  [Warning] No node found within tolerance at {targetPoint}. "
					f"Closest is {distances[minIdx]:.2e} away.")
				return None
				
			return minIdx

	def get_nTotElements(self):
		"""
		Returns the total number of domain elements.
		Excludes boundary 'line' or 'vertex' elements.
		"""
		return len(self.mesh.cells_dict.get(self.elementName, []))

	def get_nTotNodes(self):
		"""
		Returns the total number of nodes in the mesh.
		"""
		return len(self.mesh.points)

	def get_nOrder(self):
		"""
		Returns the element order.
		"""
		return self.nOrder
	
	def get_nDim(self):
		"""
		Returns the spatial dimension (2 or 3).
		"""
		return self.nDim

	def getMaterialSet(self, name: str):
		"""
		Returns Global Element IDs and Unique Nodes for a Physical Surface/Volume.
		"""
		if name not in self.field_data:
			raise KeyError(f"Material set '{name}' not found.")

		target_tag, group_dim = self.field_data[name]

		# Check if the group dimension matches the simulation dimension
		# If simulation is 2D (quad/tri), group_dim must be 2.
		# If simulation is 3D (hex), group_dim must be 3.
		if group_dim != self.nDim:
			raise ValueError(
				f"[ERROR] Physical Group '{name}' has dimension {group_dim}, "
				f"but your simulation is set to {self.nDim}D ({self.elementName})."
			)
			
		elements = []
		nodes = []
		global_offset = 0
		
		for cell_block, tags in zip(self.mesh.cells, self.mesh.cell_data["gmsh:physical"]):
			if cell_block.type == self.elementName:
				mask = (tags == target_tag)
				if np.any(mask):
					elements.extend(np.where(mask)[0] + global_offset)
					nodes.extend(cell_block.data[mask].flatten())
				global_offset += len(cell_block.data)
				
		return {
			"elements": np.array(elements, dtype=np.int64),
			"nodes": np.unique(nodes).astype(np.int64)
		}

	def WriteMesh(self, outputName: str):
		"""
		Writes `_mesh.hdf5`.
		"""
					
		with h5py.File(outputName+"_mesh.hdf5", "w") as f:
			# Node coordinates
			f.create_dataset("NodeCoordinates", data=self.mesh.points, dtype=np.float64)
			
			# Node connectivity
			main_conn = self.mesh.cells_dict[self.elementName]
			f.create_dataset("NodeConnectivity", data=main_conn, dtype=np.int64)

			# Material/Element Sets (Elements/ElementSet_1, etc.)
			# We iterate through the material names provided in the input
			for i, mat_name in enumerate(self.materialNames, start=1):
				mat_data = self.getMaterialSet(mat_name)
				
				group_path = f"Elements/ElementSet_{i}"
				grp = f.create_group(group_path)
				
				# Write Element IDs for this set
				nElems = len(mat_data["elements"])
				grp.create_dataset("nElements", data=nElems, dtype=np.int64)
				grp.create_dataset("ElementSetIDs", data=mat_data["elements"], dtype=np.int64)
				
				# Write unique Node IDs for this set
				nNodes = len(mat_data["nodes"])
				grp.create_dataset("nNodes", data=nNodes, dtype=np.int64)
				grp.create_dataset("NodeSetIDs", data=mat_data["nodes"], dtype=np.int64)
				
				# Add the name as an attribute for debugging
				grp.attrs["PhysicalName"] = mat_name