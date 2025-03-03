# Created by Thomas Lavigne
# 22/03/2023
# 
def print_info_versions():
	print("""The version of Neper used is:
		
========================    N   e   p   e   r    =======================
Info   : A software package for polycrystal generation and meshing.
Info   : Version 4.5.1-4
Info   : Built with: gsl|muparser|opengjk|openmp|nlopt|libscotch (full)
Info   : Running on 8 threads.
Info   : <https://neper.info>
Info   : Copyright (C) 2003-2022, and GNU GPL'd, by Romain Quey.
========================================================================

The python3 environment is used. Pymesh can be run within a docker environment:

docker run -ti -v $(pwd):/home/pymesh/shared -w /home/pymesh/shared pymesh/pymesh bash
	""")
	pass
# 
def read_tess_file(filename):
	"""
	Given a filename, the number of cells is identified.
	This functions extracts from the .tess file generated by neper the vertices and edges and return them in two lists.

	Inputs: 
	- number_of_cells: integer
	Outputs:
	- vertices: np.array of all the vertices (N_vertices,3)
	- edges: np.array of all the edges (N_edges,2)
	"""
	# Import libraries
	import numpy
	import re
	# Extract from the file name the number of cells.
	regex = re.compile(r'\d+')
	regex.findall(filename)
	list_numbers_in_filename=[int(x) for x in regex.findall(filename)]
	number_of_cells = list_numbers_in_filename[0]
	# .tess reading and extraction
	# 9 first lines are general info, number_of_cells+1 next are seeds, number_of_cells+2 next are ori, 
	# then 1 line header for vertices, 1 line number of vertices, number_of_vertices lines
	# then 1 line header for edges, 1 line number of edges, number_of_edges lines
	with open(filename, 'r') as f:
		content = f.readlines()
		# extract size of the matrices from the file
		number_of_vertices = int(content[13+2*number_of_cells])
		number_of_edges = int(content[15+2*number_of_cells+number_of_vertices])
		vertices = numpy.zeros((number_of_vertices,3))
		edges = numpy.zeros((number_of_edges,2), dtype=int)
		ii=0
		for x in content[13+2*number_of_cells+1:13+2*number_of_cells+number_of_vertices+1]:
			row = x.split()
			vertices[ii]=[float(row[1]),float(row[2]),float(row[3])]
			ii+=1
		jj=0
		for x in content[15+2*number_of_cells+number_of_vertices+1:15+2*number_of_cells+number_of_vertices+number_of_edges+1]:
			row = x.split()
			edges[jj]=[int(row[1])-1,int(row[2])-1]
			jj+=1
	return vertices, edges, number_of_cells
# 
def create_merged_mesh(vertices, edges, radius, Height, output_filename):
	"""
	Given a table of vertices and the corresponding edges, this functions creates a volumic network by merging spheres at the vertices positions and cylinders along the edges.

	Inputs:
	- vertices: np.array (N_vertices,3)
	- edges: np.array (N_edges,2)
	- radius: radius of the cylinder and spheres. Homogeneous for now. Must have same unit as the mesh.
	- output_filename: string

	Outputs:
	- None, write a stl file.
	"""
	# Import libraries
	import numpy 
	import pymesh
	import time
	# initialize timer
	begin_t = time.time()
	# Initialize the mesh with the first edge
	mesh  = pymesh.generate_cylinder(vertices[edges[0,0]], vertices[edges[0,1]], radius, radius, num_segments=10)
	# Compute it for all edges and merge
	count=1
	for edge in edges[1:]:
		count+=1
		print('Opération ',count,'/',len(edges))
		cylinder  = pymesh.generate_cylinder(vertices[edge[0]], vertices[edge[1]], radius, radius, num_segments=10)
		mesh=pymesh.merge_meshes([mesh,cylinder])
	# Add the spheres at the vertices locations to avoid void.
	count = 0
	for i in range(len(vertices)):
		if numpy.abs(vertices[i,2] - Height) < 1e-8 or numpy.abs(vertices[i,2]) < 1e-8:
			count+=1
			print('number of added spheres:', count)
			# print('Opération ',i+1,'/',len(vertices))
			sphere = pymesh.generate_icosphere(radius, vertices[i], refinement_order=1)
			mesh=pymesh.merge_meshes([mesh,sphere])
	pymesh.save_mesh(output_filename, mesh, ascii=True)
	# Evaluate final time
	end_t = time.time()
	t_hours, tmin, tsec = int((end_t-begin_t)//3600), int(((end_t-begin_t)%3600)//60), int(((end_t-begin_t)%3600)%60)
	print(f"Operated in {t_hours} h {tmin} min {tsec} sec")
	return mesh
# 
def create_boolean_mesh_cube(vertices, edges, radius, Lx, Ly, Lz, Ncyl, Nsph, output_filename):
	"""
	Given a table of vertices and the corresponding edges, this functions creates a volumic network by creating boolean operations on a Constructive Solid Geometry tree (CSG)
	with spheres at the vertices positions and cylinders along the edges.
	https://pymesh.readthedocs.io/en/latest/mesh_boolean.html
	Inputs:
	- vertices: np.array (N_vertices,3)
	- edges: np.array (N_edges,2)
	- radius: radius of the cylinder and spheres. Homogeneous for now. Must have same unit as the mesh.
	- output_filename: string

	Outputs:
	- None, write a stl file.
	"""
	# Import libraries
	import numpy 
	import pymesh
	import time
	# initialize timer
	begin_t = time.time()
	# Create a box to restrict the final mesh to its volume
	box = pymesh.generate_box_mesh([0,0,0], [Lx,Ly,Lz], num_samples=1, keep_symmetry=False, subdiv_order=0, using_simplex=True)
	count=1
	# create list for boolean union
	list_union = []
	for edge in edges:
		count+=1
		print('Opération ',count,'/',len(edges))
		cylinder  = pymesh.generate_cylinder(vertices[edge[0]], vertices[edge[1]], radius, radius, num_segments=Ncyl)
		list_union.append({"mesh": cylinder})
	# Add the spheres at the vertices locations to avoid void.
	count = 0
	for i in range(len(vertices)):
		if numpy.abs(vertices[i,2] - Lz) < 1e-8 or numpy.abs(vertices[i,2]) < 1e-8:
			count+=1
			print('number of added spheres:', count)
			# print('Opération ',i+1,'/',len(vertices))
			sphere = pymesh.generate_icosphere(radius, vertices[i], refinement_order=Nsph)
			list_union.append({"mesh": sphere})
	print("Create the CSG tree")
	csg = pymesh.CSGTree({
		"intersection": [{"mesh": box}, {"union": list_union}] 
		})
	print("Begin boolean operation")
	mesh = csg.mesh
	print("End boolean operation")
	pymesh.save_mesh(output_filename, mesh, ascii=True)
	# Evaluate final time
	end_t = time.time()
	t_hours, tmin, tsec = int((end_t-begin_t)//3600), int(((end_t-begin_t)%3600)//60), int(((end_t-begin_t)%3600)%60)
	print(f"Operated in {t_hours} h {tmin} min {tsec} sec")
	return mesh
# 
def fix_mesh(mesh, radius, tolerance, detail=0.1,outer_hull=True):
    """https://github.com/PyMesh/PyMesh/blob/main/scripts/fix_mesh.py
	https://pymesh.readthedocs.io/en/latest/api_geometry_processing.html
    """
    import argparse
    import numpy as np    
    # from numpy.linalg import norm
    import pymesh
    target_len = detail
    print("Target resolution: {} mm".format(target_len))
    count = 0
    tol = tolerance   
    mesh, info = pymesh.remove_duplicated_vertices(mesh, tol)
    NN=info["num_vertex_merged"]
    print(f"{NN} duplicated vertices were removed")
    mesh, info = pymesh.remove_duplicated_faces(mesh)
    mesh, __ = pymesh.collapse_short_edges(mesh, target_len,
                                           preserve_feature=True)
    mesh, __ = pymesh.split_long_edges(mesh, target_len)
    num_vertices = mesh.num_vertices
    while True:
        mesh, info = pymesh.remove_duplicated_vertices(mesh, tol)
        mesh, info = pymesh.remove_duplicated_faces(mesh)
        mesh, __ = pymesh.collapse_short_edges(mesh, target_len,
                                               preserve_feature=True)
        mesh, __ = pymesh.remove_obtuse_triangles(mesh, 150.0, 100)
        if mesh.num_vertices == num_vertices:
            print('exit_loop clean')
            break
        num_vertices = mesh.num_vertices
        print("#v: {}".format(num_vertices))
        count += 1
        if count > 10: 
            print('Max iterations reached')
            break
    mesh = pymesh.resolve_self_intersection(mesh)
    mesh, __ = pymesh.remove_duplicated_faces(mesh)
    mesh, info = pymesh.remove_duplicated_vertices(mesh, tol)
    mesh, __ = pymesh.collapse_short_edges(mesh, target_len,
                                           preserve_feature=True)
    # In case of a mesh contained within an envelop, it is mandatory to remove the outer hull
    if outer_hull == True:
    	mesh = pymesh.compute_outer_hull(mesh)
    mesh, info = pymesh.remove_duplicated_vertices(mesh, tol)
    mesh, __ = pymesh.remove_duplicated_faces(mesh)
    mesh = pymesh.resolve_self_intersection(mesh)
    print(f" Manifold mesh: {mesh.is_manifold()}, if false, clean in meshmixer or blender for instance")
    return mesh
# 
def clean_mesh(mesh, radius, level_fixing, tolerance, hull, output_filename):
	"""
	clean the mesh using the fix_mesh fucntion. and time it

	Inputs:
	- vertices: np.array (N_vertices,3)
	- edges: np.array (N_edges,2)
	- radius: radius of the cylinder and spheres. Homogeneous for now. Must have same unit as the mesh.
	- output_filename: string

	Outputs:
	- cleaned mesh, write a stl file.
	"""
	# Import libraries
	import numpy 
	import pymesh
	import time
	# initialize timer
	begin_t = time.time()
	# processing
	mesh = fix_mesh(mesh, radius, tolerance, detail=level_fixing, outer_hull=hull)
	# Export the mesh
	pymesh.save_mesh(output_filename, mesh, ascii=True)
	# Evaluate final time
	end_t = time.time()
	t_hours, tmin, tsec = int((end_t-begin_t)//3600), int(((end_t-begin_t)%3600)//60), int(((end_t-begin_t)%3600)%60)
	print(f"Operated in {t_hours} h {tmin} min {tsec} sec")
	return mesh
# 
# 
# This function can be satisfying but I was unable to make it work properly given our characteristic lengths
def create_inflated_mesh(vertices, edges, ratio, output_filename):
	"""
	Given a table of vertices and the corresponding edges, this functions creates a volumic network by inflating the network with a thickness of ratio*length(edge)

	Inputs:
	- vertices: np.array (N_vertices,3)
	- edges: np.array (N_edges,2)
	- ratio: ratio of the length used for inflating
	- output_filename: string

	Outputs:
	- None, write a stl file.
	"""
	# Import libraries
	import numpy 
	import pymesh
	import time
	# initialize timer
	begin_t = time.time()
	# Initialize the wire_network and create a debug_file
	wire_network = pymesh.wires.WireNetwork.create_from_data(vertices, edges)
	wire_network.write_to_file("debug.wire")
	# Compute the thickness list based on the distance
	thickness = []
	for edge_ in edges:
		thickness.append(ratio*numpy.linalg.norm(vertices[edge_[0]]-vertices[edge_[1]]))
	# Define the inflator function. Allow for self-intersection to create the output anyway and see the complete log file
	inflator = pymesh.wires.Inflator(wire_network)
	# Smoothing operator
	inflator.set_refinement(2, "simple")
	# inflate
	inflator.inflate(thickness, per_vertex_thickness=False, allow_self_intersection=True)
	mesh = inflator.mesh
	# Export the mesh
	pymesh.save_mesh(output_filename, mesh, ascii=True)
	# Evaluate final time
	end_t = time.time()
	t_hours, tmin, tsec = int((end_t-begin_t)//3600), int(((end_t-begin_t)%3600)//60), int(((end_t-begin_t)%3600)%60)
	print(f"Operated in {t_hours} h {tmin} min {tsec} sec")
	return mesh