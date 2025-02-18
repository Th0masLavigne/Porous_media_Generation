# Based on the toolbox developed for stltofoxel package
# See licence file https://github.com/cpederkoff/stl-to-voxel/blob/master/LICENCE.md
# 
def calculate_mesh_limits(meshes):
	"""
	Calculates the global minimum and maximum limits (bounding box) for a set of 3D meshes.
	The function evaluates the minimum and maximum coordinates for all vertices across multiple meshes
	to define the bounding box that encapsulates all the meshes.

	This is useful when you want to know the overall extent of the meshes in space, 
	for tasks such as setting up a viewing window or determining the spatial range of a collection of meshes.

	Parameters:
	- meshes : A list of meshes.

	Returns:
	- mesh_min (np.ndarray): A numpy array of shape (3,), representing the minimum (x, y, z) coordinates that cover all the meshes.
	- mesh_max (np.ndarray): A numpy array of shape (3,), representing the maximum (x, y, z) coordinates that cover all the meshes.

	Process:
	1. The function starts by calculating the minimum and maximum coordinates for the first mesh.
	2. It then iterates over the remaining meshes, updating the minimum and maximum coordinates across all meshes.
	3. Finally, it returns the global minimum and maximum coordinates that encompass all the meshes.

	Notes:
	- The input meshes should be numpy arrays where each vertex is represented by a row with 3 values (x, y, z).
	- The function assumes that the input meshes are valid and each mesh is at least one vertex in size.
	- The function uses `np.minimum` and `np.maximum` to ensure the global minimum and maximum values are computed efficiently.

	"""
	import numpy as np
	mesh_min = meshes[0].min(axis=(0, 1))
	mesh_max = meshes[0].max(axis=(0, 1))
	for mesh in meshes[1:]:
		mesh_min = np.minimum(mesh_min, mesh.min(axis=(0, 1)))
		mesh_max = np.maximum(mesh_max, mesh.max(axis=(0, 1)))
	return mesh_min, mesh_max
# 
def calculate_resolution(meshes, i_res=5):
	"""
	Calculates the resolution for a 3D mesh grid based on the bounding box of a set of meshes.
	The function determines the resolution along each axis to match the longest dimension of the mesh
	with an adjustable resolution factor. It calculates a voxel size and returns the resolution
	in terms of the number of subdivisions along each axis.

	This function is useful for tasks such as voxelization or grid-based simulations, where you need 
	to adapt the resolution of the grid based on the physical size of the object being modeled.

	Parameters:
	- meshes (list of np.ndarray): A list of meshes, where each mesh is a numpy array with shape (N_vertices, 3),
	  representing 3D points (vertices) in space.
	- i_res (int, optional): The resolution exponent (default is 5). This controls how fine the resolution 
	  will be relative to the longest dimension of the mesh.

	Returns:
	- resolution (int): The optimal resolution for the grid, i.e., the number of subdivisions along the smallest axis.
	- idx_L (int): The index of the longest axis (0 for x, 1 for y, 2 for z).
	- vx_size (float): The size of each voxel (subdivision) in the longest dimension, given the resolution.
	- res_vec (np.ndarray): A numpy array representing the resolution along each axis, where the resolution
	  is represented as 2^i_res or 2^(i_res+1) depending on the axis.

	Process:
	1. The function first calculates the bounding box (min and max coordinates) of all the meshes using the 
	   `calculate_mesh_limits` function.
	2. It then determines the longest dimension of the bounding box and adjusts the resolution accordingly, 
	   while keeping the resolution in the other dimensions proportionally smaller.
	3. The voxel size is computed based on the longest dimension and the resolution.

	Notes:
	- The function relies on `calculate_mesh_limits` to determine the global bounding box of the meshes.
	- The `i_res` parameter is an exponent used to define the resolution level. A higher value for `i_res` means a finer resolution.
	- The function calculates a resolution vector that adapts the resolution to the mesh dimensions, with a higher resolution for the longest axis.

	"""
	import numpy as np
	mesh_min, mesh_max = calculate_mesh_limits(meshes)
	ROI = mesh_max-mesh_min
	idx_L = np.argmax(ROI)
	list_potential_resolution = []
	for ii in range(len(ROI)):
		if ii != idx_L:
			list_potential_resolution.append(int(ROI[idx_L]*2**i_res/ROI[ii]))
		else:
			list_potential_resolution.append(2**(i_res+1))
	resolution = np.min(list_potential_resolution)
	vx_size  = ROI[idx_L]/resolution
	res_vec = [2**(i_res),2**(i_res), 2**(i_res)]
	res_vec[idx_L]=2**(i_res+1)
	return int(resolution), idx_L, vx_size, np.array(res_vec)
# 
def UL_convert_files(input_file_paths, output_file_path, output_file_path_voxel, output_file_path_meta, seedn, i_res=5, parallel=False):
		"""
	Converts a set of STL files into a voxelized representation and writes the result to multiple output files, including:
	- A file containing the voxelized grid data (coordinates and values),
	- A file containing only the voxel values,
	- A metadata file containing relevant information about the conversion process.
	
	The function also generates a `.tif` image of the voxelized structure.

	Parameters:
	- input_file_paths (list of str): A list of paths to the input STL files to be converted. Each STL file represents a 3D mesh.
	- output_file_path (str): Path to the output file where the voxelized coordinates and values will be saved.
	- output_file_path_voxel (str): Path to the output file where only the voxel values will be saved.
	- output_file_path_meta (str): Path to the output file where metadata information (e.g., resolution, voxel size, etc.) will be saved.
	- seedn (int): A seed number to be included in the metadata for reproducibility or traceability.
	- i_res (int, optional): The resolution exponent used to calculate the voxel grid resolution (default is 5).
	- parallel (bool, optional): If `True`, enables parallel processing for the voxelization process, otherwise, it runs sequentially (default is `False`).

	Outputs:
	- Writes voxelized grid data to `output_file_path` and `output_file_path_voxel`.
	- Writes metadata about the voxelization process to `output_file_path_meta`.
	- Saves a `.tif` image of the voxelized grid.

	Process:
	1. Loads the input STL files using the `mesh.Mesh.from_file` method from the `stl` library.
	2. Calculates the optimal resolution for the voxel grid based on the bounding boxes of the meshes using the `calculate_resolution` function.
	3. Voxelizes the meshes using the `stltovoxel.convert.convert_meshes` function, which converts the meshes into a 3D voxel grid.
	4. Pads the resulting voxel grid to match the required dimensions and writes the voxel data and coordinates to the appropriate output files.
	5. Writes metadata information, including the seed number, voxel size, resolution, and domain dimensions to the metadata file.
	6. Saves the voxelized grid as a `.tif` image using the `tifffile.imsave` function.

	Notes:
	- The voxelization process depends on the `stltovoxel` package and requires the `stl` and `numpy` libraries.
	- The `tifffile` package is used to save the voxelized structure as a `.tif` image.
	- The output files are tab-delimited text files containing the voxel data, with each line representing a voxel at a specific coordinate with its corresponding value.
	- The function applies padding to ensure that the voxel grid has consistent dimensions based on the calculated resolution.

	"""
	import stltovoxel
	from stl import mesh
	import numpy as np
	meshes = []
	for input_file_path in input_file_paths:
		mesh_obj = mesh.Mesh.from_file(input_file_path)
		org_mesh = np.hstack((mesh_obj.v0[:, np.newaxis], mesh_obj.v1[:, np.newaxis], mesh_obj.v2[:, np.newaxis]))
		meshes.append(org_mesh)
	output = open(output_file_path, 'w')
	output_voxel = open(output_file_path_voxel, 'w')
	output_meta = open(output_file_path_meta, 'w')
	resolution, idx_L,vx_size, res_vec = calculate_resolution(meshes, i_res)
	voxels, scale, shift = stltovoxel.convert.convert_meshes(meshes, resolution-1, parallel)
	# padding
	pad = np.zeros((int(res_vec[2]),int(res_vec[1]),int(res_vec[0])),dtype=np.float32)
	offset = (np.asarray(pad.shape)-np.asarray(voxels.shape))//2
	print(offset)
	pad[offset[0]:offset[0]+voxels.shape[0],offset[1]:offset[1]+voxels.shape[1],offset[2]:offset[2]+voxels.shape[2]]=voxels
	for z in range(pad.shape[0]):
		for y in range(pad.shape[1]):
			for x in range(pad.shape[2]):
					point = (np.array([x, y, z]) / scale) # + shift
					output.write('%s\t%s\t%s\t' % tuple(point))
					output_voxel.write('%s\t' % pad[z][y][x])
	output_meta.write('seeding number: ')
	output_meta.write('%s\n' % seedn)
	output_meta.write('Voxel size = L_max/resolution : ')
	output_meta.write('%s\n' % vx_size)
	output_meta.write('Resolution = ')
	output_meta.write('%s x %s x %s \n' % tuple(res_vec))
	output_meta.write('domain dims = ')
	output_meta.write('%s x %s x %s' % tuple(res_vec*vx_size))
	output.close()
	output_voxel.close()
	output_meta.close()
	from tifffile import imsave
	imsave('voxelized.tif', pad)
# 
def UL_convert_files_cyl(input_file_paths, output_file_path, output_file_path_voxel, output_file_path_meta, Diameter, internal_Diameter, seedn, i_res=5, parallel=False):
	"""
	Converts a set of STL files into a voxelized representation with cylindrical regions, writes the result to multiple output files, including:
	- A file containing the voxelized coordinates and values,
	- A file containing only the voxel values,
	- A metadata file containing relevant information about the conversion process.
	
	The function also generates a `.tif` image of the voxelized structure, where regions within the external and internal cylinder diameters are processed specifically to ensure they are marked as solid.

	Parameters:
	- input_file_paths (list of str): A list of paths to the input STL files to be converted. Each STL file represents a 3D mesh.
	- output_file_path (str): Path to the output file where the voxelized coordinates and values will be saved.
	- output_file_path_voxel (str): Path to the output file where only the voxel values will be saved.
	- output_file_path_meta (str): Path to the output file where metadata information (e.g., resolution, voxel size, etc.) will be saved.
	- Diameter (float): The diameter of the external cylindrical region.
	- internal_Diameter (float): The internal diameter of the cylindrical region, used to differentiate solid and empty regions.
	- seedn (int): A seed number to be included in the metadata for reproducibility or traceability.
	- i_res (int, optional): The resolution exponent used to calculate the voxel grid resolution (default is 5).
	- parallel (bool, optional): If `True`, enables parallel processing for the voxelization process, otherwise, it runs sequentially (default is `False`).

	Outputs:
	- Writes voxelized grid data to `output_file_path` and `output_file_path_voxel`.
	- Writes metadata about the voxelization process to `output_file_path_meta`.
	- Saves a `.tif` image of the voxelized structure.

	Process:
	1. Loads the input STL files using the `mesh.Mesh.from_file` method from the `stl` library.
	2. Calculates the optimal resolution for the voxel grid based on the bounding boxes of the meshes using the `calculate_resolution` function.
	3. Voxelizes the meshes using the `stltovoxel.convert.convert_meshes` function, which converts the meshes into a 3D voxel grid.
	4. Pads the resulting voxel grid to match the required dimensions and handles special cases for the cylindrical region.
	5. Marks the cylindrical regions as solid (1) or empty (0) based on their relative position to the internal diameter.
	6. Writes voxel coordinates and values to the appropriate output files.
	7. Writes metadata information, including the seed number, voxel size, resolution, and domain dimensions to the metadata file.
	8. Saves the voxelized grid as a `.tif` image using the `tifffile.imsave` function.

	Notes:
	- The voxelization process depends on the `stltovoxel` package and requires the `stl` and `numpy` libraries.
	- The `tifffile` package is used to save the voxelized structure as a `.tif` image.
	- The output files are tab-delimited text files containing the voxel data, with each line representing a voxel at a specific coordinate with its corresponding value.
	- The function applies padding to ensure that the voxel grid has consistent dimensions based on the calculated resolution.
	- The cylindrical regions are processed to ensure that voxels within the external diameter and outside the internal diameter are marked as solid to avoid any flow through the structure.

	"""
	import stltovoxel
	from stl import mesh
	import numpy as np
	meshes = []
	for input_file_path in input_file_paths:
		mesh_obj = mesh.Mesh.from_file(input_file_path)
		org_mesh = np.hstack((mesh_obj.v0[:, np.newaxis], mesh_obj.v1[:, np.newaxis], mesh_obj.v2[:, np.newaxis]))
		meshes.append(org_mesh)
	output = open(output_file_path, 'w')
	output_voxel = open(output_file_path_voxel, 'w')
	output_meta = open(output_file_path_meta, 'w')
	resolution, idx_L,vx_size, res_vec = calculate_resolution(meshes, i_res)
	voxels, scale, shift = stltovoxel.convert.convert_meshes(meshes, resolution=resolution-1, parallel=parallel)
	# padding
	pad = np.zeros((int(2**(i_res+1)),int(2**i_res),int(2**i_res)),dtype=np.float32)
	offset = (np.asarray(pad.shape)-np.asarray(voxels.shape))//2
	rest_offset = (np.asarray(pad.shape)-np.asarray(voxels.shape))%2
	pad[offset[0]:offset[0]+voxels.shape[0],offset[1]:offset[1]+voxels.shape[1],offset[2]:offset[2]+voxels.shape[2]]=voxels
	# stl to voxel adds an empty slice at the beginning and two at the end of the stack 
	pad[0,:,:]=pad[1,:,:]
	pad[pad.shape[0]-1,:,:]=pad[pad.shape[0]-3,:,:]
	pad[pad.shape[0]-2,:,:]=pad[pad.shape[0]-3,:,:]
	for z in range(pad.shape[0]):
		for y in range(pad.shape[1]):
			for x in range(pad.shape[2]):
					point = (np.array([x, y, z]) / scale) # + shift
					# tag external of cylinder as ssolid to avoid flow
					center_cyl = np.array([Diameter/2, Diameter/2])
					xp_m_offset = np.array([(x+ shift[0]-offset[2])/scale[0], (y+ shift[1]-offset[1])/scale[1]]) # + shift
					if np.linalg.norm(xp_m_offset-center_cyl)-internal_Diameter/2 >= 0:
						pad[z][y][x]=1
					output.write('%s\t%s\t%s\t' % tuple(point))
					output_voxel.write('%s\t' % pad[z][y][x])
	output_meta.write('seeding number: ')
	output_meta.write('%s\n' % seedn)
	output_meta.write('Voxel size = L_max/resolution : ')
	output_meta.write('%s\n' % vx_size)
	output_meta.write('Resolution = ')
	output_meta.write('%s x %s x %s \n' % tuple(res_vec))
	output_meta.write('domain dims = ')
	output_meta.write('%s x %s x %s' % tuple(res_vec*vx_size))
	output.close()
	output_voxel.close()
	output_meta.close()
	from tifffile import imsave
	imsave('voxelized_cyl.tif', pad)
# 
def calculate_resolution_2(meshes, resolution_expected):
    """
    Calculates the optimal resolution for a voxel grid based on the mesh bounding box dimensions
    and a set of expected resolutions for each axis. The function computes the voxel size and
    ensures that the resolution aligns with the expected values while preserving the aspect ratio
    defined by the mesh bounding box.

    Parameters:
    - meshes (list of np.ndarray): A list of meshes, where each mesh is represented as an array of vertices.
      Each mesh is typically a set of 3D vertices defining a surface.
    - resolution_expected (list or np.ndarray): A list or array of expected resolution values along each axis 
      (x, y, z). This defines the target resolution that should be achieved for each of the axes, guiding the
      calculation of the overall voxel grid resolution.

    Returns:
    - resolution (int): The final resolution (size of the voxel grid along the smallest axis), calculated
      based on the bounding box dimensions and expected resolutions.
    - idx_L (int): The index of the largest dimension (axis) used in the resolution calculation.
    - vx_size (float): The size of a voxel along the largest dimension (axis).
    - res_vec (np.ndarray): A 3-element array containing the resolution for each axis (x, y, z) as specified by
      the `resolution_expected` input.

    Process:
    1. The function first calculates the bounding box dimensions (ROI) for each mesh using the `calculate_mesh_limits` function.
    2. It then uses the expected resolution values (`resolution_expected`) for each axis to determine the final voxel grid resolution.
    3. The resolution is computed based on the relative sizes of the mesh bounding boxes and the expected resolution along each axis.
    4. The final resolution is chosen as the minimum resolution across all axes, and the voxel size for the largest dimension is calculated.
    5. The function returns the final resolution, the index of the largest axis, the voxel size, and the expected resolutions for each axis.

    Notes:
    - The expected resolution along each axis should be provided in the correct order corresponding to the mesh's dimensions.
    - The largest axis (`idx_L`) is used as a reference to adjust the resolution in the other dimensions while maintaining the aspect ratio.
    - This function ensures that the resolution calculation is consistent with the provided expected resolutions for each axis, making it more flexible when dealing with non-cubic meshes.

    """
	# Expectation, the idx_L is already known when asking the resolution
	import numpy as np
	mesh_min, mesh_max = calculate_mesh_limits(meshes)
	ROI = mesh_max-mesh_min
	idx_L = np.argmax(resolution_expected)
	list_potential_resolution = []
	for ii in range(len(ROI)):
		if ii != idx_L:
			list_potential_resolution.append(int(ROI[idx_L]*resolution_expected[ii]/ROI[ii]))
		else:
			list_potential_resolution.append(resolution_expected[idx_L])
	resolution = np.min(list_potential_resolution)
	vx_size  = ROI[idx_L]/resolution
	res_vec = resolution_expected
	return int(resolution), idx_L, vx_size, np.array(res_vec)
# 
def UL_convert_files_cyl_2(input_file_paths, output_file_path, output_file_path_voxel, output_file_path_meta, Diameter, internal_Diameter, seedn, resolution_expected=[64,64,128], parallel=False):
    """
    Converts a set of 3D STL files into a voxelized representation and saves the results to multiple output files.
    This version of the function specifically handles cylindrical geometry, tagging regions inside a cylinder as solid to 
    avoid flow and appropriately padding the voxel grid to align with the expected resolution.

    Parameters:
    - input_file_paths (list of str): A list of paths to the input STL files that will be converted to voxels.
    - output_file_path (str): The path where the 3D coordinates (x, y, z) will be written in a tab-delimited format.
    - output_file_path_voxel (str): The path where the voxel data (binary values) will be written in a tab-delimited format.
    - output_file_path_meta (str): The path where meta information (e.g., seed number, voxel size) will be written.
    - Diameter (float): The external diameter of the cylindrical region.
    - internal_Diameter (float): The internal diameter of the cylindrical region. Used to identify the solid region inside the cylinder.
    - seedn (int): A seed number that will be recorded in the metadata to track the voxelization process.
    - resolution_expected (list or np.ndarray, optional): Expected resolution for the voxel grid along each axis (x, y, z).
      Default is `[64, 64, 128]`.
    - parallel (bool, optional): If set to `True`, the voxel conversion process will be performed in parallel. Default is `False`.

    Returns:
    - None: The function writes output files for the voxelized representation and metadata but does not return a value.

    Outputs:
    - Voxelized 3D grid in the form of a `.tif` image file containing the voxel data.
    - Tab-delimited text files containing the 3D coordinates and corresponding voxel values.
    - A metadata file with details about the voxelization process, including the seed number, voxel size, and resolution.

    Process:
    1. Reads and processes the input STL files into a list of meshes.
    2. Uses the `calculate_resolution_2` function to determine the appropriate resolution based on the expected voxel grid size and the mesh's bounding box.
    3. Converts the meshes into a voxelized format, applying necessary padding to ensure alignment with the expected resolution.
    4. Tags the external cylindrical region as solid to prevent flow through the region, ensuring that the internal area with the specified diameter is maintained.
    5. Writes the voxelized data (including coordinates and voxel values) to the specified output files.
    6. Saves a `.tif` image file of the voxelized structure using the `tifffile.imsave` function.

    Notes:
    - The cylindrical regions are tagged as solid based on their position relative to the internal diameter, ensuring that the outer cylindrical region is correctly represented.
    - Padding is applied to ensure that the voxel grid matches the specified resolution and that the voxelized structure is centered within the grid.
    - The function can process multiple STL files simultaneously and will write the results to separate output files.
    - The expected resolution should be provided in the order corresponding to the x, y, and z axes.

    """
	import stltovoxel
	from stl import mesh
	import numpy as np
	meshes = []
	for input_file_path in input_file_paths:
		mesh_obj = mesh.Mesh.from_file(input_file_path)
		org_mesh = np.hstack((mesh_obj.v0[:, np.newaxis], mesh_obj.v1[:, np.newaxis], mesh_obj.v2[:, np.newaxis]))
		meshes.append(org_mesh)
	output = open(output_file_path, 'w')
	output_voxel = open(output_file_path_voxel, 'w')
	output_meta = open(output_file_path_meta, 'w')
	resolution, idx_L,vx_size, res_vec = calculate_resolution_2(meshes, resolution_expected)
	voxels, scale, shift = stltovoxel.convert.convert_meshes(meshes, resolution=resolution-1, parallel=parallel)
	# padding
	pad = np.zeros((int(res_vec[2]),int(res_vec[1]),int(res_vec[0])),dtype=np.float32)
	offset = (np.asarray(pad.shape)-np.asarray(voxels.shape))//2
	rest_offset = (np.asarray(pad.shape)-np.asarray(voxels.shape))%2
	pad[offset[0]:offset[0]+voxels.shape[0],offset[1]:offset[1]+voxels.shape[1],offset[2]:offset[2]+voxels.shape[2]]=voxels
	# stl to voxel adds an empty slice at the beginning and two at the end of the stack 
	if offset[0]==0:
		pad[0,:,:]=pad[1,:,:]
		pad[pad.shape[0]-1,:,:]=pad[pad.shape[0]-3,:,:]
		pad[pad.shape[0]-2,:,:]=pad[pad.shape[0]-3,:,:]
	elif rest_offset[0]==0:
		for jj in range(offset[0]+1):
			pad[jj,:,:]=pad[offset[0]+1,:,:]
		for jj in range(offset[0]+voxels.shape[0]-2,pad.shape[0]):
			pad[jj,:,:]=pad[offset[0]+voxels.shape[0]-3,:,:]
	elif rest_offset[0]!=0:
		for jj in range(offset[0]):
			pad[jj,:,:]=pad[offset[0]+1,:,:]
		for jj in range(offset[0]+voxels.shape[0]-2,pad.shape[0]):
			pad[jj,:,:]=pad[offset[0]+voxels.shape[0]-3,:,:]
	for z in range(pad.shape[0]):
		for y in range(pad.shape[1]):
			for x in range(pad.shape[2]):
					point = (np.array([x, y, z]) / scale) # + shift
					# tag external of cylinder as ssolid to avoid flow
					center_cyl = np.array([Diameter/2, Diameter/2])
					xp_m_offset = np.array([(x+ shift[0]-offset[2])/scale[0], (y+ shift[1]-offset[1])/scale[1]]) # + shift
					if np.linalg.norm(xp_m_offset-center_cyl)-internal_Diameter/2 >= 0:
						pad[z][y][x]=1
					output.write('%s\t%s\t%s\t' % tuple(point))
					output_voxel.write('%s\t' % pad[z][y][x])
	output_meta.write('seeding number: ')
	output_meta.write('%s\n' % seedn)
	output_meta.write('Voxel size = L_max/resolution : ')
	output_meta.write('%s\n' % vx_size)
	output_meta.write('Resolution = ')
	output_meta.write('%s x %s x %s \n' % tuple(res_vec))
	output_meta.write('domain dims = ')
	output_meta.write('%s x %s x %s' % tuple(res_vec*vx_size))
	output.close()
	output_voxel.close()
	output_meta.close()
	from tifffile import imsave
	imsave('voxelized_cyl_2.tif', pad)
# 
def cut_before_inlet(voxels, voxel_size, chamber_height, expected_inlet_vx):
    """
    Cuts the voxelized representation of a structure to exclude the portion before the inlet,
    based on the given criteria for identifying the inlet region. This function trims the voxel grid
    to start at the first slice that represents the inlet, taking into account the expected size 
    of the inlet region in terms of the number of voxels.

    Parameters:
    - voxels (np.ndarray): A 3D numpy array representing the voxelized structure. The shape of the array
      is (Z, Y, X), where Z is the number of slices, Y is the number of rows, and X is the number of columns.
    - voxel_size (float): The physical size of each voxel. This parameter is not used in this specific implementation 
      but could be useful if further adjustments based on the physical dimensions are needed.
    - chamber_height (float): The height of the chamber. This is not directly used in the function but could be 
      useful for determining the overall size or constraints of the voxelized structure.
    - expected_inlet_vx (int): The number of voxel slices expected to represent the inlet region. This value determines 
      how many slices will be included after cutting the voxel grid.

    Returns:
    - voxels_cut (np.ndarray): A 3D numpy array representing the portion of the voxelized structure 
      starting from the inlet region and extending through the expected number of inlet voxel slices.
    - res_vec (list): A list representing the shape of the resulting voxel grid after cutting, in the order 
      [X, Y, Z] (columns, rows, slices).

    Process:
    1. Identifies the first slice (Z-index) that represents the inlet by comparing voxel values.
    2. Cuts the voxel array to retain only the portion starting from the inlet region, extending 
       `expected_inlet_vx` slices beyond it.
    3. Returns the cut voxel array and its new dimensions.

    Notes:
    - The function assumes that the inlet is located at the first slice where the sum of the voxels in the 
      slice increases, indicating the start of the porous region.
    - The function does not account for potential cases where the inlet is not found based on the given criteria.
    - The returned `voxels_cut` contains only the region starting from the detected inlet, with a fixed number 
      of expected inlet voxels included in the output.

    """
    # Expectation, the idx_L is already known when asking the resolution
	import numpy as np
	for ii in range(2,voxels.shape[0]-1):
		if np.sum(voxels[ii,:,:])<np.sum(voxels[ii+1,:,:]):
			begin_porous = ii
			break
	voxels_cut = voxels[begin_porous-(expected_inlet_vx-1):,:,:]
	res_vec = [voxels_cut.shape[2],voxels_cut.shape[1],voxels_cut.shape[0]]
	return voxels_cut, res_vec
# 
def UL_convert_files_cyl_3(input_file_paths, output_file_path, output_file_path_voxel, output_file_path_meta, Diameter, internal_Diameter, seedn, chamber_height, expected_inlet_vx=6, voxel_size__=1e-3, parallel=False):
    """
    This function processes a set of STL files representing 3D models, converts them into a voxel grid, and 
    outputs the voxel data, metadata, and voxelized file in TIFF format. It specifically targets cylindrical 
    geometries, tagging external parts of the model as solid to avoid flow. The voxel grid is padded and adjusted 
    based on the given inlet parameters and resolution.

    Parameters:
    - input_file_paths (list of str): List of file paths to the STL files that need to be processed.
    - output_file_path (str): File path for the output file where the voxel coordinates will be saved.
    - output_file_path_voxel (str): File path for the output file where the voxel values will be saved.
    - output_file_path_meta (str): File path for the output file where the metadata will be saved.
    - Diameter (float): The outer diameter of the cylindrical model. This value is used to define the outer limits 
      for the flow region.
    - internal_Diameter (float): The internal diameter of the cylindrical model. This value helps define the 
      boundary for the flow within the model.
    - seedn (int): A number used for seeding or tracking purposes. It is written to the metadata file.
    - chamber_height (float): The height of the chamber represented by the voxel grid.
    - expected_inlet_vx (int, optional): The number of voxel slices expected to represent the inlet region 
      (default is 6). This value affects the voxel slicing and grid adjustment.
    - voxel_size__ (float, optional): The size of each voxel in meters (default is 1e-3 meters).
    - parallel (bool, optional): Flag indicating whether parallel processing should be used when converting 
      the STL files to voxels (default is False).

    Returns:
    - None: The function writes the resulting voxel grid, metadata, and voxelized image to the specified output files.
    
    Process:
    1. Loads and processes the STL files to obtain the mesh data.
    2. Calculates the mesh's bounding box and resolution based on the specified voxel size.
    3. Converts the mesh data into a voxelized grid using the `stltovoxel` library.
    4. Applies padding to ensure the voxel grid has the correct dimensions.
    5. Tags the external parts of the cylinder as solid based on the given internal and external diameters.
    6. Writes the resulting voxel grid, coordinates, and metadata to the output files.
    7. Saves the voxelized model as a TIFF image file.

    Notes:
    - The voxel grid is generated based on the calculated resolution derived from the bounding box of the model.
    - Padding is applied to the voxel grid to ensure the correct dimensions are maintained during the conversion.
    - The cylinder's internal and external diameters are used to mark the solid region of the voxel grid, which is 
      then used for flow simulations or further analysis.
    - The function assumes the STL models represent cylindrical geometries and adjusts the voxel grid accordingly.

    """
    import stltovoxel
	from stl import mesh
	import numpy as np
	meshes = []
	for input_file_path in input_file_paths:
		mesh_obj = mesh.Mesh.from_file(input_file_path)
		org_mesh = np.hstack((mesh_obj.v0[:, np.newaxis], mesh_obj.v1[:, np.newaxis], mesh_obj.v2[:, np.newaxis]))
		meshes.append(org_mesh)
	# 
	output = open(output_file_path, 'w')
	output_voxel = open(output_file_path_voxel, 'w')
	output_meta = open(output_file_path_meta, 'w')
	# 
	mesh_min, mesh_max = calculate_mesh_limits(meshes)
	ROI = mesh_max-mesh_min
	idx_L = np.argmax(ROI)
	# 
	resolution = int(round(ROI[idx_L]/voxel_size__))
	voxels, scale, shift = stltovoxel.convert.convert_meshes(meshes, resolution=resolution-1, parallel=parallel)
	# padding
	voxels_cut[0,:,:]=voxels_cut[1,:,:]
	voxels_cut[voxels_cut.shape[0]-2:,:,:]=voxels_cut[voxels_cut.shape[0]-3,:,:]
	for z in range(voxels_cut.shape[0]):
		for y in range(voxels_cut.shape[1]):
			for x in range(voxels_cut.shape[2]):
					point = (np.array([x, y, z]) / scale) # + shift
					# tag external of cylinder as ssolid to avoid flow
					center_cyl = np.array([Diameter/2, Diameter/2])
					if np.linalg.norm([x/scale[0]+shift[0],y/scale[1]+shift[1]]-center_cyl)-internal_Diameter/2 >= 0:
						voxels_cut[z][y][x]=1
					output.write('%s\t%s\t%s\t' % tuple(point))
					output_voxel.write('%s\t' % voxels_cut[z][y][x])
	# 
	output_meta.write('seeding number: ')
	output_meta.write('%s\n' % seedn)
	output_meta.write('Voxel size = vx/m : ')
	output_meta.write('%s\n' % np.mean(scale))
	output_meta.write('Resolution = ')
	output_meta.write('%s x %s x %s \n' % tuple(res_vec))
	output_meta.write('domain dims = ')
	output_meta.write('%s x %s x %s' % tuple(res_vec/scale))
	output.close()
	output_voxel.close()
	output_meta.close()
	from tifffile import imsave
	tif_matrix = np.asarray(voxels_cut,dtype=np.float32)
	imsave('voxelized_cyl_3.tif', tif_matrix)
# 
def enforce_size_x_y(voxels):
    """
    This function enforces a square shape for the x and y dimensions of the input voxel grid by truncating the 
    larger dimension to match the smaller one. The function modifies the voxel grid by cutting it along the x 
    and y axes to ensure the resulting shape is square in those dimensions.

    Parameters:
    - voxels (ndarray): A 3D numpy array representing the voxel grid with dimensions (z, x, y). The function 
      modifies the x and y dimensions to make them equal by trimming the larger dimension.

    Returns:
    - voxels_cut (ndarray): The resulting 3D numpy array with the x and y dimensions cut to the same size.
    - res_vec (list): A list containing the dimensions (x, y, z) of the cut voxel grid, where x and y are equal.

    Process:
    1. The function compares the x and y dimensions of the input voxel grid.
    2. It truncates the larger dimension to match the smaller one, ensuring a square shape for the x and y axes.
    3. It returns the modified voxel grid and a list of its new dimensions.

    Notes:
    - This function does not affect the z-dimension, which remains unchanged.
    - The cut voxel grid will have the same depth (z-dimension) as the original, but the x and y dimensions will be 
      made equal by cutting the larger dimension.
    """
	# Expectation, the idx_L is already known when asking the resolution
	voxels_cut = voxels[:,:min(voxels.shape[0],voxels.shape[1]),:min(voxels.shape[0],voxels.shape[1])]
	res_vec = [voxels_cut.shape[2],voxels_cut.shape[1],voxels_cut.shape[0]]
	return voxels_cut, res_vec
# 
def UL_convert_files_cyl_4(input_file_paths, output_file_path, output_file_path_voxel, output_file_path_meta, Diameter, external_Diameter, seedn, chamber_height, expected_inlet_vx=6, voxel_size__=1e-3, parallel=False):
    """
    This function processes input STL files representing a cylindrical geometry, converts them into a voxel grid, 
    and exports the voxelized data along with metadata to specified output files. The function handles external 
    and internal cylinder boundaries, adjusts the size of the resulting voxel grid, and includes optional 
    parallel processing.

    Parameters:
    - input_file_paths (list): A list of file paths to the input STL files that define the geometry.
    - output_file_path (str): The file path where the 3D coordinates of the voxelized geometry will be saved.
    - output_file_path_voxel (str): The file path where the voxel values will be saved.
    - output_file_path_meta (str): The file path where metadata about the voxelization will be saved.
    - Diameter (float): The diameter of the internal region of the cylindrical geometry.
    - external_Diameter (float): The outer diameter of the cylindrical geometry, used for marking the external region.
    - seedn (int): The seed number to be included in the metadata.
    - chamber_height (float): The height of the chamber, used to determine the voxel region to retain.
    - expected_inlet_vx (int, optional): The expected number of inlet voxels. Defaults to 6.
    - voxel_size__ (float, optional): The voxel size in meters. Defaults to 1e-3 meters.
    - parallel (bool, optional): Whether to use parallel processing. Defaults to False.

    Returns:
    - None: The function writes voxel data and metadata to the specified output files.

    Process:
    1. Converts input STL meshes into voxelized representation using the specified resolution.
    2. Pads the voxel grid and cuts the portion before the inlet.
    3. Marks the external region of the cylinder as solid (value 3) to prevent fluid flow.
    4. Enforces a square shape for the x and y dimensions of the voxel grid.
    5. Saves the 3D coordinates and voxel values to the specified output files.
    6. Writes metadata, including voxel size, resolution, and domain dimensions.
    7. Optionally, parallel processing can be enabled for faster conversion of meshes.

    Notes:
    - The function assumes that the input STL meshes represent cylindrical geometries.
    - The final voxel grid is adjusted to ensure proper representation of internal and external cylinder regions.
    - The `voxels_cut` array is adjusted to match the chamber height and expected inlet region.

    Side Effects:
    - The function writes voxelized data and metadata to disk in the specified output file paths.
    """
    import stltovoxel
	from stl import mesh
	import numpy as np
	meshes = []
	for input_file_path in input_file_paths:
		mesh_obj = mesh.Mesh.from_file(input_file_path)
		org_mesh = np.hstack((mesh_obj.v0[:, np.newaxis], mesh_obj.v1[:, np.newaxis], mesh_obj.v2[:, np.newaxis]))
		meshes.append(org_mesh)
	# 
	output = open(output_file_path, 'w')
	output_voxel = open(output_file_path_voxel, 'w')
	output_meta = open(output_file_path_meta, 'w')
	# 
	mesh_min, mesh_max = calculate_mesh_limits(meshes)
	ROI = mesh_max-mesh_min
	idx_L = np.argmax(ROI)
	# 
	resolution = int(round(ROI[idx_L]/voxel_size__))
	voxels, scale, shift = stltovoxel.convert.convert_meshes(meshes, resolution=resolution-1, parallel=parallel)
	voxels[0,:,:]=voxels[1,:,:]
	voxels[voxels.shape[0]-2:,:,:]=voxels[voxels.shape[0]-3,:,:]
	# padding
	voxels_cut, res_vec = cut_before_inlet(voxels, scale, chamber_height, expected_inlet_vx)
	in_x, in_y, in_z, out_x, out_y, out_z = [], [], [], [], [], []
	for z in range(voxels_cut.shape[0]):
		for y in range(voxels_cut.shape[1]):
			for x in range(voxels_cut.shape[2]):
					point = (np.array([x, y, z]) / scale) # + shift
					# tag external of cylinder as ssolid to avoid flow
					center_cyl = np.array([Diameter/2, Diameter/2])
					if not (np.linalg.norm([x/scale[0]+ shift[0],y/scale[1]+ shift[1]]-center_cyl)-(external_Diameter)/2 >= voxel_size__/2 and voxels_cut[z][y][x]==0):
						output.write('%s\t%s\t%s\t' % tuple(point))
						output_voxel.write('%s\t' % voxels_cut[z][y][x])
					else:
						voxels_cut[z][y][x]=3
	# 
	voxels_cut, res_vec = enforce_size_x_y(voxels_cut)
	output_meta.write('seeding number: ')
	output_meta.write('%s\n' % seedn)
	output_meta.write('Voxel size = vx/m : ')
	output_meta.write('%s\n' % np.mean(scale))
	output_meta.write('Resolution = ')
	output_meta.write('%s x %s x %s \n' % tuple(res_vec))
	output_meta.write('domain dims = ')
	output_meta.write('%s x %s x %s' % tuple(res_vec/scale))
	output.close()
	output_voxel.close()
	output_meta.close()
	from tifffile import imsave
	tif_matrix = np.asarray(voxels_cut,dtype=np.float32)
	imsave('voxelized_cyl_4.tif', tif_matrix)
# 
def UL_convert_files_cyl_4_debug(input_file_paths, output_file_path, output_file_path_voxel, output_file_path_meta, Diameter, external_Diameter, seedn, chamber_height, expected_inlet_vx=6, voxel_size__=1e-3, parallel=False):
    """
    This debug version of the UL_convert_files_cyl_4 function processes input STL files representing a cylindrical geometry,
    converts them into a voxel grid, and exports the voxelized data along with metadata to specified output files.
    It also visualizes the voxel grid and the separation of internal and external regions for debugging purposes.

    Parameters:
    - input_file_paths (list): A list of file paths to the input STL files that define the geometry.
    - output_file_path (str): The file path where the 3D coordinates of the voxelized geometry will be saved.
    - output_file_path_voxel (str): The file path where the voxel values will be saved.
    - output_file_path_meta (str): The file path where metadata about the voxelization will be saved.
    - Diameter (float): The diameter of the internal region of the cylindrical geometry.
    - external_Diameter (float): The outer diameter of the cylindrical geometry, used for marking the external region.
    - seedn (int): The seed number to be included in the metadata.
    - chamber_height (float): The height of the chamber, used to determine the voxel region to retain.
    - expected_inlet_vx (int, optional): The expected number of inlet voxels. Defaults to 6.
    - voxel_size__ (float, optional): The voxel size in meters. Defaults to 1e-3 meters.
    - parallel (bool, optional): Whether to use parallel processing. Defaults to False.

    Returns:
    - None: The function writes voxel data and metadata to the specified output files and displays a 3D plot.

    Process:
    1. Converts input STL meshes into voxelized representation using the specified resolution.
    2. Pads the voxel grid and cuts the portion before the inlet.
    3. Marks the external region of the cylinder as solid (value 3) to prevent fluid flow.
    4. Visualizes the internal and external regions of the cylinder in a 3D scatter plot.
    5. Saves the voxelized data and metadata.

    Debugging:
    - The function visualizes the distribution of "in" (internal) and "out" (external) points in a 3D plot using Matplotlib.
    - The scatter plot will show "in" points in red and "out" points in black, helping to verify the accuracy of the voxelization process.
    - The function terminates immediately after displaying the plot for debugging purposes.

    Notes:
    - This is a debug version of the `UL_convert_files_cyl_4` function. It includes additional visualization steps for debugging.
    - The function assumes that the input STL meshes represent cylindrical geometries.
    - The final voxel grid is adjusted to ensure proper representation of internal and external features.
    """
    import stltovoxel
	from stl import mesh
	import numpy as np
	meshes = []
	for input_file_path in input_file_paths:
		mesh_obj = mesh.Mesh.from_file(input_file_path)
		org_mesh = np.hstack((mesh_obj.v0[:, np.newaxis], mesh_obj.v1[:, np.newaxis], mesh_obj.v2[:, np.newaxis]))
		meshes.append(org_mesh)
	# 
	output = open(output_file_path, 'w')
	output_voxel = open(output_file_path_voxel, 'w')
	output_meta = open(output_file_path_meta, 'w')
	# 
	mesh_min, mesh_max = calculate_mesh_limits(meshes)
	ROI = mesh_max-mesh_min
	idx_L = np.argmax(ROI)
	# 
	resolution = int(round(ROI[idx_L]/voxel_size__))
	voxels, scale, shift = stltovoxel.convert.convert_meshes(meshes, resolution=resolution-1, parallel=parallel)
	voxels[0,:,:]=voxels[1,:,:]
	voxels[voxels.shape[0]-2:,:,:]=voxels[voxels.shape[0]-3,:,:]
	# padding
	voxels_cut, res_vec = cut_before_inlet(voxels, scale, chamber_height, expected_inlet_vx)
	in_x, in_y, in_z, out_x, out_y, out_z = [], [], [], [], [], []
	for z in range(voxels_cut.shape[0]):
		for y in range(voxels_cut.shape[1]):
			for x in range(voxels_cut.shape[2]):
					point = (np.array([x, y, z]) / scale) # + shift
					# tag external of cylinder as ssolid to avoid flow
					center_cyl = np.array([Diameter/2, Diameter/2])
					if not (np.linalg.norm([x/scale[0]+ shift[0],y/scale[1]+ shift[1]]-center_cyl)-(external_Diameter)/2 >= 0 and voxels_cut[z][y][x]==0):
						if voxels_cut[z][y][x]!=0:
							in_x.append(x / scale[0])
							in_y.append(y/ scale[1])
							in_z.append(z / scale[2])
						else:
							out_x.append(x / scale[0])
							out_y.append(y/ scale[1])
							out_z.append(z / scale[2])
						output.write('%s\t%s\t%s\t' % tuple(point))
						output_voxel.write('%s\t' % voxels_cut[z][y][x])
	# 
	voxels_cut, res_vec = enforce_size_x_y(voxels_cut)
	import matplotlib.pyplot as plt
	import numpy as np
	fig = plt.figure()
	ax = fig.add_subplot(projection='3d')
	ax.set_aspect('equal', 'box')
	ax.scatter(in_x, in_y, in_z, color='red')
	ax.scatter(out_x, out_y, out_z,color='black')
	plt.show()
	exit()
# 
def UL_convert_files_cyl_5(input_file_paths, output_file_path, output_file_path_voxel, output_file_path_meta, Diameter, external_Diameter, seedn, chamber_height, expected_inlet_vx=6, voxel_size__=1e-3, parallel=False):
    """
    This function processes input STL files representing a cylindrical geometry, converts them into a voxel grid, and exports
    the voxelized data along with metadata to specified output files. It also tags the external regions of the cylinder as solid
    to prevent flow and saves the voxelized data as a TIFF image.

    Parameters:
    - input_file_paths (list): A list of file paths to the input STL files that define the geometry.
    - output_file_path (str): The file path where the 3D coordinates of the voxelized geometry will be saved.
    - output_file_path_voxel (str): The file path where the voxel values will be saved.
    - output_file_path_meta (str): The file path where metadata about the voxelization will be saved.
    - Diameter (float): The diameter of the internal region of the cylindrical geometry.
    - external_Diameter (float): The outer diameter of the cylindrical geometry, used for marking the external region.
    - seedn (int): The seed number to be included in the metadata.
    - chamber_height (float): The height of the chamber, used to determine the voxel region to retain.
    - expected_inlet_vx (int, optional): The expected number of inlet voxels. Defaults to 6.
    - voxel_size__ (float, optional): The voxel size in meters. Defaults to 1e-3 meters.
    - parallel (bool, optional): Whether to use parallel processing. Defaults to False.

    Returns:
    - None: The function writes voxel data and metadata to the specified output files and saves the voxelized data as a TIFF image.

    Process:
    1. Converts input STL meshes into voxelized representation using the specified resolution.
    2. Pads the voxel grid and cuts the portion before the inlet.
    3. Tags the external region of the cylinder as solid (value 3) to prevent fluid flow.
    4. Saves the voxelized data and metadata in the specified output files.
    5. Saves the voxelized data as a TIFF image for further processing or visualization.

    Notes:
    - The function assumes that the input STL meshes represent cylindrical geometries.
    - The function writes out both 3D coordinates and voxel values to the output files.
    - The function ensures that the x and y dimensions are squared for consistency in the voxel grid.
    - The voxelized data is saved in a TIFF file format using the `tifffile` library for potential use in image processing or visualization.

    Example:
    ```python
    UL_convert_files_cyl_5(
        input_file_paths=['path/to/mesh1.stl', 'path/to/mesh2.stl'],
        output_file_path='output/coordinates.txt',
        output_file_path_voxel='output/voxels.txt',
        output_file_path_meta='output/meta.txt',
        Diameter=1.0,
        external_Diameter=2.0,
        seedn=42,
        chamber_height=10.0,
        expected_inlet_vx=6,
        voxel_size__=1e-3,
        parallel=True
    )
    ```
    """
    import stltovoxel
	from stl import mesh
	import numpy as np
	meshes = []
	for input_file_path in input_file_paths:
		mesh_obj = mesh.Mesh.from_file(input_file_path)
		org_mesh = np.hstack((mesh_obj.v0[:, np.newaxis], mesh_obj.v1[:, np.newaxis], mesh_obj.v2[:, np.newaxis]))
		meshes.append(org_mesh)
	# 
	output = open(output_file_path, 'w')
	output_voxel = open(output_file_path_voxel, 'w')
	output_meta = open(output_file_path_meta, 'w')
	# 
	mesh_min, mesh_max = calculate_mesh_limits(meshes)
	ROI = mesh_max-mesh_min
	idx_L = np.argmax(ROI)
	# 
	resolution = int(round(ROI[idx_L]/voxel_size__))
	voxels, scale, shift = stltovoxel.convert.convert_meshes(meshes, resolution=resolution-1, parallel=parallel)
	voxels[0,:,:]=voxels[1,:,:]
	voxels[voxels.shape[0]-2:,:,:]=voxels[voxels.shape[0]-3,:,:]
	# padding
	voxels_cut, res_vec = cut_before_inlet(voxels, scale, chamber_height, expected_inlet_vx)
	in_x, in_y, in_z, out_x, out_y, out_z = [], [], [], [], [], []
	for z in range(voxels_cut.shape[0]):
		for y in range(voxels_cut.shape[1]):
			for x in range(voxels_cut.shape[2]):
					point = (np.array([x, y, z]) / scale) # + shift
					# tag external of cylinder as ssolid to avoid flow
					center_cyl = np.array([Diameter/2, Diameter/2])
					if (np.linalg.norm([x/scale[0]+ shift[0],y/scale[1]+ shift[1]]-center_cyl)-(external_Diameter)/2 >= 0 and voxels_cut[z][y][x]==0):
						voxels_cut[z][y][x]=3
					output.write('%s\t%s\t%s\t' % tuple(point))
					output_voxel.write('%s\t' % voxels_cut[z][y][x])
	# 
	voxels_cut, res_vec = enforce_size_x_y(voxels_cut)
	output_meta.write('seeding number: ')
	output_meta.write('%s\n' % seedn)
	output_meta.write('Voxel size = vx/m : ')
	output_meta.write('%s\n' % np.mean(scale))
	output_meta.write('Resolution = ')
	output_meta.write('%s x %s x %s \n' % tuple(res_vec))
	output_meta.write('domain dims = ')
	output_meta.write('%s x %s x %s' % tuple(res_vec/scale))
	output.close()
	output_voxel.close()
	output_meta.close()
	from tifffile import imsave
	tif_matrix = np.asarray(voxels_cut,dtype=np.float32)
	imsave('voxelized_cyl_5.tif', tif_matrix)
# EoF