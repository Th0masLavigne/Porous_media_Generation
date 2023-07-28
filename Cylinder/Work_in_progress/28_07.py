# Extract Vertices and Edges from .tess file Created using neper.
# Support : https://neper.info/
# https://pymesh.readthedocs.io/en/latest/
# Done in meters
# 
from custom_functions_26_07 import *
from custom_stl2voxel import *
import numpy 
import subprocess
import os
import stltovoxel
# 
# 
if __name__ == "__main__"  :
		# Define the parameters of the problem:
		number_of_seeds =5
		# creating a cylinder with pymesh we specify the number of points per 2*pi: num_segments=Ncyl; and the radius of the cylinder
		Ncyl = 15
		radius = 0.5e-3 #[m]
		# creating a sphere requires in pymesh the refinement order: refinement_order=Nsph and the same radius as for the cylinder
		Nsph = 1
		# define the cylinder coordinates in m
		Diameter, Height = 24e-3, 40e-3
		# Define the parameters of the vascular tube
		thickness_tube = 0.3e-3 #[m]
		rad_i_tube = 0.5e-3
		rad_e_tube = rad_i_tube+thickness_tube
		Length_tube = 60e-3
		# Define the parameters for the connectors of the permeabilitty sample
		Diameter_connector = 1/8 #pouce externe
		Diameter_connector_mm = Diameter_connector*25.4e-3 #[m]
		Diameter_connector_i = Diameter_connector_mm / 2
		Length_connector = 60e-3 #[mm]
		chamber_height = 10e-3
		chamber_width = 0.75e-3
		# Define the outputs file names for the merged mesh, the boolean (CSG) mesh and the cleaned one
		output_filename_merged = str(number_of_seeds)+'_seed_merged.stl'
		output_filename_perm = str(number_of_seeds)+'_seed_permeabilitty.stl'
		# Define the Neper command. Several examples are provided here-after. 
		cmd = f"neper -T -n {number_of_seeds} -morpho \"diameq:dirac(1),1-sphericity:lognormal(0.145,0.03)\" -domain \"cylinder({Height},{Diameter})\""
		#
		filename = 'n'+str(number_of_seeds)+'-id1.tess'
		if os.path.isfile(filename)==False:
			print('Beginning of neper')
			subprocess.run(cmd, shell=True, check=True)
			print('Neper ended')
		# 
		# Extract vertices and edges from neper output
		vertices, edges, number_of_cells = read_tess_file(filename)
		# 
		# vertices, edges, number_of_cells =read_tess_file_mm_to_m(filename)
		#
		# Check that the Neper process has created the good output before continuing
		try:
			assert(number_of_cells == number_of_seeds)
		except:
			print("inconsistency between file and tess command")
			exit()
		# 
		# Generate a merged mesh to ensure having a geometry fulfilling the expectations
		# create_merged_mesh(vertices, edges, radius, Height, Ncyl, Nsph, output_filename_merged)
		#
		# Here is a quick interactive loop which allows you to stop the code at that point if you are unsatisfied by your geometry 
		mesh_perm = create_permeabilitty_sample_ANAS(vertices, edges, radius, Height, rad_i_tube, rad_e_tube, Length_tube, Diameter, Diameter_connector_mm, Diameter_connector_i, Length_connector,chamber_height,chamber_width, Ncyl, Nsph, output_filename_perm)
		#
		#
		voxel_size_ = radius
		resolution_order = int((Height+2*chamber_height)/voxel_size_)
		fold_name = 'n'+str(number_of_seeds)+'seed_'+str(resolution_order)+'_vx_mm'
		try: 
			os.mkdir(fold_name)
		except:
			pass
		# UL_convert_files([output_filename_perm],fold_name+'/xp', fold_name+'/chip', fold_name+'/meta', number_of_seeds, i_res=resolution_order, parallel= False)
		# UL_convert_files_cyl([output_filename_perm],fold_name+'/xp', fold_name+'/chip', fold_name+'/meta', 0.9*Diameter, number_of_seeds, i_res=resolution_order, parallel= False)
		UL_convert_files_cyl_3([output_filename_perm],fold_name+'/xp', fold_name+'/chip', fold_name+'/meta', Diameter, 0.9*Diameter, number_of_seeds, chamber_height, expected_inlet_vx=6, voxel_size__=voxel_size_, parallel=False)