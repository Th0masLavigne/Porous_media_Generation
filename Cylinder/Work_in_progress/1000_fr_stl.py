# Extract Vertices and Edges from .tess file Created using neper.
# Support : https://neper.info/
# https://pymesh.readthedocs.io/en/latest/
# Done in meters
# 
from custom_stl2voxel import *
import numpy
import os
import stltovoxel
# 
# 
if __name__ == "__main__"  :
		# Define the parameters of the problem:
		number_of_seeds =1000
		# Define the outputs file names for the merged mesh, the boolean (CSG) mesh and the cleaned one
		output_filename_perm = str(number_of_seeds)+'_seed_permeabilitty.stl'
		#
		#resolution_order = 7
		#UL_convert_files([output_filename_perm],fold_name+'/xp', fold_name+'/chip', fold_name+'/meta', number_of_seeds, i_res=resolution_order, parallel= False)
		Diameter = 24e-3
		chamber_height = 10e-3
		chamber_width  = 0.75e-3
		# voxel_size = 1e-3
		voxel_size = 2e-4
		res_ex = int(1e-3/voxel_size)
		fold_name = 'n'+str(number_of_seeds)+'seed_'+ str(res_ex)+'_vx_per_mm_4'
		try: 
			os.mkdir(fold_name)
		except:
			pass
		#UL_convert_files_cyl([output_filename_perm],fold_name+'/xp', fold_name+'/chip', fold_name+'/meta', Diameter, 0.9*Diameter, number_of_seeds, i_res=resolution_order, parallel= False)
		# UL_convert_files_cyl_2([output_filename_perm],fold_name+'/xp', fold_name+'/chip', fold_name+'/meta', Diameter, 0.9*Diameter, number_of_seeds, resolution_expected=res_ex, parallel=False)
		# UL_convert_files_cyl_3([output_filename_perm],fold_name+'/xp', fold_name+'/chip', 
		# 											fold_name+'/meta', Diameter, 0.9*Diameter, number_of_seeds, chamber_height, expected_inlet_vx=6, voxel_size__=1e-3, parallel=False)
		# UL_convert_files_cyl_4_debug([output_filename_perm],fold_name+'/xp', fold_name+'/chip',
		# 											fold_name+'/meta', Diameter, Diameter+2*chamber_width, number_of_seeds, chamber_height, expected_inlet_vx=6, voxel_size__=voxel_size, parallel=False)
		UL_convert_files_cyl_4([output_filename_perm],fold_name+'/xp', fold_name+'/chip',
													fold_name+'/meta', Diameter, Diameter+2*chamber_width, number_of_seeds, chamber_height, expected_inlet_vx=6, voxel_size__=voxel_size, parallel=False)