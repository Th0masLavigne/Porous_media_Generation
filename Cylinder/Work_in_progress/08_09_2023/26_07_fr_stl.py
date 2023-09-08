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
		number_of_seeds =5
		# Define the outputs file names for the merged mesh, the boolean (CSG) mesh and the cleaned one
		output_filename_perm = str(number_of_seeds)+'_seed_permeabilitty.stl'
		#
		#resolution_order = 7
		res_ex = [80,80,160]
		fold_name = 'n'+str(number_of_seeds)+'seed_'+ str(res_ex[0])+'x'+ str(res_ex[1])+'x'+ str(res_ex[2])+'_vx_2'
		try: 
			os.mkdir(fold_name)
		except:
			pass
		#UL_convert_files([output_filename_perm],fold_name+'/xp', fold_name+'/chip', fold_name+'/meta', number_of_seeds, i_res=resolution_order, parallel= False)
		Diameter = 24e-3
		#UL_convert_files_cyl([output_filename_perm],fold_name+'/xp', fold_name+'/chip', fold_name+'/meta', Diameter, 0.9*Diameter, number_of_seeds, i_res=resolution_order, parallel= False)
		UL_convert_files_cyl_2([output_filename_perm],fold_name+'/xp', fold_name+'/chip', fold_name+'/meta', Diameter, 0.9*Diameter, number_of_seeds, resolution_expected=res_ex, parallel=False)