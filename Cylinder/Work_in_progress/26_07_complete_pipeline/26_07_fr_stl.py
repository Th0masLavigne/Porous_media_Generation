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
		resolution_order = 5
		fold_name = 'n'+str(number_of_seeds)+'seed_'+str(resolution_order)+'_vx'
		try: 
			os.mkdir(fold_name)
		except:
			pass
		UL_convert_files([output_filename_perm],fold_name+'/xp', fold_name+'/chip', i_res=resolution_order, parallel= False)