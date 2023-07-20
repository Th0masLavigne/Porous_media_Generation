# Extract Vertices and Edges from .tess file Created using neper.
# Support : https://neper.info/
# https://pymesh.readthedocs.io/en/latest/
# 
from custom_functions19_07 import *
import numpy 
import subprocess
# 
# 
if __name__ == "__main__"  :
	seeds_list = [10]
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
	chamber_width = 0.5e-3
	# Specify the expected element size after the cleaning operation, the tolerance and the maximum of authorized execution of the cleaning function/ hull is a boolean to allow or not outer_hull computation. 
	# To be avoided in case of an envelop
	level_fixing = thickness_tube/3
	tolerance = 1e-8
	maxiter = 5
	hull = True
	hull_perm = True
	for number_of_seeds in seeds_list:
		# Define the outputs file names for the merged mesh, the boolean (CSG) mesh and the cleaned one
		output_filename_merged = str(number_of_seeds)+'_seed_merged.stl'
		output_filename_CSG = str(number_of_seeds)+'_seed_CSG.stl'
		output_filename_clean = str(number_of_seeds)+'_seed_clean.stl'
		output_filename_perm = str(number_of_seeds)+'_seed_permeabilitty.stl'
		output_filename_perm_clean = str(number_of_seeds)+'_seed_permeabilitty_clean.stl'
		# Extract vertices and edges from neper output
		filename = 'n'+str(number_of_seeds)+'-id1.tess'
		# vertices, edges, number_of_cells = read_tess_file(filename)
		vertices, edges, number_of_cells =read_tess_file_mm_to_m(filename)
		# Check that the Neper process has created the good output before continuing
		try:
			assert(number_of_cells == number_of_seeds)
		except:
			print("inconsistency between file and tess command")
			exit()
		mesh_perm, mesh_perm_top, mesh_perm_bottom = create_permeabilitty_sample_ANAS(vertices, edges, radius, Height, rad_i_tube, rad_e_tube, Length_tube, Diameter, Diameter_connector_mm, Diameter_connector_i, Length_connector,chamber_height,chamber_width, Ncyl, Nsph, output_filename_perm)
		print("Clean middle part")
		mesh_inter = mesh_perm
		manifold=False
		iter_=0
		while manifold==False:
			print(f'Cleaning, iter_: {iter_}')
			mesh_clean = clean_mesh(mesh_inter, radius, level_fixing, tolerance, hull_perm, output_filename_perm_clean)
			mesh_inter = mesh_clean
			manifold = mesh_inter.is_manifold()
			iter_+=1
			if iter_>maxiter:
				break
