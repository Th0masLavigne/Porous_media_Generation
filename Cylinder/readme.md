# Synthetic porous cylinders
## General

This folder contains the samples which were generated for [Lavigne *et al.*](url)[^1]. From a cylindrical geometry, the tessalation is generated using neper. The porous scafold is created using Pymesh CSG tree and exported as a stl. The stl is then voxelised for EDAC (particle) simulations and permeability evaluation.

<img 
    style="display: block; 
           margin-left: auto;
           margin-right: auto;
           width: 100%;"
src=./cylinder_porous_pipeline.png
alt="pipeline">
</img>

[^1]: Lavigne, Thomas *et al.*, Titre, url


## Organisation
```
.
├── Launcher_threads_fr_stl_userRes.sh: *HPC bacth file*
├── custom_functions_cylinder.py: *Contains geometry -> stl functions*
├── custom_stl2voxel.py: *Contains stl -> voxel functions*
├── main_geo_to_voxel.py: *Main file from the geometry*
├── main_stl_to_voxel.py: *Main file from existing stls*
└── readme.md
```

### Available functions
Some unused functions are available in the custom_functions and custom_stl2voxel files:
- create_inflated_mesh(vertices, edges, ratio, output_filename): create a wireframe and inflate it. It was not retained as we had difficulty in properly controlling the inflation when the number of seeds was high.
- create_boolean_mesh_cyl_sph(vertices, edges, radius, Height, rad_i_tube, rad_e_tube, Length_tube, output_filename): create a tube in the middle of the domain (for my PhD experiment). Has a difference operator.
- fix_mesh(mesh, radius, tolerance, detail=0.1,outer_hull=True): adapted from [example](https://github.com/PyMesh/PyMesh/blob/main/scripts/fix_mesh.py) and [example](https://pymesh.readthedocs.io/en/latest/api_geometry_processing.html)
- create_boolean_mesh_cube(vertices, edges, radius, Lx, Ly, Lz, Ncyl, Nsph, output_filename): create the case for the cubic domain
- create_merged_mesh(vertices, edges, radius, Height, output_filename): cylinders and spheres are just merged. Allows to have a quick look on the expected result for the CSG tree.
- read_tess_file(filename): extract vertices and edges.
- create_permeabilitty_sample(vertices, edges, radius, Height, rad_i_tube, rad_e_tube, Length_tube, Diameter, Diameter_connector_e, Diameter_connector_i, Length_connector, Ncyl, Nsph, output_filename): **BE CAREFUL with the chamber, some values are still to be changed directly inside the function for the geometry**
- create_permeabilitty_connector(vertices, edges, radius, Height, rad_i_tube, rad_e_tube, Length_tube, Diameter, Diameter_connector_e, Diameter_connector_i, Length_connector, Ncyl, Nsph, output_filename):**BE CAREFUL with the chamber, some values are still to be changed directly inside the function for the geometry**