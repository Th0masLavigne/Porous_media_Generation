# Porous_media_Generation_and_properties

**To do update the help in custom functions codes**

This repository contains the codes used for:

- the generation of [porous voronoi cubes](./Cube/reamde.md) for the training of a neural network computing meniscii permeability values, 
- the creation of [porous cylinders](./Cylinder/reamde.md) thought to be the samples for my first PhD experiment, as well as the ones for the experimental and numerical evaluation of the permeabilitty,
- the [permeabilitty evaluation](./Permeabilitty_Evaluation/reamde.md) procedure (c and python codes)

## Generation of a Porous media: pipeline

The next sub-sections briefly present the overall method for the presented pipeline. At first, a domain is tesselated. Vertices and Edges are then used to create 

![image](https://github.com/Th0masLavigne/Porous_media_Generation_and_properties/blob/main/Pipeline.png)
> *Graphical representation of the route followed for the creation of the porous scaffold. Blue elements correspond to the neper environment and red elements to the pymesh environment.*

**Be careful, Pymesh requires numpy.__version__ <1.25**

### Creation of a tesselation

The creation of the tesselation is based on the use of [Neper](https://neper.info/) open-source software. See the [link](https://neper.info/) for further examples. 

```bash
========================    N   e   p   e   r    =======================
Info   : A software package for polycrystal generation and meshing.
Info   : Version 4.5.1-4
Info   : Built with: gsl|muparser|opengjk|openmp|nlopt|libscotch (full)
Info   : Running on 20 threads.
Info   : <https://neper.info>
Info   : Copyright (C) 2003-2022, and GNU GPL'd, by Romain Quey.
========================================================================
```

An example of a command is:

```bash
neper -T -n $N -morpho $Law -domain $DOMAIN
```

where `$N` is the required number of seeds to fill a specified `$DOMAIN` given a method `$LAW`. Here after are example of:
- `LAW`=
	- `"diameq:dirac(1),1-sphericity:lognormal(0.145,0.03)"` specifies the sphericity, it is a regularization of the pores' size
	- `"voronoi"` is a voronoi tesselation
- `$DOMAIN`=
	- `"cube({Lx},{Ly},{Lz})"` creates a cubic domain
	- `"cylinder({Lz},{2*Lx})"` computes a cylinder
	- boolean operations can be applied. See the [documentation](https://neper.info/).

The considered output for the next part of the procedure is the \*.tess file which contains the information about vertices and connectivity of the tesselated volume.

### Creation of a porous medium from a tesselation

The [Pymesh](https://pymesh.readthedocs.io/en/latest/) environment is used for the computation of the porous medium from the tesselation. The method consist in creating a [CSG tree](https://pymesh.readthedocs.io/en/latest/mesh_boolean.html?highlight=CSG#csg-tree) given the high number of required boolean operations needed to obtain a conform mesh.

The vertices and edges are extracted from the \*.tess file. A cylinder is created according to the connectivity table. If the egde is at the top or bottom interface, spheres were added to ensure continuity of the medium. In the provided example the radii of cylinders and spheres was kept constant but its transformation to variable values is straightforward.

Let's take the example of the Cube case:

1. Create a cylinder around each element (using the edges and vertices lists from the \*.tess file):
```python
list_union = []
for edge in edges:
	cylinder  = pymesh.generate_cylinder(vertices[edge[0]], vertices[edge[1]], radius, radius, num_segments=Ncyl)
	list_union.append({"mesh": cylinder})
for i in range(len(vertices)):
	if numpy.abs(vertices[i,2] - Lz) < 1e-8 or numpy.abs(vertices[i,2]) < 1e-8:
		sphere = pymesh.generate_icosphere(radius, vertices[i], refinement_order=Nsph)
		list_union.append({"mesh": sphere})
``` 
1. Add a box to "cut" the part of the geometry that is now outside of the initial/expected domain:
```python
box = pymesh.generate_box_mesh([0,0,0], [Lx,Ly,Lz], num_samples=1, keep_symmetry=False, subdiv_order=0, using_simplex=True)
```
1. Create the CSG tree (note that differences could be added in case of need):
```python
csg = pymesh.CSGTree({
		"intersection": [{"mesh": box}, {"union": list_union}] 
		})
```

1. Compute the boolean operation and export the result:
```python
mesh = csg.mesh
pymesh.save_mesh(output_filename, mesh, ascii=True)
```

Further information can be found in [Pymesh documentation](https://pymesh.readthedocs.io/en/latest/).



### Available functions
Other functions are available in the custom functions:
- create_inflated_mesh(vertices, edges, ratio, output_filename): create a wireframe and inflate it. It was not retained as we had difficulty in properly controlling the inflation when the number of seeds was high.
- create_boolean_mesh_cyl_sph(vertices, edges, radius, Height, rad_i_tube, rad_e_tube, Length_tube, output_filename): create a tube in the middle of the domain (for my PhD experiment). Has a difference operator.
- fix_mesh(mesh, radius, tolerance, detail=0.1,outer_hull=True): adapted from [example](https://github.com/PyMesh/PyMesh/blob/main/scripts/fix_mesh.py) and [example](https://pymesh.readthedocs.io/en/latest/api_geometry_processing.html)
- create_boolean_mesh_cube(vertices, edges, radius, Lx, Ly, Lz, Ncyl, Nsph, output_filename): create the case for the cubic domain
- create_merged_mesh(vertices, edges, radius, Height, output_filename): cylinders and spheres are just merged. Allows to have a quick look on the expected result for the CSG tree.
- read_tess_file(filename): extract vertices and edges.
- create_permeabilitty_sample(vertices, edges, radius, Height, rad_i_tube, rad_e_tube, Length_tube, Diameter, Diameter_connector_e, Diameter_connector_i, Length_connector, Ncyl, Nsph, output_filename): **BE CAREFUL with the chamber, some values are still to be changed directly inside the function for the geometry**
- create_permeabilitty_connector(vertices, edges, radius, Height, rad_i_tube, rad_e_tube, Length_tube, Diameter, Diameter_connector_e, Diameter_connector_i, Length_connector, Ncyl, Nsph, output_filename):**BE CAREFUL with the chamber, some values are still to be changed directly inside the function for the geometry**


One can try to give a look to [blender tesselation module](https://ryomizutagraphics.gumroad.com/l/TPMS_V1). A module exists in Fusion360 too but needs to be paid for.


## From a surface mesh to a particle volume

Camilo: cleaning and transforming