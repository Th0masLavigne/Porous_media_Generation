# Porous_media_Generation

This repository contains the codes used for:

- the generation of [porous voronoi cubes](./Cube/reamde.md) for the training of a neural network computing meniscii permeability values ([Rajabi *et al.*](https://hdl.handle.net/10993/57127)[^1]), 
- the creation of [porous cylinders](./Cylinder/reamde.md) for synthetic porous media structures (3D printing) and numerical evaluation of the permeabilitty ([Lavigne *et al.*](url)[^2]).


```
.
├── ContGeoGen_V4.sif
├── Cube
│   ├── Readme.md
│   ├── cube_voronoi.png
│   ├── custom_functions_cube.py
│   └── main.py
├── Cylinder
│   ├── Launcher_threads_fr_stl_userRes.sh
│   ├── custom_functions_cylinder.py
│   ├── custom_stl2voxel.py
│   ├── cylinder_porous_pipeline.png
│   ├── main_geo_to_voxel.py
│   ├── main_stl_to_voxel.py
│   └── readme.md
├── Pipeline.png
└── README.md
```


## Generation of a Porous media: pipeline

The next sub-sections briefly present the overall method for the presented pipeline. At first, a domain is tesselated. Vertices and Edges are then used to create 

<img 
    style="display: block; 
           margin-left: auto;
           margin-right: auto;
           width: 100%;"
src=./Pipeline.png
alt="pipeline">
</img>
<em>*Graphical representation of the route followed for the creation of the porous scaffold. Blue elements correspond to the neper environment and red elements to the pymesh environment.*</em>


**Be careful, Pymesh requires numpy.__version__ <1.25. The use of a container might be useful.**

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




[^1]: RAJABI, Mohammadmahdi *et al.*, Physics-informed Dynamic Graph Convolutional Neural Network with Curriculum Learning for Pore-scale Flow Simulations, https://hdl.handle.net/10993/57127
[^2]: Lavigne, Thomas *et al.*, Titre, url
