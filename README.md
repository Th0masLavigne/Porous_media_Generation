# Porous_media_Generation_and_properties

This repository contains the codes used for:

- the generation of [porous voronoi cubes](./Cube/reamde.md) for the training of a neural network computing meniscii permeability values, 
- the creation of [porous cylinders](./Cylinder/reamde.md) thought to be the samples for my first PhD experiment, as well as the ones for the experimental and numerical evaluation of the permeabilitty,
- the [permeabilitty evaluation](./Permeabilitty_Evaluation/reamde.md) procedure (c and python codes)

## Generation of a Porous media: pipeline

The next sub-sections briefly present the overall method for the presented pipeline. At first, a domain is tesselated. Vertices and Edges are then used to create 

![image](https://github.com/Th0masLavigne/Porous_media_Generation_and_properties/blob/main/Pipeline.png)
> *Graphical representation of the route followed for the creation of the porous scaffold. Blue elements correspond to the neper environment and red elements to the pymesh environment.*


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


### Creation of a porous medium from a tesselation


Further information can be found in [Pymesh documentation](https://pymesh.readthedocs.io/en/latest/).



### Other methods
inflate 
blender truc


Further information can be found in [Pymesh documentation](https://pymesh.readthedocs.io/en/latest/).

## From a surface mesh to a particle volume