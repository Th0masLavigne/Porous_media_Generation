# Porous_media_Generation_and_properties

This repository contains the codes used for:

- the generation of [porous voronoi cubes](./Cube/reamde.md) for the training of a neural network computing meniscii permeability values, 
- the creation of [porous cylinders](./Cylinder/reamde.md) thought to be the samples for my first PhD experiment, as well as the ones for the experimental and numerical evaluation of the permeabilitty,
- the [permeabilitty evaluation](./Permeabilitty_Evaluation/reamde.md) procedure (c and python codes)

## Creation of a tesselation

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

```


## Creation of a porous medium from a tesselation