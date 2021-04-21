===================
CosmoMC-relike
===================
:CosmoMC-relike: Extended CosmoMC with RELIKE (Reionization Effective Likelihood)
:Homepage: https://github.com/chenheinrich/CosmoMC-relike

Description and installation
=============================

CosmoMC-relike uses the generic sampler of CosmoMC to sample the fortran implementation of the `relike` likelihood. 

For more information on CosmoMC and getdist (the plotting package), see https://cosmologist.info/cosmomc/readme.html.

- Install/Load MPI (optional)

  It is recommended to load/install MPI for running chains. To do so
  
  - On a cluster: Find and load the MPI module (e.g. `openmpi`, `mpich` or `pmi`) on the cluster using `module avail` and `module load XX`; consult the cluster’s user guidelines).
  - On a laptop: Install OpenMPI (https://www.open-mpi.org/) using your system’s package manager (`sudo apt install libopenmpi` in Debian-based systems)

- Use `git clone --recurse-submodules` if you are cloning this repo for the first time.

  Use `git submodule update --init --recursive` if you already used `git clone` without getting the submodules. 

- Compile the code: 

  `cd CosmoMC-relike/cosmomc`
  
  `make`
  
- Untar the chain files used for KDE:

  `tar -zxvf relike_data/pl18_zmax30/chains.tar.gz -C relike_data/pl18_zmax30/`

- Run an example by outputting a single point: 

  `./cosmomc relike_example_tanh_kde_single_point.ini`

- Run an example of tanh chains in Gaussian mode: 

  `./cosmomc relike_example_tanh_gauss_chains.ini`
  
  or with MPI:
  
  `mpirun -np 4 ./cosmomc relike_example_tanh_gauss_chains.ini`
  

Using the code
==================

<more description goes here.>

Algorithm details
==================

See the latest `paper <http://arxiv.org/abs/...>`_. <to be added>

Related code
==================

The Python package `relike <https://github.com/chenheinrich/RELIKE>`_ is a python 
version of the reionization effective likelihood code used in CosmoMC-relike. It 
includes only the Gaussian approximation mode for now (it does not have the KDE mode).
You can easily incorporate the python likelihood with other samplers such as 
`Cobaya <https://github.com/CobayaSampler/cobaya>` or 
`COSMOSIS <https://bitbucket.org/joezuntz/cosmosis/wiki/Home>`. 

Branches
=============================

The master branch contains latest changes to the main release version.

The develop branch is a development branch.

=============
