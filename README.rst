===================
CosmoMC-RELIKE
===================
:CosmoMC-RELIKE: Extended CosmoMC with Reionization Effective Likelihood (RELIKE)
:Homepage: https://github.com/chenheinrich/CosmoMC-RELIKE

If you use this code, please cite `Heinrich & Hu 2021 <https://arxiv.org/abs/2104.13998>`_.

CosmoMC-RELIKE uses the generic sampler of CosmoMC to sample the fortran implementation of the RELIKE likelihood. 

For more information on CosmoMC and getdist (the plotting package), see `here <https://cosmologist.info/cosmomc/readme.html>`_. 

Installation
=============================

- Install/Load MPI (optional)

  It is recommended to load/install MPI for running chains (so you can calculate convergence statistics during the run). To do so
  
  - On a cluster: Find and load the MPI module (e.g. `openmpi`, `mpich` or `pmi`) on the cluster using `module avail` and `module load XX`; consult the cluster’s user guidelines).
  - On a laptop: Install `OpenMPI <https://www.open-mpi.org/>`_ using your system’s package manager (`sudo apt install libopenmpi` in Debian-based systems)

- Make sure all submodules are updated during cloning (skip if you already cloned it recursively as a submodule) ::

      git clone --recursive https://github.com/chenheinrich/CosmoMC-RELIKE.git 
      cd CosmoMC-RELIKE
      
- If you already cloned without using the --recurse-submodules flag, you can still update the submodules::

      cd CosmoMC-RELIKE
      git submodule update --init --recursive
  
- Compile the code

  - **With MPI** (this is the default option, which is recommended for running chains)::
  
      make
  
  - **Without MPI**::

      export BUILD=NOMPI
      make
  
- Untar the chain files used for KDE::

     tar -zxvf relike_data/pl18_zmax30/chains.tar.gz -C relike_data/pl18_zmax30/


Running Examples
=============================
- Run an example for a single tanh model in KDE mode (note that most of the time reported for this single point computation is for loading the chain, which is a one-time upfront cost that is not repeated when running chains):: 

    ./cosmomc relike_example_tanh_kde_single_point.ini

- Run an example chain for the tanh model in Gaussian mode without MPI (or with MPI). This will start a chain running, you can stop it at any time using Ctrl+C:: 

    ./cosmomc relike_example_tanh_gauss_chains.ini (mpirun -np 4 ./cosmomc relike_example_tanh_gauss_chains.ini)

- Because we only modified the generic sampler in the original CosmoMC, you should be able to run the original test as well::

    ./cosmomc test.ini

Using the code
==================

- Choose a mode:


  - The code has two modes: The KDE and the Gaussian mode. Both are tested to be sufficiently accurate, while the KDE mode captures skewness in distributions slightly better. While KDE takes longer to run, both are very fast in comparison to a sampling of the exact likelihoods. 

  - When using the KDE likelihood, we suggest using the default value of f = 0.14 to avoid over-smoothing pa-rameter posteriors while maintaining accuracy during the KDE operation. 
  

- Define your xe(z) function:


  - Open relike_xe.f90 define your model inside the function custom_xe(z) function along with additional model parameters. Do not forget to modify the header file relike_xe.h as well accordingly. Note that the tanh model example appears more complicated than it needs to be (in order to sample from a flat tau prior, we choose for xe(z) to take in as model parameter the optical depth tau instead of redshift of transition zre, so the extra code is used to calculate zre given tau.
  
  - Copy the file relike_example_tanh_gauss_chains.ini and add parameter names and priors for this model. 
  
  - Copy the paramnames/params_relike_tanh_tau.paramnames file and add the relevant parameter names and latex labels.
  
  - Note that the model parameters or priors must be arranged to explicitly satisfy fully-ionized hydrogen and singly-ionized helium for z≤6.


Algorithm details
==================

Please see the latest paper `Heinrich & Hu 2021 <https://arxiv.org/abs/2104.13998>`_ for more details.

Related code
==================

The Python package `relike <https://github.com/chenheinrich/RELIKE>`_ is a python 
version of the reionization effective likelihood code used in CosmoMC-relike. 

- It includes only the Gaussian approximation mode for now (it does not have the KDE mode).

- You can easily incorporate the python likelihood with other samplers such as `Cobaya <https://github.com/CobayaSampler/cobaya>`_ or `CosmoSIS <https://bitbucket.org/joezuntz/cosmosis/wiki/Home>`_ 

Branches
=============================

The master branch contains latest changes to the main release version.

The develop branch is a development branch.

=============
