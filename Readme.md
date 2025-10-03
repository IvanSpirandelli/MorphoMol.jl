# MorphoMol.jl

The main purpose of the code in this package is to simulate molecules in a solvent utilizing the [*morphometric approach to solvation free energy*](https://pubmed.ncbi.nlm.nih.gov/36638318/) and a combined approach using it in tandem with a *topological biasing potential*. To this end it implements versions of *Random Walk Metropolis*, *Simulated Annealing* and *Hamiltonian Monte Carlo*, which can be found under 'src/modules/Algorithms/src/'. 

Publications using this code: 
* [Topological potentials guiding protein self-assembly](https://arxiv.org/abs/2508.15321)
* [Solvation, geometry, and assembly of the tobacco mosaic virus](https://academic.oup.com/pnasnexus/article/4/3/pgaf065/8042116)
* [Exotic self-assembly of hard spheres in a morphometric solvent](https://www.pnas.org/doi/10.1073/pnas.2314959121)

Calculations of the energy can be found in 'src/modules/Energies/src/' and some Utilities for simulation setup and evaluation under 'src/simulation_setup.jl'. 

There is an an accompanying repository: [MorphoMolHPC](https://github.com/IvanSpirandelli/MorphoMolHPC) containing code to set up simulations on a High Performance Cluster.

If you have any questions don't hesitate to reach out to spirandelli@uni-potsdam.de
