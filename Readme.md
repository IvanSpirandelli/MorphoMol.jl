# MorphoMol.jl

The main purpose of the code in this package is to simulate molecules in a solvent utilizing the [*morphometric approach to solvation free energy*](https://pubmed.ncbi.nlm.nih.gov/36638318/). To this end it implements versions of *Random Walk Metropolis* as well as *Hamiltonian Monte Carlo*, which can be found under 'src/modules/Algorithms/src/'. 

Furthermore, calculations of the energy can be found in 'src/modules/Energies/src/' and some Utilities for simulation setup and evaluation under 'src/simulation_setup.jl'.

Under 'examples' you will find simple to complex examples utilizing both algorithms implemented in this package. 

There is an an accompanying repository: [MorphoMolHPC](https://github.com/IvanSpirandelli/MorphoMolHPC) containing code to setup simulations on a High Performance Cluster.

If you have any questions don't hesitate to reach out to spirandelli@uni-potsdam.de
