# Spectral Energy Density analysis of phonon modes in 2D materials 

This code performs the Spectral Energy Density (SED) analysis of 2D materials using the trajectories from MD simulations as input. SED are calculated from the velocities in LAMMPS dump files. 

The theory behind the method implemented here is the same as the one formulated by [Thomas et. al](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.81.081411). 

This code was used to study the SED of phonon modes, their lifetime and mean free path, in Graphene/MoS2/Graphene van der Waals heterostructure. The results are published [here](https://pubs.acs.org/doi/abs/10.1021/acs.langmuir.7b03974)  

## Dependencies 

Numpy, scipy 

`dump.py` from the Pizza tools for lammps must be placed in the path 


## Description of files 

`SED.py` - the main code used to run the analysis 

`tauk.py` - this module contains methods to perform lorentzian fit to the spectral energy density peaks 

`eigen.py` - this module contains the class `eigen` with in-built methods to perform spectral energy decomposition 

`dump.py` - module that comes with [pizza](https://pizza.sandia.gov/)

`PowerSpectrum-HeatCapacity.py` - code to compute the Heat Capacity


## Restrictions 

The atom IDs in the dump file must be sorted. 

## Input: 

1. dump file containing velocity 
2. mapping of atoms into the unit cells and basis positions 
3. eigen vectors for the potential model used - can be computed using GULP 


## To do 
Generalize code for dispersion along other symmetry lines 
(currently does only along gamma-M) 

## Citation 

Please cite the following work if you use this code - *Srilok Srinivasan, Ganesh Balasubramanian, Reduced Thermal Transport in the Graphene/MoS2/Graphene Heterostructure: A Comparison with Freestanding Monolayers,Langmuir 2018, 34, 10, 3326–3335*. DOI - https://doi.org/10.1021/acs.langmuir.7b03974
