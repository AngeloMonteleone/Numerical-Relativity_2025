# Numerical-Relativity_2025
Repository for the 2024-25 course in Numerical Relativity. The repository is organized as follows
## Homework 1
There is one repository for each exercise of the Assignement (Homework 1.1, Homework 1.2, Homework 1.3).
### Homework 1.1 and 1.2
- Scripts
  - "study_advection.py" This is used to plot the evolution of the function at some given instants and the evolution of the norm.
  - "repeated_advection.py" This is used to compare the results for different courant factors. The output is a single plot with the evolution of the norm for the different $$C_f$$ values.
  - "makeani.py" Script used to animate the evolution of the function. It can be used either to plot a live animation while the integration is going or to produce a GIF with the evolution of the function over time. 
- There are also different directories containing the images used for the assignments and some example animations.
### Homework 1.3
- Scripts
  - "study_burgers.py" This is used to plot the evolution of the function at some given instants.
  - "repeated_burgers.py" This is used to compare the results for different number of points $$J$$. The output is a single plot with the evolution of the norm for the different $$J$$ values.- 
  - "makeani.py" Script used to animate the evolution of the function. It can be used either to plot a live animation while the integration is going or to produce a GIF with the evolution of the function over time. 
- There are also different directories containing the images used for the assignments and some example animations.
## Homework 2
There are two repositories, one for the first excercise (the SOD shock tube) and one for the second (the simulation of a stable TOV solution)
### SOD
The repository contains the parameter files used for the simulations ("Sod_1d_200.par","Sod_1d_400.par","Sod_1d_800.par") and some plots. The images containing the word "subplot" consist of three plots, one for the initial state, one for an intermediate state at half the simulated time and one at the final instant of the simulation (compared also with the exact solution)
- "SOD_subplots" for the evolution of the density
- "SOD_subplots_vel" for the evolution of the velocity
- "SOD_subplots_press" for the evolution of the pressure
### TOV
I chose the second option, introducing a perturbation in pressure and evolving the system with three different resolutions. For comparison I also performed just one unperturbed simulation. There are two repositories:
- "unperturbed simulation" for the files concerning the one unperturbed simulation (PAR file: "TOV_20_unperturbed.par") 
- "perturbed simulations" for the files concerning the three simulations with the pressure perturbation (PAR files: "TOV_20_perturbed.par", "TOV_15_perturbed.par", "TOV_10_perturbed.par")

The names "TOV20","TOV15","TOV10" refer, respectively, to the simulations with spacings 2.0,1.5,1.0. There are GIF's of the whole evolutions and also the plots for the evolutions maximum of the rest mass density.
