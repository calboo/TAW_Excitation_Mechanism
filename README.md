# TAW_Mode_Conversion

Torsional oscillations in a magnetic pore in the photosphere were detected by Marco Stangalini using data from the Interferometric Bidimensional Spectropolarimeter (IBIS) instrument at the Dunn Solar Telescope. These oscillations were interpreted as torsional Alfv&egrave;n waves (TAWs). To support this interpretation I ran simulations in Lare3d, https://warwick.ac.uk/fac/sci/physics/research/cfsa/people/tda/larexd/, that demonstrated a possible formation mechanism for m=1 TAWs in a tube like structure. 

The simulations showed that any shear Alfv&egrave;n wave interacting with a tubular structure will induced kink waves that will subsequently mode couple to generate TAWs. This interaction was shown to be ubiquitous acroos a wide range of parameters for the tube and the wave driving. In this simulation we embedded a high density tube structure within a magnetic field with exponentially diverging field lines. We used a plasma beta of one to simulate conditions in the lower chromosphere. More details of the detection and of these simulations can be found in this paper: https://www.nature.com/articles/s41550-021-01354-8.

This repository contains all the files needed to reproduce the simulations described and the figures shown in the main body and supplementary material for this paper.

## Simulation files

The folders Original and Wide each contain four files:

- control.f90
- initial_conditions.f90
- boundary.f90
- shared_data.F90

To utilise these files in Lare3d the first three should replace files of the same name in the src directory and the fourth should replace a file of the same name in the core subdirectory which is found within the src directory. Information on how to run Lare3d is given in the Lare3d manual which can be accessed from the webpage https://warwick.ac.uk/fac/sci/physics/research/cfsa/people/tda/larexd/.

### Original

These files are used to reproduce the original simulation described in the paper. This simulation was used to produce all figures and videos with the exception of figure 1c in the main body and figure 5 in the supplementary material. The simulation parameters that can be changed can all be found in the shared_data.F90 file, under TAW mode conversion simulation parameters, and are as follows:

- H - the magnetic scale height
- rampsteep - the gradient at which density increases across the tube boundary
- vortexrad - the radius of the central high density tube
- omega - the fundamental frequency for wave driving at the lower boundary
- a0 - the amplitude of wave driving at the lower boundary
- t0 - the rampup time for wave driving at the lower boundary

Note that all parameters are given in normalised units, see Lare3d manual for more info.

### Wide

These files are used to reproduce a wider variantion of the original simulation. The purpose of this simulation was to show the TAWs in the context of the wider magnetic structure. This simulation was used to produce figure 1c in the main body and figure 5 in the supplementary material. The only real difference for this simulation is that the size of the domain is doubled in the horizontal x and y directions.

## Visualisation files

The visualisation script files are all located in the Visualisation_scripts folder. These IDL scripts were used to produce the figures and videos in https://www.nature.com/articles/s41550-021-01354-8 and its supplementary material. The purpose of each script is as follows:

- density_ring.pro - plots a coloured contour of density variation around a ring of radius 2 Mm as it varies with time throughout the simulation, used to produce the third panel in Figure S6 in the supplementary material.
- fields.pro - used to confirm the equivalence of the analytic magnetic coordinates and the magnetic field used in the simulation. Also plots fields lines over a coloured contour of magnetic potential, used to produce Figure S2 in the supplementary material.
- stream_b_mag.pro - plots partial magnetic field lines over a coloured contour of the total magnetic field strength at a particular height and time, used with outputs from the wide simulation to produce Figure 1C in the main body and Figure S5 in the supplementary material.
- stream_b_rho.pro - plots partial magnetic field lines over a coloured contour of the density at a particular height and time, used with outputs from the original simulation to produce Figure S4 in the supplementary material.
- stream_diff.pro - plots velocity streamlines over a coloured contour of density contrast, relative to the initial density, at a particular height and time, used with outputs from the original simulation to produce Figure 4 in the main body and the first two panels in Figure S6 in the supplementary material.
- stream_v_mag.pro - plots velocity streamlines over a coloured contour of the total magnetic field strength at a particular height and time, used with outputs from the wide simulation to produce Figure S5 in the supplementary material.

The scripts are all written in IDL and require the Start.pro script distributed with Lare3d in order to run, in addition to the SDF and IDL directories distributed with Lare3d that allow the Lare3d output files to be read within IDL. The user must initiate idl with the command idl Start.pro and can then proceed to run the chosen visualisation script. 

The user must also be sure that the simulation data is available for the simulation and the data is being read from the right location within the visualisation script.


