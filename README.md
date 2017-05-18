# README #
These scripts reconstruct canonical flux tubes and calculate their helicity from experimental measurements of
magnetic field, density, electron temperature and Mach numbers of a plasma.
I used these scripts to analyse data from the Reconnection Scaling Experiment (RSX) for my PhD dissertation.

Canonical flux tubes trace out cross-sections of constant flux of circulation of a species’ canonical momentum.
In the zero-flow limit, canonical flux tubes are magnetic flux tubes, but in full form, 
present the distinct advantage of reconciling all plasma regimes (single particle, kinetic, two-fluid, and MHD) 
with topological concepts of helicity, twists and linkages.
The theory of canonical helicity transport predicts the threshold for the observed bifurcation of end states in 
compact toroid merging experiments and hypothesizes a mechanism for converting destabilizing magnetic twist 
into stabilizing helical shear flows in astrophysical jets.


### Canonical flux tubes and their helicity in RSX  ###
Ion and electron canonical flux tubes are visualized from a volumetric dataset of Mach, triple, and $\dot{B}$ probe 
measurements at over $10,000$ spatial locations of a gyrating kinked plasma column in a sub-volume of the RSX 
experiment. The flux tubes gyrate into and out of the sub-volume.
The flux tubes twist with the ion canonical flux tube writhing around the electron canonical flux tube.
Relative canonical helicity is calculated by solving for reference fields with discrete cosine transforms and 
choosing gauges so that vector potentials can be calculated from simple integrals.
The cross helicity between the magnetic and flow vorticity flux tubes dominates the ion canonical helicity.
The magnetic and cross helicity are anti-correlated.


### Dependencies ###
* python 2.7.12
* numpy 1.11.2
* scipy 0.18.1
* matplotlib 1.5.3
* seaborn 0.7.1
* sqlite 3.13
* vtk 6.3.0
* MDSplus 7.7.1 and mdsplus python module (alpha7.0.144)
* VisIt 2.10.3 and corresponding visit python module


The dependencies can be installed with the anaconda python distribution.
The MDSplus can be found at http://www.mdsplus.org/index.php/Introduction.
The mdsplus python module can be found at https://pypi.python.org/pypi/MDSplus.
VisIt can be found at https://visit.llnl.gov.
Instructions for installing the VisIt https://visit.llnl.gov/manuals

### RSX data and analysis scripts ###
This code processes both raw data and processed data from the RSX experiment.
The RSX data and RSX data analysis codes have not been released yet. 
Without these inputs the code can not currently be run, however, it could be adapted to run on other datasets.
When they are a note will be added to the Github repository and Zenodo on how to obtain them. 

### Directory structure ###
1. Create two directories `output`, `source`. 
2. Git clone or copy the repo into `source`.

### Important files ###
Each of these files has documentation that can be accessed with `python filename --help`

`write_measurements_to_unstructured_grid.py` reads processed Bdot and triple probe, and raw Mach probe measurements 
and writes them to unstructured vtk files files.

`fit_field_null.py` reads unstructured Bx and By fields and finds the field nulls in a measurement plane.

`interpolate_measurements.py` interpolates fields from unstructured measurements to rectilinear grid.

`calculate_dependent_quantities.py` calculates the dependent quantities needed to plot flux tubes and 
calculate helicity from interpolated measurements.

`calculate_helicity.py` calculates gauge-dependent and relative helicities.

`vis_canonical_flux_tubes.py` plots frames of canonical flux tube animations.
Plot options include ion, electron canonical flux tubes, temperature and density isosurfaces, current contours
in an x-y plane.

`thesis_figures/thesis_figures.ipynb` a jupyter notebook generating the figures for the 5th chapter of my thesis.

### Reference ###
von der Linden. 2017. Investigating the Dynamics of Canonical Flux Tubes. (PhD thesis) (In preparation.) 

### Contact ###
Jens von der Linden jensv@uw.edu

