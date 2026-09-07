# Condensin Reverse Stepping Simulations
MATLAB code to simulate the effect of reverse stepping by condensin.

The code uses a matlab plot recoloring module [linespecer](https://www.mathworks.com/matlabcentral/fileexchange/42673-beautiful-and-distinguishable-line-colors-colormap) as a dependency.

### To fit parameters:
Open `Three_Step_Loop_Extrusion/DataFit_wDwellTimeGraph.m` and run each section. The section `Fit with 4 parameters` will perform the fit and output four predicted parameters. These parameters are then used for the rest of the code.

### To plot step-size distributions:
Open `Three_Step_Loop_Extrusion/threestepLEsimConstForce.m`, change the `currForce = ` value to the appropriate force and run each section.

