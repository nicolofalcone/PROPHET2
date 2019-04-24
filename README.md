# PROPHET2Plot
## PROPHET2UncertaintyPlot.py
PROPHET2UncertaintyPlot.py is a Python script to plot RELAP5-3D results obtained by uncertainty analysis performed with RAVEN on the PROPHET2 test facility located at Politecnico di Torino.

It has been mainly conceived to plot the numerical results obtained performing uncertainty analysis using the software RAVEN toghether with the thermal hydraulic code RELAP5-3D, both developed at INL.

It consist of a class and different functions. The functions are:

- DataPlot(variable, t_start, t_end):This function experimental data of the "variable" paramenter against the numerical results obtained for the same parameters through the different RELAP5-3D simulazion performed during the uncertainty analysis.
  
- DataPlotPercentile(variable, t_start, t_end): Similar to the previous function, but instead of plotting the numerical results of all the numerical simulation, it just compare the experimental data with the percentile 5% and 95% computed by RAVEN during the basic statistics calculation.
  
- Pearson(): This function plots the correlation coefficients computed by RAVEN during the basic statistics calculation and stored in the Pearson matrix.
  
- Cobweb(): This function normalized the parameters perturbed durign the uncertainty analysis in the range [0,1] for every numerical simulation. Then it plots the normalized values of all the parameters used in every simulation in a plot in which the x-coordinate  is the n-parameters and the y-axis is the normalized value of the parameters. The values related to the same simulation are plotted  with the same color.

## UncertaintyRange.xlsx
UncertaintyRange.xlsx is an Excel table that contains the lower and upper bounds of the range of variation for all the selected parameters for the uncertainty analysis. The data stored in this file are used in the normalization of the input values used in the numerical simulations.

## UncertaintyAnalysis_PROPHET2_test18.xlsx
UncertaintyAnalysis_PROPHET2_test18.xlsx is an Excel table thet contains the list of the parameters perturbed for the uncertainty analysis. For each parameters the following information are reported:

- the name of the parameter

- the type of distribution adopted (normal, uniform, triangular, etc.)

- the reference value i.e. the value used as input in the reference case

- the range of variation with respect to the reference value (-%/+% of the reference value)

- the upper and lower bound of the range of variation

- the card and the words modified in the RELAP5-3D input, written in the form "card:word"

- reference, if any.