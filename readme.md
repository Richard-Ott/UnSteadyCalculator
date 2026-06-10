# Unsteady Calculator

This toolbox has two main purposes: (1) provide forward models to calculate
nuclide concentrations for unsteady erosion scenarios, and (2) provide inversion
schemes to analyze data for these scenarios. The code is currently adapted to work
for 10Be, 14C, and 26Al nuclide measurements (not all samples need all nuclides).
All models can be run for a single sample as well as a suite of data. The code 
package also includes function to calculate standard apparent erosion
rates.

## erosion scenarios
Currently 9 erosion scenarios are supported. There are two main scenario types:
step-changes in erosion and erosion spikes. In step change models, the erosion rate
changes at certain times, whereas in spike models a sudden soil loss is simulated
by removing the upper part of the production profile. For these main scenarios, there
are four sub-scenarios with varying degrees of freedom. All models allow multiple
changes in erosion through time. However, the number of parameters should be smaller
than the number of data points (your nuclide measurements). Therefore, if you want
to resolve a model with more than one change in erosion through time, choose a simpler
erosion scenario. In all scenarios, the timing of change is assumed to be uniform
across all catchments. The scenarios are:

* 'step': a (multi-)step-change in erosion rates at one or more times, where erosion
rates are allowed to vary between catchments. 
* 'samestep': a (multi-)step-change in erosion rates, where background erosion varies 
between catchments, but erosion increases/decreases by a common factor.
* 'samebackground_step': a (multi-)step-change in erosion rates, where all background
erosion rates are the same, but catchments have different erosion rate changes
* 'samebackground_samestep': a (multi-)step-change in erosion rates, where all background
erosion rates are the same and all change factors.
* 'spike': a (multi-)spike in soil loss at one or more times, where erosion
rates are allowed to vary between catchments, as well as, the soil loss heights. 
* 'samespike': a (multi-)spike soil loss, where background erosion varies 
between catchments, but the amount of soil loss is the same.
* 'samebackground_spike': a (multi-)spike soil loss, where all background
erosion rates are the same, but catchments have varying soil loss.
* 'samebackground_samespike': a (multi-)spike soil loss, where all background
erosion rates are the same and all soil loss heights.
* 'curve' add a curve that mimicks the changes in erosion and that gets scaled by 
a scaling factor. In the test data I take a pollen curve, calculate the mean tree
pollen for certain times and inversely scale this with erosion rate. The more steps
you allow in your curve the longer the computation will take.

## Production rates
Productions rates in the Unsteady Calculator are scaled using Stone (2000). 
Rock density was set to 2.65 g/cm3 and the effective attenuation length for 
spallation to 160 g/cm3 (Balco et al., 2008). Air pressure was calculated for
the median latitude, longitude, and elevation of every catchment using ERA40 
reanalysis data. Reference scaling constants for production at sea level were 
taken from ‘the online calculator formerly known as Cronus-Earth online calculators’
(Balco et al., 2008). Muon production rates were calculated using functions provided
 by Balco (2017) for model 1A, evaluated at the median elevation of every catchment
 to a depth of 25 m. The resulting  muon production profile was subsequently 
approximated with two exponential functions, which were used for forward modelling 
cosmogenic nuclide inventories. 

## Forward model
The forward model Nforward_discretized simulates cosmogenic nuclide concentrations as 
a set of exponentials that are being modified to depth cut-offs depending on the total 
amount of erosion occuring during an erosion event or time period. 
A separte forward model calculates nuclide concentrations for steady state erosion (Nforward_steadystate.m).  
 
## Soil mixing
The code can be run with and without mixing in the soil. Sol mixing is modelled as a well-mixed
layer on top of the bedrock with a constant thickness through time. For soil mixing, switch to the
'soil_mixing' branch of this repository.

## Example scripts
* 'Test_inversion': Unified test runner for synthetic scenarios. In the user input
	section of this script, select 
	(1) algorithm profile ('quick', 'balanced', 'robust'), (2) priors, and
	(3) scenario flags (true/false for each scenario).

* 'WC_inversion': Unified runner for the Western Crete dataset. The same
	algorithm and profile selection applies here, while priors and scenario flags are
	set directly in the user input section of the script.

* 'Isolines_10Be_14C': Calculates the analytical solution for step and spike
scenarios. Remember: 1 samples --> 2 equations and three unknowns (e1,e2,t), 
and hence all parameter combinations follow a line, which this script calculates.

* 'Apparent_erosion_rates': Calculates the apparent erosion rates for 10Be, 14C, and 26Al.

* 'Curve_solver_10Be_14C': this script uses the curve scenario and the provided pollen data
set to find the erosion rate history that scales with the provided curve and fits the data.

## Inversion sampler
This toolbox currently uses the MATLAB Hamiltonian Monte Carlo (HMC):

* MATLAB Hamiltonian Monte Carlo ('hmcSampler' / 'drawSamples'). This implementation
	uses unconstrained parameter transforms internally and supports profile-based tuning
	settings via 'inversion_build_config'.

## License
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see
<https://www.gnu.org/licenses/>.

How to cite:
Ott, R. F., & Scherler, D. (2026). Unsteady Calculator v1.0 (1.0). Zenodo. 
https://doi.org/10.5281/zenodo.19942751

References:

Balco, G, Stone, JO, Lifton, NA, Dunai, TJ. 2008. A complete and easily accessible means 
of calculating surface exposure ages or erosion rates from 10Be and 26Al measurements. 
Quaternary Geochronology 3(3):174–195. https://doi.org/10.1016/j.quageo.2007.12.001 

Balco, G, 2017, Production rate calculations for cosmic-ray-muon-produced 10Be and 26Al 
benchmarked against geological calibration data, Quaternary Geochronology, 39: 150-173,
https://doi.org/10.1016/j.quageo.2017.02.001.

