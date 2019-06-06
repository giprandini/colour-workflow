# ColourWorkflow
This repository contains the workflow and the tools necessary in order to calculate the reflectivity and colour of metals from first principles within the independent-particle approximation (IPA).

#### How to use

- Install AiiDA (version 0.7) following the instructions in the documentation (https://aiida-core.readthedocs.io/en/v0.7.0/) and the plugins for Quantum ESPRESSO and SIMPLE.

#### Test example

We include an example to test the ColourWorkflow on elemental silver.
In order to run the example you need to:

- Run the script run_example.py within the example directory to launch the ColourWorkflow in order to calculate the IPA dielectric function of elemental silver.
- Run the script optical_constants.py within the tools directory to calculate and write to file reflectivity and colour giving as input the pk of the previously finished ColourWorkflow.

#### Content

scripts:
- run_elemental-metals.py          --> example of input script used to launch the ColourWorkflow (for elemental metals)
- run_intermetallics.py            --> example of input script used to launch the ColourWorkflow (for binary compounds of gold)
- optical_constants.py             --> post-processing AiiDA script to calculate several optical constants, such as the reflectivity and the colour coordinates, starting from the dielectric function stored in the AiiDA database from a previously run ColourWorkflow
- two-phase_alloys.py              --> additional python script to calculate the optical properties of two-phase binary alloys (Bruggeman model)
- mitsuba_spectra.py               --> additional python script to calculate 30 samples of the complex refractive index to be used as input for the Mitsuba renderer

#### Photorealistic rendering

- install and compile the Mitsuba renderer following the official documentation (https://www.mitsuba-renderer.org/) with the compilation flag SPECTRUM_SAMPLES = 30 in order to perform spectral rendering with 30 samples within the visible range
- use the script mitsuba_spectra.py to obtain the 30 spectral samples starting from the previously computed IPA dielectric function 
- perform the rendering using one of the example scenes freely available for the Mitsuba renderer (https://www.mitsuba-renderer.org/download)


#### References

G. Prandini, Predicting the reflectivity and colour of metals from first principles (Ph.D. thesis), École Polytechnique Fédérale de Lausanne, 2019

