# drapmodules
Implementation of NOAAs DRAP model.
For an extensive explanation of DRAP model [see this link](https://www.swpc.noaa.gov/content/global-d-region-absorption-prediction-documentation)
### Introduction

DRAP stands for D Region Absortion Prediction. It is an empirical model developped by the National Oceanic and Atmospheric Administration (NOAA) to predict the global radio affectation of HF radiocommunications caused by solar flares. These solar X rays emissions cause a Sudden Ionospheric Disturbance (SID) which absorbs radio waves in the HF range.

This model was developped over the years by making radio measurements of radio absorbtions during the ocurrence of solar flares over the HF range. These radio probes were done in northern latitudes so a more extensive measurements in different latitudes are necessary.  With this measurements, the Highest Frequency Affected by 1 db was related with the solar X ray flux and described in thew next formula.
[aqui va la fórmula.]
### Implementation
With this relation the DRAP  model was implemented making a 2x2 degree  dot array of the world map. Each dot in the array is evaluted by its distance to the subsolar point during a solar flare and the X ray flux with the HAF formula. This implementation was developped in Python3.
### Running
To run the python code (drapmodulesv0-1.py), edit in the script the input values for the model (flux,time, lines 132 and 133) which are the datetime of the solar flare in UT and the solar xray flux detected in W/m^2s.
There are already example values in the code.
The output is a global map with the dot grid painted in colors showing the HAF with colors. 
![DrapMap](https://user-images.githubusercontent.com/19211938/140848105-55649d72-6621-4ac0-afd4-9740d47568fd.png)
