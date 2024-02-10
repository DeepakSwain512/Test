# Shift-T: Automated Variable-Temperature NMR Data Analysis
This repository contains the Python code that serves as the back end of the Shift-T web server (https://meieringlab.uwaterloo.ca).

NMR chemical shifts are extremely sensitive probes of molecular structure. The temperature dependences of protein chemical shifts are of interest because they offer insight into structure and dynamics. The Shift-T web server  is designed to automate the analysis of protein solution NMR 1H-15N correlation spectra collected over a range of temperatures. Services include:
 - Automatic tracking of cross peaks over temperature (ShiftTrack)
 - Calculation of amide proton and amide nitrogen chemical shift temperature coefficients (ShiftTrack)
 - Detection and analysis of curvature in the temperature dependences of amide proton chemical shifts (Curvalyzer)
 - Figure generation (both ShiftTrack and Curvalyzer)
