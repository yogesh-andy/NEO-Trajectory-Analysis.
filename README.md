# NEO-Trajectory-Analysis.
This repository contains the MATLAB codebase developed for the Mission Analysis work package of the ESCUT (Early Survey, Characterisation and Updating of Trajectories) Group Design Project at Cranfield University (MSc Astronautics & Space Engineering, 2025-2026).

## Overview
The primary script (`NEO_Yogesh2.m`) performs an independent grid-search verification to determine the optimal Earth-to-NEO direct launch windows for Asteroid 292220 (2006 SU49). 

The analysis specifically maps the performance capabilities of the **Falcon Heavy (Expendable)** launch vehicle, evaluating characteristic energy ($C_3$), time of flight (TOF), and arrival braking $\Delta V$ to determine the maximum deliverable spacecraft mass.

## Features
* **Lambert Problem Solver:** Evaluates both short-way (prograde) and long-way (retrograde) heliocentric transfer arcs.
* **Launch Vehicle Performance Scaling:** Analytically derives payload capacity using the Tsiolkovsky rocket equation based on published Falcon Heavy Mars-injection reference data.
* **Porkchop Contour Generation:** Automatically generates high-resolution porkchop plots mapping the feasible mission space and optimal delivery windows.

## Files Included
* `NEO_Yogesh.m`: The main executable MATLAB script.
* `All_NEOS_ATA&TD_2018_2019.csv`: The Near-Earth Object database file required to pull the target ephemeris data.

## Dependencies
This script requires the following standard astrodynamics functions (developed during the Cranfield ATD module) to be present in your MATLAB path:
* `LambertArc_ATD2026.m`
* `EphSS_car.m`
* `ephNEO.m`
* `kep2car.m`
* `date2mjd2000.m` & `mjd20002date.m`
* `getAstroConstants.m`

## How to Run
1. Ensure all dependencies and the database file are in your active MATLAB directory and make sure the tool box is added to the path
2. Open `NEO_Yogesh2.m`.
3. Run the script. The console will output the grid search progress, summarize the optimum trajectories, and generate the corresponding porkchop plots.

---
*Author: Yogesh Andiyappan*
