# PIA-Navier-Stokes

This repository contains the implementation of the Semi-Implicit Method for Pressure-Linked Equations (SIMPLE) algorithm for solving the lid-driven cavity problem in computational fluid dynamics (CFD). The repository includes both the C++ code for the simulation and the Python script for generating contour plots of the velocity fields.

**Overview**:

The lid-driven cavity problem is a well-known test case for CFD algorithms. It involves fluid contained in a square cavity with a moving lid, which creates complex internal fluid flows. This repository provides a C++ implementation of the SIMPLE algorithm to solve the Navier-Stokes equations for this problem. Additionally, a Python script is included for visualizing the results through contour plots.


**Content**:

Simulation.cpp: C++ code implementing the SIMPLE algorithm for the lid-driven cavity problem. It calculates the velocity and pressure fields within the cavity.
plot_contours.py: Python script for generating contour plots from the C++ simulation output. It visualizes the u and v velocity components.
u_data.csv and v_data.csv: Example output files from the C++ simulation, used by the Python script for plotting.
