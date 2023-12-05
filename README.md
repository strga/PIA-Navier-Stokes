This repository contains an implementation of the Semi-Implicit Method for Pressure-Linked Equations (SIMPLE) algorithm for solving the lid-driven cavity problem in computational fluid dynamics (CFD). It includes both the C++ code for simulation and a Python script for generating contour plots of the velocity fields.

**Overview:**

The lid-driven cavity problem is a classic test case for CFD algorithms. It involves fluid within a square cavity, set in motion by a moving lid, creating complex fluid flows. This repository offers a C++ implementation of the SIMPLE algorithm to solve the Navier-Stokes equations for this scenario. Additionally, there is a Python script provided for visualizing the results via contour plots. The implementation of the SIMPLE algorithm for the Finite Volume Method (FVM) is based on the following sources:

- [Math Online](https://mathonline.fme.vutbr.cz/download.aspx?id_file=1305)
- [University of Sydney Digital Repository](https://ses.library.usyd.edu.au/handle/2123/376)

**Content:**

Simulation.cpp: This C++ code implements the SIMPLE algorithm for the lid-driven cavity problem, calculating the velocity and pressure fields within the cavity.
plot_contours.py: A Python script used for generating contour plots from the C++ simulation output. It visualizes the components of velocity, u and v.
u_data.csv and v_data.csv: These are example output files from the C++ simulation that the Python script uses for plotting.
