# ERT_toolbox_CURE
ERT toolbox CURE project

This repository is part of an ongoing project which aims at characterizing the water content in landfills. The toolbox is written in python and does not require any special package to run. At the moment, the toolbox contains:

- Function to generate a dipole-dipole array that can be input in a Syscal Pro unit
- Function to generate a survey file for the Kragge landfill
- Example

The following is work in progress:

- Function to generate a Wenner-Schlumber array that can be input in a Syscal Pro unit
- Function to read the output binary file
- Function to plot raw data

Another toolbox needs to be created to process the data. This toolbox needs the python packages pygimli and pybert. The toolbox should contain:

- Function to reorganize the raw data into lines
- Function to invert the data
- Function to plot the data in 2D and maybe 3D
