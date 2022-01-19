# BTMS-Design

The battery thermal management system (BTMS) of a lithium-ion battery aims to prevent accelerated battery aging at elevated temperatures and reduced operability at low temperatures. Cooling or heating the battery prevents it from being operated outside the preferred temperature window but increases energy consumption, increases maintenance costs and requires an additional investment. Therefore, for a given use case, battery designers need to decide whether installing a heating system is required, if it is cost effective to install a cooling system, and how the battery should be thermally connected to the ambient air.

This repository provides the source code to a method for the techno-economic assessment of BTMS in different climates. The detailed documentation and a case study can be found in the following publication: https://doi.org/10.1016/j.est.2021.103832

## Getting Started
These instructions will get a copy of the project up and running on your local machine for development and testing purposes. 
 
### Prerequisites
Running the software requires a Julia distribution version 1.6.5 or above, which can be downloaded here: https://julialang.org/downloads/#long_term_support_release

Although Julia can be run from the commandline, we advise using Visual Studio (available here: https://code.visualstudio.com/download) with the Julia extension (via the VS Code marketplace: https://marketplace.visualstudio.com/items?itemName=julialang.language-julia). 
  
### Installing

All required code is obtained by cloning the repository. Additionally, the following packages need to be installed via the Pkg REPL, that can be entered by pressing <code>]</code>. Packages can then be installed with the following command <code>add <var>Package</var></code>.
  - DelimitedFiles
  - Statistics
  - Interpolations
  - JLD
  - Plots
  - PGFPlotsX
  - StatsPlots
  - NamedArrays

## Running the Model/Code
The repository containts two executable scripts: 
- main_single_run.jl starts the simulation of a single BTMS configuration in two climates and generates Figures 4,5,S1,S2,S3 and S4 from the research article.
- main.jl starts the batch simulation of all BTMS configurations defined in the design of experiment and generates figures 6,7,8 and S5 from the research article

## Contributing

If you would like to contribute to this work or have any feedback, please do not hesitate to contact me at olaf.teichert@tum.de

## Authors
Olaf Teichert
 
## License
This project is licensed under the LGPL License - see the LICENSE.md file for details
