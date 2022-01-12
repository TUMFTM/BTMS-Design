#=
Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
Created on: 11.01.2022
Version: Julia 1.7.0

Description: 
This script starts the batch simulations of different battery thermal management 
systems (BTMS) in different climates. The results correspond to Figures 6,7 and 8 in the 
research article entitled "Techno-economic Design of Battery Thermal Management Systems"
and figure S5 in the supplement. 
=#

## Set path and required modules
if Sys.iswindows()
    push!(LOAD_PATH, string(pwd(), "\\Modules"))
else
    push!(LOAD_PATH, string(pwd(), "/Modules"))
end
using PreProcessing, Experiment, Visualization

## Load inputs
Usecase = loadPowerProfile("BLEL50_Kentridge") #Load power demand profile
Climate_Sgp = loadTemperature("Singapore_2015", Usecase["dt"]) #Load ambient temperature and solar irradiation data for Singapore
Climate_Muc = loadTemperature("Munich_2015", Usecase["dt"]) #Load ambient temperature and solar irradiation data for Munich

## Define Design of Experiment
DOE = Dict()
DOE["Ebat"] = 98e3 #Battery size in Wh
DOE["Pthreshold"] = 300e3 #Peak shaving threshold in W
DOE["Cooling_thresholds"] = 25:60 #Cooling thresholds
DOE["Heating_thresholds"] = 10:20 #Heating threshold in °C
DOE["Passive"] = ["base", "fins"] #Thermal designs with passive cooling
DOE["Active"] = ["base", "fins", "ins"] #Thermal designs with active cooling

## Run experiment
@time res_Sgp = experiment(DOE, Usecase, Climate_Sgp) #Run experiment for Singapore
@time res_Muc = experiment(DOE, Usecase, Climate_Muc) #Run experiment for Munich

## Plot results
Visualization.pgfplotsx() #Plotting backend for producing tikz figures
Theatplot(res_Muc) #Required heating threshold, Fig. 6 in article
Costplots(res_Muc, res_Sgp) #cost components over cooling threshold for different BTMS, Fig. 7 in article
cost_comparison(res_Muc, res_Sgp) #Cost comparison for Singapore, Fig. 8 in article
restable(res_Muc, res_Sgp) #Print detailed results, Tab. 3 in article
supp_plots(res_Muc, res_Sgp) #plot additional parameters, Fig. S5 in supplement