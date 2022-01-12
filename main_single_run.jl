#=
Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
Created on: 11.01.2022
Version: Julia 1.7.0

Description: 
This script starts the simulation of a single battery thermal management systems 
for a given climate, cooling threshold and heating threshold. The result are used 
to generate Figure 4 and 5 in the research article entitled "Techno-economic Design 
of Battery Thermal Management Systems" and Fig. S1, S2, S3 and S4 in the supplement.
=#

if Sys.iswindows()
        push!(LOAD_PATH, string(pwd(), "\\Modules"))
else
        push!(LOAD_PATH, string(pwd(), "/Modules"))
end
using PreProcessing, BatteryModel, Visualization

## Load inputs
Usecase = loadPowerProfile("BLEL50_Kentridge") #Load bus charging demand power profile
Climate_Sgp = loadTemperature("Singapore_2015", Usecase["dt"]) #Load ambient temperature and solar irradiation
Climate_Muc = loadTemperature("Munich_2015", Usecase["dt"]) #Load ambient temperature and solar irradiation

## Define configuration
Ebat = 98e3 #Battery size in Wh
Pthreshold = 300e3 #Peak shaving threshold in W
Thermal_design = "base" #Thermal designs
Tcool = 35 #Cooling threshold in °C
Theat = 13 #Heating threshold in °C

## Simulate battery operation
@time res_Sgp = sim(Ebat, Pthreshold, Thermal_design, Tcool, Theat, Usecase["Pbus"], 
        Usecase["dt"], Climate_Sgp["Ta"]; fulloutput = true)
@time res_Muc = sim(Ebat, Pthreshold, Thermal_design, Tcool, Theat, Usecase["Pbus"], 
        Usecase["dt"], Climate_Muc["Ta"]; fulloutput = true)

## Visualize operation
Visualization.pgfplotsx() #Plotting backend for producing tikz figures
Plim_plot() #cell power limits, Fig. 4 in article
BatOperation_plot(Usecase, Climate_Sgp, res_Sgp) #Figure 5 in article
Ta_plot(Climate_Muc, Climate_Sgp) #Figure S1 in article supplement
Pbus_plot(Usecase) #Figure S2 in article supplement
BatOperation_plot(Usecase, Climate_Muc, res_Muc) #Figure S3 in article supplement
BatOperation_long(res_Muc, ["SOC", "Tc", "Pmin"], [(0,1), (10,55), (-55,-35)]) #Figure S4 in article supplement
