module ControlAlgorithm
using ThermalModel
export control

function control(Pthreshold, Pdem, Pmin, Pmax, T, Theat, Tcool, ncells)
    #=
    Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
    Created on: 11.01.2022
    Version: Julia 1.7.0

    Description: 
    This control algorithm ensures that the power drawn from the grid does not 
    exceed the peak shaving threshold by discharging the battery when the power 
    demand exceeds this threshold. 

    Input:    
    - Pthreshold (Float64): peak-shaving threshold in W
    - Pdem (Float64): power demand in w
    - Pmin (Float64): maximum discharging power in W
    - Pmax (Float64): maximum charging power in W
    - T (Float64): cell temperature in °C
    - Tcool (Float64): cooling threshold in °C
    - Theat (Float64): heating threshold in °C
    - ncells (Float64): number of cells in the pack

    Output:   
    - Pgrid (Float64): power drawn from grid in W
    - Pcell (Float64): power drawn from or supplied to cell in W
    - Pheat (Float64): applied heating power in W
    - Pcool (Float64): applied cooling power in W
    =#

    Pavail = max(0, Pthreshold - Pdem - ncells*Pmin) #Power available for cooling and heating
    if T < Theat
        Pheat = min(ThermalModel.Pheater, Pavail) #applied heating power
        Pcool = 0 #applied cooling power
    elseif T > Tcool
        Pheat = 0 #applied heating power
        Pcool = min(ThermalModel.Pcooler, Pavail) #applied cooling power
    else 
        Pheat = 0 #applied heating power
        Pcool = 0 #applied cooling power
    end

    Pbat = clamp(Pthreshold - Pdem - Pheat - Pcool, ncells*Pmin, ncells*Pmax) #power drawn from or supplied to the battery
    Pgrid = Pdem + Pheat + Pcool + Pbat #Power drawn from grid
    Pcell = Pbat/ncells #power drawn from or supplied to a single cell

    return Pgrid, Pcell, Pheat, Pcool
end
end
