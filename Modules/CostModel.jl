module CostModel
using ThermalModel
export cost

## Constants
const c_bat = 420.0 # Battery cost in €/kWh [VDE]
const c_energy = 0.1909 # Electricity cost in €/Wh [BDEW]
const c_cooler_inv = 13000 # Investment of active cooling [Danish Energy Agency]
const c_cooler_om = 135 #operation and maintenance cost of cooler [Danish Energy Agency]
const t_cooler = 20 #20 years lifetime [Danish Energy Agency]
const c_heater_inv = 8000 # Investment of heater [Danish Energy Agency]
const c_heater_om = 100 # Operation and maintenance [Danish Energy Agency]
const t_heater = 30 # [Danish Energy Agency]
const c_fins_var = 25*0.86633 # EUR/m/m [Peng et al.]
const c_fins_fix = 187.5*0.86633 # Eur [Peng et al.]
const c_ins = 83*0.8633 #Eur/m/m/m [Alsayed et al.]
const t_fins_ins = 30 #lifetime of fins and insulation [assumption]
const r = 0.05 # Discount factor [Trocker et al.]

## Functions
function cost(Ebat, design, Tcool, t_eol, Ecool, Eheat, Eloss)
    #=
    Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
    Created on: 11.01.2022
    Version: Julia 1.7.0

    Description: 
    This function calculates all cost components that are influenced by the choice of BTMS

    Input:    
    - Ebat (Float64): Battery size in Wh
    - design (String): name of battery thermal management system
    - Tcool (Float64): cooling threshold in °C
    - t_eol (Float64): battery life in years
    - Ecool (Float64): Annual cooling energy consumption in kWh
    - Eheat (Float64): Annual heating energy consumption in kWh
    - Eloss (Float64): Annual energy consumption due to ohmic losses in kWh

    Output:   
    - res (Dict): dictionary containing the following parameters: 
        - Cinv (Float64): Investment costs of the BTMS in €/year
        - Cbat (Float64): Battery cost in €/year
        - Cene_cool (Float64): Cooling energy consumption cost in €/year
        - Cene_heat (Float64): Heating energy consumption cost in €/year
        - Cene_loss (Float64): Cost of energy consumption due to ohmic losses in €/year
        - Ctot (Float64): Sum of cost components in €/year
    =#

    res = Dict()

    #Investment cost
    res["Cinv"] = investment_costs(design, Tcool, Eheat)

    #Battery costs
    q_bat = (r*(1+r)^t_eol)/((1+r)^t_eol-1)
    res["Cbat"] = c_bat * Ebat/1000 * q_bat

    #Cooling system energy consumption
    res["Cene_cool"] = Ecool*c_energy

    #Heating system energy consumption
    res["Cene_heat"] = Eheat*c_energy

    #Ohmic losses energy consumption
    res["Cene_loss"] = Eloss*c_energy

    #Total cost
    res["Ctot"] = res["Cinv"] + res["Cbat"] + res["Cene_cool"] + res["Cene_heat"] + res["Cene_loss"]

    return res
end

function investment_costs(design, Tcool, Eheat)
    #=
    Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
    Created on: 11.01.2022
    Version: Julia 1.7.0

    Description: 
    This function calculates the investmentcost of the BTMS.

    Input:
    - design (String): name of battery thermal management system
    - Tcool (Float64): cooling threshold in °C
    - Eheat (Float64): Annual heating energy consumption in kWh

    Output:
    - Cbtms (Float64): Investment costs in €/year
    =#
    
    Cbtms = 0
    #cost for fins or insulation
    if design == "fins"
        nfins = 4+(2ThermalModel.l+2ThermalModel.w)/ThermalModel.sfins
        Afins = nfins*ThermalModel.h*2ThermalModel.lfins["fins"]
        Cbtms += (c_fins_fix+c_fins_var*Afins)*(r*(1+r)^t_fins_ins)/((1+r)^t_fins_ins-1)
    elseif design == "ins"
        V_ins = ThermalModel.h*ThermalModel.t_ins["ins"]*(2ThermalModel.l+2ThermalModel.w+4ThermalModel.t_ins["ins"])
        Cbtms += V_ins*c_ins*(r*(1+r)^t_fins_ins)/((1+r)^t_fins_ins-1)
    end

    #If a cooling system is installed
    if Tcool<=60
        Cbtms += c_cooler_om + c_cooler_inv*(r*(1+r)^t_cooler)/((1+r)^t_cooler-1)
    end
    
    #If installing a heater is required: 
    if Eheat>0
        Cbtms += c_heater_om + c_heater_inv*(r*(1+r)^t_heater)/((1+r)^t_heater-1)
    end

    return Cbtms
end
end
