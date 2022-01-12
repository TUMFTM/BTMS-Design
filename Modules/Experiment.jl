module Experiment
using JLD, Dates
using BatteryModel, CostModel
export experiment

## Functions
function experiment(DOE, Usecase, Climate)
    #=
    Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
    Created on: 11.01.2022
    Version: Julia 1.7.0

    Description: 
    This function executes the method for all configurations defined by the DOE, postprocesses the results and saves the data

    Input:    
    - DOE (Dict): containing the following parameters: 
        - Ebat (Float64): Battery size in Wh
        - Pthreshold (Float64): Peak shaving threshold in W
        - Cooling_thresholds (Range of Float64): investigated cooling thresholds in °C
        - Heating_thresholds (Range of Float64): investigated heating thresholds in °C
        - Passive (Vector of Strings): thermal designs of passive BTMS
        - Active (Vector of Strings): thermal designs of active BTMS
    - Usecase (Dict): containing the following parameters: 
        - name (String): powerprofile name
        - tbus (Vector of Float64): time of power profile in seconds
        - Pbus (Vector of Float64): power demand from charging buses in W 
        - dt (Float64): timestep of powerprofile in seconds
    - Climate (Dict): containing the following parameters: 
        - name (String): powerprofile name
        - Ta (Vector of Float64): ambient temperature in °C


    Output: 
    - res (Dict): containing all results from the battery simulation and the cost calcluation
    =#

    data = Dict()
    #Evaluate passive configurations
    data["Passive"] = Dict()
    for design in DOE["Passive"]
        data["Passive"][design] = Dict()
        data["Passive"][design]["Theat_req"] = NaN #Set default value
        for Theat in DOE["Heating_thresholds"]
            data["Passive"][design][Theat] = method(DOE["Ebat"],
                DOE["Pthreshold"], design, 1000, Theat, Usecase, Climate)  #The large cooling threshold represents passive cooling
            if data["Passive"][design][Theat]["pass_operability"] #Heating threshold does not need to be increased further
                data["Passive"][design]["Theat_req"] = Theat #Store required heating threshold
                break #and exit
            end
        end
    end

    #Evaluate actively cooled configurations
    data["Active"] = Dict()
    for design in DOE["Active"]
        data["Active"][design] = Dict()
        for Tcool in DOE["Cooling_thresholds"]
            data["Active"][design][Tcool] = Dict()
            data["Active"][design][Tcool]["Theat_req"] = NaN #Set default value
            for Theat in DOE["Heating_thresholds"]
                data["Active"][design][Tcool][Theat] = method(DOE["Ebat"], DOE["Pthreshold"], design, Tcool, Theat, Usecase, Climate)
                if data["Active"][design][Tcool][Theat]["pass_operability"] #Heating threshold does not need to be increased further
                    data["Active"][design][Tcool]["Theat_req"] = Theat #Store required heating threshold
                    break #and exit
                end
            end
        end
    end
    println("Experiment completed")

    #Add usecase information
    res = Dict()
    res["Data_raw"] = data
    res["DOE"] = DOE
    res["Climate"] = Climate["name"]
    res["Usecase"] = Usecase["name"]

    #Postprocess data
    postprocess(res)

    #save results
    timestamp = Dates.format(now(),"YYYYmmdd_HHMM")
    if Sys.iswindows()
        save("Results\\Data\\res_$timestamp.jld", res)
    else
        save("Results/Data/res_$timestamp.jld", res)
    end
    
    return res
end

function method(Ebat, Pthreshold, design, Tcool, Theat, Usecase, Climate)
    #=
    Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
    Created on: 11.01.2022
    Version: Julia 1.7.0

    Description: 
    The method combines the battery simulation with the cost model

    Input:    
    - Ebat (Float64): Battery size in Wh
    - Pthreshold (Float64): Peak shaving threshold in W
    - design (String): name of thermal design
    - Cooling_threshold (Float64): investigated cooling thresholds in °C
    - Heating_threshold (Float64): investigated heating thresholds in °C
    - Usecase (Dict): containing the following parameters: 
        - name (String): powerprofile name
        - tbus (Vector of Float64): time of power profile in seconds
        - Pbus (Vector of Float64): power demand from charging buses in W 
        - dt (Float64): timestep of powerprofile in seconds
    - Climate (Dict): containing the following parameters: 
        - name (String): powerprofile name
        - Ta (Vector of Float64): ambient temperature in °C

    Output: 
    - res (Dict): containing all results from the battery simulation and the cost calcluation
    =#

    #Print simulation status
    println("Design: $design, Tcool: $Tcool, Theat: $Theat")

    #Simulate battery operation
    res_sim = sim(Ebat, Pthreshold, design, Tcool, Theat,
                  Usecase["Pbus"], Usecase["dt"], Climate["Ta"])

    #Calculate system costs
    res_cost = cost(Ebat, design, Tcool, res_sim["teol"], res_sim["Ecool"], res_sim["Eheat"], res_sim["Eloss"])

    #Merge results
    res = merge(res_sim, res_cost)

    return res
end

function postprocess(res)
    #=
    Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
    Created on: 11.01.2022
    Version: Julia 1.7.0

    Description: 
    This function selects the results corresponding to the required heating threshold, 
    removes simulations that violate the safety constraint and restructures 
    the simulations into a dict where the values of the parameters are given for 
    the investigated cooling thresholds.  

    Input:    
    - res (Dict): contains the results of all simulated configurations

    Output: 
    - res (Dict): valid and restructured results 
    =#

    #Determine variables
    res["parameters"] = keys(res["Data_raw"]["Passive"]["base"][
                            res["DOE"]["Heating_thresholds"][1]])

    #Passively cooled designs
    res["Data"] = Dict()
    res["Data"]["Passive"] = Dict()
    for design in res["DOE"]["Passive"]
        res["Data"]["Passive"][design] = Dict()
        Theat_req = res["Data_raw"]["Passive"][design]["Theat_req"]

        #Limit solutions to those that satisfy the operability and safety constraints
        if isnan(Theat_req)
            map!(x->NaN, res["Data"]["Passive"][design])
        elseif !res["Data_raw"]["Passive"][design][Theat_req]["pass_safety"]
            map!(x->NaN, res["Data"]["Passive"][design])
        else
            res["Data"]["Passive"][design] = copy(res["Data_raw"]["Passive"][design][Theat_req])
            res["Data"]["Passive"][design]["Tcool"] = NaN
            res["Data"]["Passive"][design]["Theat_req"] = copy(Theat_req)
        end
    end

    #Actively cooled designs
    res["Data"]["Active"] = Dict()
    for design in res["DOE"]["Active"]
        res["Data"]["Active"][design] = Dict()

        #Limit solutions to those that satisfy the operability and safety constraints
        Theat_req = []
        for Tcool in res["DOE"]["Cooling_thresholds"]
            Theat = res["Data_raw"]["Active"][design][Tcool]["Theat_req"]
            if isnan(Theat)
                append!(Theat_req, NaN) #None of the heating thresholds were feasible
            elseif !res["Data_raw"]["Active"][design][Tcool][Theat]["pass_safety"]
                append!(Theat_req, NaN) #The maximum temperature was exceeded
            else
                append!(Theat_req, Theat)
            end
        end

        #restructure results
        res["Data"]["Active"][design]["Theat_req"] = copy(Theat_req)
        for param in res["parameters"]
            res["Data"]["Active"][design][param] = copy(
                            [isnan(Theat_req) ? NaN :
                            res["Data_raw"]["Active"][design][Tcool][Theat_req][param]
                            for (Tcool, Theat_req)
                            in zip(res["DOE"]["Cooling_thresholds"], Theat_req)
                            ])
        end

        #Deterime optimal cooling threshold
        Ctot = res["Data"]["Active"][design]["Ctot"]
        Ctot[isnan.(Ctot)] .= Inf
        _, Iopt = findmin(Ctot)
        Tcool_opt = res["DOE"]["Cooling_thresholds"][Iopt]
        Theat_opt =  res["Data"]["Active"][design]["Theat_req"][Iopt]
        res["Data"]["Active"][design]["opt"] = copy(res["Data_raw"]["Active"][design][Tcool_opt][Theat_opt])
        res["Data"]["Active"][design]["opt"]["Tcool"] = copy(Tcool_opt)
        res["Data"]["Active"][design]["opt"]["Theat_req"] = copy(Theat_opt)
    end

    return res
end
end
