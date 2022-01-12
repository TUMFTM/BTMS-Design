module PreProcessing
using DelimitedFiles, Statistics, Interpolations
export loadPowerProfile, loadTemperature

function loadPowerProfile(PowerProfile)
    #=
    Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
    Created on: 11.01.2022
    Version: Julia 1.7.0

    Description: 
    This function reads in the powerprofile, determines the timestep, and 
    cuts and pastes the given profile to ensure it starts at midnight. 

    Input:    
    - Powerprofile (String): name of demand power profile

    Output:   
    - Usecase (Dict): containing the following parameters: 
        - name (String): powerprofile name
        - tbus (Vector of Float64): time of power profile in seconds
        - Pbus (Vector of Float64): power demand from charging buses in W 
        - dt (Float64): timestep of powerprofile in seconds
    =#

    ## Load bus charging power profile data
    if Sys.iswindows()
        MyPowerprofile = readdlm("Inputs\\Powerprofile_$PowerProfile.txt", ',');
    else
        MyPowerprofile = readdlm("Inputs/Powerprofile_$PowerProfile.txt", ',');
    end
    t_bus_raw = MyPowerprofile[:,1];
    P_bus_raw = MyPowerprofile[:,2];

    #Determine profile characteristics
    ndays = floor(t_bus_raw[end]/3600/24);
    dt_raw = unique(diff(t_bus_raw));
    if length(dt_raw)>1
        println("Error: Powerprofile does not have constant timestep")
        return
    else
        dt = Int8(dt_raw[1]); #Timestep
    end

    #Convert into profile starting at midnight
    t_before = collect(0:dt:t_bus_raw[1]-dt); #Create missing start of the day
    P_before = zeros(size(t_before)); #Create missing start of the day
    P_after = P_bus_raw[t_bus_raw.>=ndays*3600*24]; #Extract overhang after midnight of the last day

    t_bus = [t_before; t_bus_raw[t_bus_raw.<ndays*3600*24]]; #Add zeros before first bus departs and remove uncompleted day
    P_bus = [P_before; P_bus_raw[t_bus_raw.<ndays*3600*24]]; #Add zeros before first bus departs and remove uncompleted day
    P_bus[1:length(P_after)] += P_after; #Add overhang to the start of the profile

    #Write to dict
    Usecase = Dict()
    Usecase["name"] = PowerProfile
    Usecase["tbus"] = t_bus
    Usecase["Pbus"] = P_bus
    Usecase["dt"] = dt

    return Usecase
end

function loadTemperature(City, dt)
    #=
    Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
    Created on: 11.01.2022
    Version: Julia 1.7.0

    Description: 
    This function reads in the ambient temperature profile and interpolates it 
    to the timestep of the power demand profile.  

    Input:    
    - City (String): name of the ambient temperature profile
    - dt (Float64): timestep in seconds

    Output:   
    - Climate (Dict): containing the following parameters: 
        - name (String): powerprofile name
        - Ta (Vector of Float64): ambient temperature in °C
    =#

    #Load data
    if Sys.iswindows()
        TempData = readdlm("Inputs\\$City.csv",',')
    else
        TempData = readdlm("Inputs/$City.csv",',')
    end

    #Read out variables
    t = TempData[:,1]
    Ta = TempData[:,2]

    #Append start value of next year for correct interpolation
    push!(t,t[end]+3600)
    push!(Ta,Ta[1])

    #Map ambient temperature profile to timestep
    fT = LinearInterpolation(t,Ta)

    #Interpolate and write to Dict
    t_sim = 0:dt:t[end]-dt
    Climate = Dict()
    Climate["name"] = City
    Climate["Ta"] = fT(t_sim)

    return Climate
end

end
