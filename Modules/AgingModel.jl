module AgingModel
using Statistics
export aging_cal, aging_cyc

function aging_cal(Qloss::Real, Rinc::Real, T::Array{Float64,1}, SOC::Array{Float64,1}, dt::Real)
    #=
    Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
    Created on: 11.01.2022
    Version: Julia 1.7.0

    Description: 
    This function calculates the calendaric-aging-induced capacity loss and 
    internal resistance increase over a given interval based on the semi-empirical 
    aging model developed by Naumann et al.: https://doi.org/10.1016/j.est.2018.01.019

    Input:    
    - Qloss (Float64): aging-induced capacity loss at the start of the aging interval in pu
    - Rinc (Float64): aging-induced increase at the start of the aging interval in pu
    - T (Vector of Float64): temperature profile over aging interval in °C
    - SOC (Vector of Float64): SOC profile over aging interval in pu
    - dt (Float64): duration of aging interval in seconds

    Output:   
    - Qloss_new (Float64): aging-induced capacity loss after aging-interval in pu
    - Rinc_new (Float64): aging-induced internal resistance increase after aging-interval in pu
    =#

    Taging = max.(T,25).+273.15 #Modification: Calendar aging is constant below 25°C, where the aging model is not valid
    k_temp_Q = mean(map(T -> 1.2571e-05exp(-17126/8.3145*(1/T - 1/298.15)),Taging))
    k_temp_R = mean(map(T -> 3.4194e-10exp(-71827/8.3145*(1/T - 1/298.15)),Taging))

    k_SOC_Q = mean(map(SOC -> 2.85750(SOC-0.5)^3 + 0.60225, SOC))
    k_SOC_R = mean(map(SOC -> -3.3903(SOC-0.5)^2 + 1.56040, SOC))

    t_eq = (Qloss/k_SOC_Q/k_temp_Q)^2
    Qloss_new = k_SOC_Q*k_temp_Q*sqrt(t_eq+dt*length(SOC))
    Rinc_new = Rinc + k_SOC_R*k_temp_R*dt*length(SOC)

    return Qloss_new, Rinc_new
end

function aging_cyc(Qloss::Real, Rinc::Real, Crate::Array{Float64,1}, SOC::Array{Float64,1}, dt::Real)
    #=
    Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
    Created on: 11.01.2022
    Version: Julia 1.7.0

    Description: 
    This function calculates the cyclic-aging-induced capacity loss and 
    internal resistance increase over a given interval based on the semi-empirical 
    aging model developed by Naumann et al.: https://doi.org/10.1016/j.jpowsour.2019.227666

    Input:    
    - Qloss (Float64): aging-induced capacity loss at the start of the aging interval in pu
    - Rinc (Float64): aging-induced increase at the start of the aging interval in pu
    - Crate (Vector of Float64): Crate profile over aging interval in 1/h
    - SOC (Vector of Float64): SOC profile over aging interval in pu
    - dt (Float64): duration of aging interval in seconds

    Output:   
    - Qloss_new (Float64): aging-induced capacity loss after aging-interval in pu
    - Rinc_new (Float64): aging-induced internal resistance increase after aging-interval in pu
    =#
    
    DODs = rfcounting(SOC) #Rainflow counting to determine effective DOD cycles
    FEC = sum(DODs)
    k_DOD_Q = sum(map(DOD -> 4.0253(DOD-0.6)^3 + 1.09230, DODs).*DODs)/FEC
    k_DOD_R = sum(map(DOD -> 6.8477(DOD-0.5)^3 + 0.91882, DODs).*DODs)/FEC

    Caging = abs.(Crate)
    k_C_Q = sum(map(Crate-> 0.0971 + 0.063Crate       , Caging).*Caging*dt/3600/2)/FEC
    k_C_R = sum(map(Crate-> 0.0023 - 0.0018min(Crate,1), Caging).*Caging*dt/3600/2)/FEC#According to Naumann fig. 6c

    FEC_eq_Q = (100Qloss/k_C_Q/k_DOD_Q)^2
    Qloss_new = 0.01k_C_Q*k_DOD_Q*sqrt(FEC_eq_Q+FEC)
    Rinc_new = Rinc + 0.01k_C_R*k_DOD_R*FEC

    return Qloss_new, Rinc_new
end

function rfcounting(SOC::Array{Float64,1})
    #=
    Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
    Created on: 11.01.2022
    Version: Julia 1.7.0

    Description: 
    This function calculates the effective depth of discharge (DOD) cycles over
    a given aging interval using the rainflow counting algorithm Algorithm I 
    described by Downing et al. https://doi.org/10.1016/0142-1123(82)90018-4. 
    The results are used to calculate capacity loss and internal resistance 
    increase resulting from cyclic aging. 
    
    Input:    
    - SOC (Vector of Float64): SOC profile over aging interval in pu

    Output:   
    - DODs (Vector of Float64): Effective DOD cycles over given aging interval in pu
    =#    

    #Remove constant SOC phases
    SOC_nozero = SOC[vcat(diff(SOC).!=0,true)]

    #Slice to start and end with the maximum SOC
    _, I = findmax(SOC_nozero)
    SOC_sorted = vcat(SOC_nozero[I:end],SOC_nozero[1:I])

    #Find extremas
    slope = diff(SOC_sorted)
    is_extremum = vcat(true, (slope[1:end-1].*slope[2:end]).<=0., true)
    SOCs = SOC_sorted[is_extremum]

    #Find DODs
    DODs = []
    index = 2
    while length(SOCs) > index
        prevDOD = abs(SOCs[index]-SOCs[index-1])
        nextDOD = abs(SOCs[index]-SOCs[index+1])
        if nextDOD < prevDOD
            index += 1
        else
            push!(DODs, prevDOD)
            deleteat!(SOCs, (index-1):index)
            index = 2
        end
    end
    return DODs
end

end
