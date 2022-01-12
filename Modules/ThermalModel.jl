module ThermalModel
using ElectricModel
export TMparam, thermalmodel

## constants
#Battery
const Ccell = 76.27 #Cell heat capacity J/K [Forgez et al.]
const Rth_in = 3.3 #Thermal resistance between cell and cell housing in K/W [Forgez et al.]
const mcell = 0.0845 # Cell mass in kg [Naumann et al.]
const m_celltopack = 1/0.552 #Increase in gravimetric energy density from cell to pack level [Löbberding et al.]
const cphousing = 896 #Specific heat capacity of 6061-T6 Aluminum at 20°C in J/kg/K [Lienhard et al. A heat transfer textbook p.714]
const Pcooler = 10e3 #Cooling power in W [Schimpe et al.]
const Pheater = 11.2e3 #Heating power in W [Schimpe et al.]
const COPcool = -3 #Coefficient of performance of cooling system in pu [Schimpe et al.]
const COPheat = 1 #Coefficient of performance of heating system in pu [Schimpe et al., Danish Energy Agency 2012]
const h = 2.2 #battery pack height in meters [Heliox]
const w = 0.8 #battery pack width in meters [Heliox]
const l = 0.8 #battery pack length in meters [Heliox]
const ϵ = 0.92 #Emissity of a painted surface [Lienhard et al. p.542]
const sfins = 0.05 #Fin spacing [Assumption]
const lfins = Dict("base"=>0.0, "fins"=>0.1, "ins"=>0.0) #fin length in meter [Assumption]
const t_ins = Dict("base"=>0.0, "fins"=>0.0, "ins"=>0.02) #insulation thickness in m [Assumption]
const kfoam = 0.035 #Thermal conductivity of expanded polystyrene in W/K/m [Lienhard et al. A heat transfer textbook p.718]

#Environment
const Pr = 0.707 #Prandtl number for air [Lienhard et al. A heat transfer textbook p.730]
const ν = 1.575e-5 #Kinematic viscosity of air at 300K in m/s^2 [Lienhard et al. A heat transfer textbook p.730]
const kair = 0.0264 #Thermal conductivity of air at 20°C in W/K/m [Lienhard et al. A heat transfer textbook p.730]
const g = 9.80665 #Standard acceleration of gravity [Lienhard et al. A heat transfer textbook p. 735]
const σ = 5.67e-8 #Stefan-Boltzmann constant [Lienhard et al. A heat transfer textbook p. 735]

## Functions
function TMparam(design, ncells)
    #=
    Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
    Created on: 11.01.2022
    Version: Julia 1.7.0

    Description: 
    This function calculates the top level thermal model parameters

    Input:    
    - design (String): name of thermal design
    - ncells (Float64): number of cells in the pack

    Output:   
    - mhousing (Float64): housing mass in kg
    - Atop (Float64): housing top surface area in m^2
    - ltop (Float64): Characteristic length of top surface in m 
    - Asides_conv (Float64): side surface area for natural convection in m^2 
    - Asides_rad (Float64): side surface area for radiation in m^2 
    - Rins (Float64): thermal resistance of insulation in k/W
    =#

    #Housing weight
    mhousing = ncells*mcell*(m_celltopack-1) #Weight of battery housing in kg

    #SES dimensions
    Atop = (l+2t_ins[design])*(w+2t_ins[design]) #Battery pack top surface in m^2
    Ptop = (2l+2w+8t_ins[design]) #Battery pack top surface perimeter in m
    ltop = Atop/Ptop #Battery pack top surface characteristic length for natural convection in m
    nfins = Ptop/sfins+4 #Number of cooling fins
    Asides_conv = h*(Ptop+2lfins[design]*nfins) #Side surface area for natural convection in m^2
    Asides_rad =  h*(Ptop+8lfins[design]) #Side surface for radiation in m^2

    #Insulation resistance
    Rins = t_ins[design]/kfoam/(Atop + Asides_rad) #Thermal resistance of insulation in K/W

    return mhousing, Atop, ltop, Asides_conv, Asides_rad, Rins
end

function thermalmodel(Tcell, Th, Ploss, Pheat, Pcool, Ta, Rins, ncells, mhousing, Atop, ltop, Asides_conv, Asides_rad, dt)
    #=
    Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
    Created on: 11.01.2022
    Version: Julia 1.7.0

    Description: 
    Thermal model with two thermal masses, corresponding to the cells and the housing.

    Input:    
    - Tcell (Float64): initial cell temperature in °C
    - Th (Float64): initial housing tempreature in °C
    - Ploss (Float64): ohmic losses in W
    - Pheat (Float64): applied heating power in W
    - Pcool (float64): applied cooling power in W
    - Ta (Float64): ambient temperature in °C
    - Rins (Float64): thermal resistance of insulation in k/W
    - ncells (Float64): number of cells
    - mhousing (Float64): housing mass in kg
    - Atop (Float64): housing top surface area in m^2
    - ltop (Float64): Characteristic length of top surface in m 
    - Asides_conv (Float64): side surface area for natural convection in m^2 
    - Asides_rad (Float64): side surface area for radiation in m^2 
    - dt (Float64): timestep duration in seconds

    Output:   
    - Tcell_new (Float64): updated cell temperature
    - Th_new (Float64): updated housing temperature
    - Rout (Float64): thermal resistance between housing and ambient in K/W
    =#

    #Cell temperature
    Tcell_new = Tcell + (Ploss + (Th-Tcell)/Rth_in)*dt/Ccell

    #Housing temperature
    Rout = ConvRad(Th, Ta, Atop, ltop, Asides_conv, Asides_rad, Rins) #Thermal resistance between surface and ambient
    Th_new = Th +  (Pheat*COPheat+Pcool*COPcool +
                    ncells*(Tcell-Th)/Rth_in +
                    (Ta-Th)/Rout
                    )*dt/mhousing/cphousing

    return Tcell_new, Th_new, Rout
end

function ConvRad(Th_Celsius, Ta_Celsius, Atop, ltop, Asides_conv, Asides_rad, Rins)
    #=
    Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
    Created on: 11.01.2022
    Version: Julia 1.7.0

    Description: 
    Function to calculate the thermal resistance due to convection and radiation between the housing and ambient air

    Input:    
    - Th_Celsius (Float64): Housing temperature in °C
    - Ta_Celsius (Float64): Ambient temperature in °C
    - Atop (Float64): housing top surface area in m^2
    - ltop (Float64): Characteristic length of top surface in m 
    - Asides_conv (Float64): side surface area for natural convection in m^2 
    - Asides_rad (Float64): side surface area for radiation in m^2 

    Output:   
    - Rout (Float64): thermal resistance between housing and ambient in K/W
    =#

    #Convert temperature to Kelvin
    Th = Th_Celsius + 273.15
    Ta = Ta_Celsius + 273.15

    #Convection top
    Ra_top = g/Ta*abs(Th-Ta)*ltop^3/ν^2 * Pr
    if Th>Ta #heated plate
        Nu_top = 0.15*(Ra_top*(1+(0.322/Pr)^(11/20))^(-20/11))^(1/3)
    else #cooled plate
        Nu_top = 0.6*(Ra_top*(1+(0.492/Pr)^(9/16))^(-16/9))^0.2
    end
    R_top = ltop/Nu_top/kair/Atop

    #Convection sides
    Ra_sides = g/Ta*abs(Th-Ta)*h^3/ν^2 * Pr
    Nu_sides = (0.825 + 0.387(Ra_sides^(1/6) / (1+(0.492/Pr)^(9/16))^(8/27)))^2
    R_sides = h/Nu_sides/kair/Asides_conv

    #Radiation
    R_rad = 1/(ϵ*σ*(Th^2+Ta^2)*(Th+Ta)*(Asides_rad+Atop))

    #Total
    if Th==Ta
        Rout = Rins + R_rad
    else
        Rout = Rins + (1/R_top+1/R_sides+1/R_rad)^-1
    end

    return Rout
end

end
