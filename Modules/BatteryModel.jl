module BatteryModel
using Statistics
using ControlAlgorithm, ElectricModel, ThermalModel, AgingModel
export sim

## Constants
const Tstart = 25 #Starting temperature in °C [Assumption]
const SOCstart = ElectricModel.SOC_max #Starting SOC [Assumption]
const taging_eval = 3600*24 #Time interval at which aging status is updated [Assumption]
const Qloss_EOL = 0.3 # EOL criterium [Assumption]
const Tmax = 60 #temperature safety limit in °C [Naumann et al.]
const tmax = 3600*24*365*20 # Maximum simulation duration [Assumption]

## Simulation
function sim(Ebat::Real, Pthreshold::Real, Design::String, Tcool::Real,
	Theat::Real, Pdem::Array{Float64,1}, dt::Real, Ta::Array{Float64,1}; fulloutput = false)
	#=
	Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
	Created on: 11.01.2022
	Version: Julia 1.7.0

	Description: 
	This function simulates the battery operation until the EOL criterium is met. 
	The structure of the model is visualized in Fig. 1 of the research article 
	entitled "Techno-economic Design of Battery Thermal Management Systems".

	Input:    
	- Ebat (Float64): battery size in Wh
	- Pthreshold (Float64): peak-shaving threshold in W
	- Design (String): name of battery thermal management system
	- Tcool (Float64): cooling threshold in °C
	- Theat (Float64): heating threshold in °C
	- Pdem (Vector of Float64): power demand in w
	- Ta (Vector of Float): ambient temperature in °C
	- dt (Float64): simulation timestep
	- fulloutput (Boolean): if true, the battery states at each timestep are 
			returned, if false, only metadata is returned to reduce the memory load

	Output:   
	- res (Dict): Dictionary containting the following parameters: 
			- pass_operability (boolean): flags if the operability constraint is met
			- pass_safety (boolean): flags if the safety constraint is met
			- teol (float64): battery life in years
			- Eloss (float64): energy consumption due to parasitic losses in kWh
			- Eheat (float64): energy consumption due to heating in kWh
			- Ecool (float64): energy consumption due to cooling in kWh
			- FEC (float64): full equivalent cycles per day
			- SOC_min (float64): minimum SOC reached in pu
			- SOC_avg (float64): average SOC during operation in pu
			- Tmin (float64): minimum cell temperature 
			- Tmax (float64): maximum cell temperature
			- Tavg (float64): average cell temperature
			- Crate_avg (float64): average C-restable
			- Ri_avg (float64): average internal resistance of a cell in Ω
			- Rout_avg (float64): average thermal resistance between housing and ambient in K/W
			- eta_avg (float64): average cell efficiency in pu
			- Qloss_cal_end (float64): capacity loss due to calendaric aging at simulation end in pu
			- Qloss_cyc_end (float64): capacity loss due to cyclic aging at simulation end in pu
			- Rinc_cal_end (float64): internal resistance increase due to calendaric aging at simulation end in pu
			- Rinc_cyc_end (float64): internal resistance increase due to cyclic aging at simulation end in pu
			- Rinc_tot_end (float64): total aging-induced internal resistance increase at simulation end in pu

			if the input flag "fulloutput" is True, the following output will be included as well: 
			- Tc (Vector of float64): cell temperature in °C
			- Th (Vector of float64): housing temperature in °C
			- SOC (Vector of float64): SOC in pu
			- Crate (Vector of float64): Crate in 1/h
			- Ri (Vector of float64): internal cell resistance in Ω
			- Rout (Vector of float64): thermal resistance between housing and ambient in K/W
			- Pbat (Vector of float64): battery power in W
			- Ploss (Vector of float64): ohmic losses in  a cell in W
			- Pheat (Vector of float64): heating power in W
			- Pcool (Vector of float64): cooling power in W
			- Pgrid (Vector of float64): power drawn from grid in W
			- Pmin (Vector of float64): maximum discharging power of a cell in W
			- Qloss_cal (Vector of float64): capacity loss due to calendaric aging in pu
			- Qloss_cyc (Vector of float64): capacity loss due to cyclic aging in pu
			- Rinc_cal (Vector of float64): internal resistance increase due to calendaric aging in pu
			- Rinc_cyc (Vector of float64): internal resistance increase due to cyclic aging in pu
			- eta (Vector of float64): efficiency in pu
	=#

	#Calculate top level parameters
	ncells = Ebat/ElectricModel.Unom/ElectricModel.Qnom #number of cells in the pack
	mhousing, Atop, ltop, Asides_conv, Asides_rad, Rins = TMparam(Design, ncells) #Top level thermal model parameters

	#Generate internal resistance interpolation functions
	fRi_ch, fRi_dch = ElectricModel.Ri_interpolants()

	#Preallocate variables that are logged every simulation timestep
	maxsteps = Int(tmax/dt) #Maximum number of steps in simulation
	SOC = Vector{Float64}(undef,maxsteps) #SOC
	Tcell = Vector{Float64}(undef,maxsteps) #Cell temperature
	Th = Vector{Float64}(undef,maxsteps) #Housing temperature
	Ri = Vector{Float64}(undef,maxsteps) #Cell internal resistance
	Rout = Vector{Float64}(undef,maxsteps) #Thermal resistance housing to ambient
	Crate = Vector{Float64}(undef,maxsteps) #Preallocate Crate
	Pgrid = Vector{Float64}(undef,maxsteps) #power drawn from grid
	Pcell = Vector{Float64}(undef,maxsteps) #power drawn from grid
	Pheat = Vector{Float64}(undef,maxsteps) #thermal system energy consumption
	Pcool = Vector{Float64}(undef,maxsteps) #thermal system energy consumption
	Ploss = Vector{Float64}(undef,maxsteps) #ohmic losses energy consumption
	Pmin = Vector{Float64}(undef,maxsteps) #maximum discharging power

	#Preallocate variables that are logged every aging timestep
	iaging = Int(taging_eval/dt) #Index range over which aging is evaluated
	maxagingsteps = Int(maxsteps/iaging) #Maximum number of aging steps in simulation
	Qloss_cal = Vector{Float64}(undef,maxagingsteps) #calendaric capacity loss
	Qloss_cyc = Vector{Float64}(undef,maxagingsteps) #cyclic capacity loss
	Rinc_cal = Vector{Float64}(undef,maxagingsteps) #calendaric internal resistance increase
	Rinc_cyc = Vector{Float64}(undef,maxagingsteps) #cyclic internal resistance increase

	#Set starting conditions of all states
	Tcell[1] = Tstart #Starting temperature cells
	Th[1] = Tstart #Starting temperature housing
	SOC[1] = SOCstart #Starting SOC
	Rinc_cal[1] = 0 #Starting calendaric resistance increase
	Rinc_cyc[1] = 0 #Starting cyclic resistance increase
	Qloss_cal[1] = 0 #Starting calendaric capacity loss
	Qloss_cyc[1] = 0 #Starting cyclic capacity loss
	j = 1 #Initial index of variables logged every aging timestep

	for i = 1:maxsteps-1
			#Get corresponding index of power and ambient temperature profiles
			iP = (i-1)%length(Pdem)+1 #power profile iterator
			iT = (i-1)%length(Ta)+1 #ambient temperature and irradiation iterator

			#Determine battery characteristics based on current state
			Uocv = ElectricModel.ocv(SOC[i])
			Ri_ch = fRi_ch(SOC[i],Tcell[i]) * (1+Rinc_cal[j]+Rinc_cyc[j])
			Ri_dch = fRi_dch(SOC[i],Tcell[i]) * (1+Rinc_cal[j]+Rinc_cyc[j])
			Q = ElectricModel.Qnom*(1-Qloss_cal[j]-Qloss_cyc[j])

			#Calculate power limits
			Pmin[i], Pmax = powerlim(Uocv, Ri_ch, Ri_dch, Q, SOC[i], Tcell[i], dt)

			#Apply peakshaving and thermal control algorithm
			Pgrid[i], Pcell[i], Pheat[i], Pcool[i] = control(
					Pthreshold, Pdem[iP], Pmin[i], Pmax, Tcell[i], Theat, Tcool, ncells)

			#Electric model
			Crate[i], SOC[i+1], Ploss[i], Ri[i] = electricmodel(Pcell[i], Uocv, Ri_ch, Ri_dch, Q, SOC[i], dt)

			#Thermal model
			Tcell[i+1], Th[i+1], Rout[i] = thermalmodel(Tcell[i], Th[i], Ploss[i], Pheat[i], Pcool[i],
																					Ta[iT], Rins, ncells, mhousing,
																					Atop, ltop, Asides_conv, Asides_rad, dt)

			#Aging model, updated at every aging timestep
			if i%iaging == 0
					#Extract aging influence factors for time window
					Taging = Tcell[(j-1)*iaging+1:j*iaging]
					SOCaging = SOC[(j-1)*iaging+1:j*iaging]
					Caging = Crate[(j-1)*iaging+1:j*iaging]

					#Update aging status
					Qloss_cal[j+1], Rinc_cal[j+1] = aging_cal(Qloss_cal[j], Rinc_cal[j], Taging, SOCaging, dt)
					Qloss_cyc[j+1], Rinc_cyc[j+1] = aging_cyc(Qloss_cyc[j], Rinc_cyc[j], Caging, SOCaging, dt)

					#Check if EOL was reached
					if (Qloss_cal[j] + Qloss_cyc[j]) > (Qloss_EOL)
							break
					else
							j += 1
					end
			end
	end

	#Crop outputs
	iend = Int(j*taging_eval/dt)
	Tcell = Tcell[1:iend]
	Th = Th[1:iend]
	SOC = SOC[1:iend]
	Crate = Crate[1:iend]
	Ri = Ri[1:iend]
	Rout = Rout[1:iend]
	Pcell = Pcell[1:iend]
	Ploss = Ploss[1:iend]
	Pheat = Pheat[1:iend]
	Pcool = Pcool[1:iend]
	Pgrid = Pgrid[1:iend]
	Pmin = Pmin[1:iend]
	Qloss_cal = Qloss_cal[1:j+1]
	Qloss_cyc = Qloss_cyc[1:j+1]
	Rinc_cal = Rinc_cal[1:j+1]
	Rinc_cyc = Rinc_cyc[1:j+1]

	#Compute results
	pass_operability = all(ncells*Pmin.<(Pthreshold-maximum(Pdem))) #operability constraint
	pass_safety = maximum(Tcell)<Tmax #safety constraint
	teol = j*taging_eval #battery lifetime in seconds
	teol_a = teol/3600/24/365 #battery lifetime in years
	Eloss = ncells*8.76*sum(Ploss)*dt/teol #Annual energy losses in kWh
	Eheat = 8.76*sum(Pheat)*dt/teol #Annual energy consumption of the heater in kWh
	Ecool = 8.76*sum(Pcool)*dt/teol #Annual energy consumption of the cooler in kWh
	eta = 1 .- Ploss./Pcell #Charging efficiency
	eta[eta.>1] = 1 ./eta[eta.>1] #Discharging efficiency
	FEC = sum(abs.(Crate))*dt/3600/(teol_a*365)/2 #Average full equivalent cycles per day

	#write top level outputs to result
	res = Dict()
	res["pass_operability"] = pass_operability
	res["pass_safety"] = pass_safety
	res["teol"] = teol_a
	res["Eloss"] = Eloss
	res["Eheat"] = Eheat
	res["Ecool"] = Ecool
	res["FEC"] = FEC
	res["SOC_min"] = minimum(SOC)
	res["SOC_avg"] = mean(SOC)
	res["Tmin"] = minimum(Tcell)
	res["Tmax"] = maximum(Tcell)
	res["Tavg"] = mean(Tcell)
	res["Crate_avg"] = mean(abs.(Crate))
	res["Ri_avg"] = mean(Ri)
	res["Rout_avg"] = mean(Rout)
	res["eta_avg"] = mean(filter(!isnan, eta))
	res["Qloss_cal_end"] = Qloss_cal[end]
	res["Qloss_cyc_end"] = Qloss_cyc[end]
	res["Rinc_cal_end"] = Rinc_cal[end]
	res["Rinc_cyc_end"] = Rinc_cyc[end]
	res["Rinc_tot_end"] = Rinc_cal[end]+Rinc_cyc[end]

	#full output (not included by default due to large amount of data)
	if fulloutput
			res["Tc"] = Tcell
			res["Th"] = Th
			res["SOC"] = SOC
			res["Crate"] = Crate
			res["Ri"] = Ri
			res["Rout"] = Rout
			res["Pbat"] = Pcell.*ncells
			res["Ploss"] = Ploss
			res["Pheat"] = Pheat
			res["Pcool"] = Pcool
			res["Pgrid"] = Pgrid
			res["Pmin"] = Pmin
			res["Qloss_cal"] = Qloss_cal
			res["Qloss_cyc"] = Qloss_cyc
			res["Rinc_cal"] = Rinc_cal
			res["Rinc_cyc"] = Rinc_cyc
			res["eta"] = eta
	end

	return res
end

end
