module Visualization
using Dates, Plots, StatsPlots, NamedArrays
using ElectricModel
export BatOperation_plot, Plim_plot, Costplots, Theatplot, cost_comparison, restable, Ta_plot, Pbus_plot, BatOperation_long, supp_plots

if Sys.iswindows()
        const resultpath = "Results\\Figures\\"
else
        const resultpath = "Results/Figures/"
end
const linestyle = Dict("base"=>:solid, "fins"=>:dot, "ins"=>:dash)

function Plim_plot()
	#=
	Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
	Created on: 11.01.2022
	Version: Julia 1.7.0

	Description: 
	Generates plot of cell power limits over SOC for different cell temperatures

	Output:   
	- Figure 4 in the research article entitled "Techno-economic Design of Battery Thermal Management Systems"
	=#

	Q = ElectricModel.Qnom
	fRi_ch, fRi_dch = ElectricModel.Ri_interpolants()
	SOC = 0:0.01:1
	T = [10 20 40]
	dt = 10
	Tlabels = ["$Temp  in deg C" for Temp in T]
	Pmin = zeros(length(SOC), length(T))
	Pmax = zeros(length(SOC), length(T))
	for i in 1:length(SOC)
		for j in 1:length(T)
			Uocv = ElectricModel.ocv(SOC[i])
			Ri_ch = fRi_ch(SOC[i],T[j])
			Ri_dch = fRi_dch(SOC[i],T[j])
			Pmin[i,j], Pmax[i,j] = ElectricModel.powerlim(Uocv, Ri_ch, Ri_dch, Q, SOC[i], T[j], dt)
		end
	end

	mycolors = range(colorant"blue", stop=colorant"red",length=length(T))
	p = plot(SOC, Pmin, lc=mycolors', labels =  Tlabels)
	plot!(p, SOC, Pmax, lc=mycolors', labels = "")
	xlabel!(p, "SOC")
	ylabel!(p, "P_cell in Watt")
	plot!(p, legend=:outertopright)

	savefig(p, string(resultpath, "Plims.tex"))
	return p
end

function BatOperation_plot(Usecase, Climate, res)
	#=
	Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
	Created on: 11.01.2022
	Version: Julia 1.7.0

	Description: 
	Generates plot of battery operation

	Input: 
	- Usecase (Dict): containing the following parameters: 
		- name (String): powerprofile name
		- tbus (Vector of Float64): time of power profile in seconds
		- Pbus (Vector of Float64): power demand from charging buses in W 
		- dt (Float64): timestep of powerprofile in seconds
	- Climate (Dict): containing the following parameters: 
			- name (String): powerprofile name
			- Ta (Vector of Float64): ambient temperature in °C
	- res (Dict): containing results from the battery simulation

	Output:   
	- Figure 5 or Figure S3 in the research article entitled "Techno-economic Design of Battery Thermal Management Systems"
	=#

	default(lw=0.5) #set line width

	#Define start and end time of plot
	tstart = DateTime(1,1,1,6)
	tend = DateTime(1,1,2,6)
	plotinterval = 6 #only plot every 6th value to reduce fig. size

	#Crop to plotting range
	t = Dates.value(tstart):Usecase["dt"]*1000*plotinterval:Dates.value(tend)
	I1 = Int(Dates.value(tstart)/1000/Usecase["dt"])
	I2 = Int(Dates.value(tend)/1000/Usecase["dt"])
	Pbus = Usecase["Pbus"][I1:plotinterval:I2]
	Pgrid = res["Pgrid"][I1:plotinterval:I2]
	Pbat = res["Pbat"][I1:plotinterval:I2]
	SOC = res["SOC"][I1:plotinterval:I2]
	Tc = res["Tc"][I1:plotinterval:I2]
	Th = res["Th"][I1:plotinterval:I2]
	Ta = Climate["Ta"][I1:plotinterval:I2]

	#Plot power
	p1 = plot(t,Pbus./1000,label="Pbus")
	plot!(p1,t,Pgrid./1000,label="Pgrid")
	plot!(p1,t,Pbat./1000,label="Pbat")
	ylims!(p1, (-439,750))
	ylabel!("kW")

	#Plot SOC
	p2 = plot(t, SOC, label="SOC", legend=:right)
	ylims!(p2, (0,1))

	#Plot Temperature
	p3 = plot(t,Tc,label="T_c")
	plot!(p3,t,Th,label="T_h")
	plot!(p3,t[1:Int(360/plotinterval):end],Ta[1:Int(360/plotinterval):end],label="T_a")
	ylabel!(p3,"deg C")
	#ylims!(p3, (24,40))

	#Add common xaxislabels
	tstep = Int(2*3600/Usecase["dt"]/plotinterval)
	xtickspos = t[1:tstep:end]
	xticklabels = Dates.format.(tstart:Second(tstep*Usecase["dt"]*plotinterval):tend,"Ip")
	xlims!(p1, (xtickspos[1], xtickspos[end]))
	xlims!(p2, (xtickspos[1], xtickspos[end]))
	xlims!(p3, (xtickspos[1], xtickspos[end]))
	xticks!(p1, xtickspos, repeat([""], length(xtickspos)))
	xticks!(p2, xtickspos, repeat([""], length(xtickspos)))
	xticks!(p3, xtickspos, xticklabels)

	#Combine plots
	ptot = plot(p1, p2, p3, layout = (3,1), link = :x)

	#Save
	savefig(ptot, string(resultpath, "$(Climate["name"])_Batops.tex"))
	display(ptot)
	return nothing
end

function Theatplot(res)
	#=
	Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
	Created on: 11.01.2022
	Version: Julia 1.7.0

	Description: 
	Generates plot of the required heating threshold of different thermal designs

	Input: 
	- res (Dict): results from the batch simulation

	Output:   
	- Figure 5 in the research article entitled "Techno-economic Design of Battery Thermal Management Systems"
	=#

	p = paramplot(res, "Theat_req")
	plot!(p,
				xlabel = "T_cool in deg C",
				ylabel = "Required heating threshold in deg C",
				ylims  = (minimum(res["DOE"]["Heating_thresholds"])-1,
									maximum(res["DOE"]["Heating_thresholds"])+1)
				)

	#Save fig
	savefig(p, string(resultpath, "Theat.tex"))
	display(p)
	return nothing
end

function Costplots(res_Muc, res_Sgp)
	#=
	Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
	Created on: 11.01.2022
	Version: Julia 1.7.0

	Description: 
	Generates plots of cost components over the cooling threshold for all thermal designs

	Input: 
	- res_Muc (Dict): results from the batch simulation in the Munich climate
	- res_Sgp (Dict): results from the batch simulation in the Singapore climate

	Output:   
	- Figure 7 in the research article entitled "Techno-economic Design of Battery Thermal Management Systems"
	=#

	params = ["Cbat", "Cene_cool", "Cene_heat", "Cene_loss", "Ctot"]
	ymaxs =  [9000,750,1500,1500,12000]

	ptot = paramplots(res_Muc, res_Sgp, params, ymaxs)

	savefig(ptot, string(resultpath, "costplot.tex"))
	display(ptot)
	return nothing
end

function paramplots(res_Muc, res_Sgp, params, ymaxs)
	#=
	Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
	Created on: 11.01.2022
	Version: Julia 1.7.0

	Description: 
	Generates plot of multiple parameters for 2 climates

	Input: 
	- res_Muc (Dict): results from the batch simulation in the Munich climate
	- res_Sgp (Dict): results from the batch simulation in the Singapore climate
	- params (Vector of Strings): parameters to be plotted
	- ymaxs (Vector of Float64): upper ylims for plots

	Output:   
	- Figure
	=#

	plots_Muc = [paramplot(res_Muc, param; ymax = y) for (param, y) in zip(params,ymaxs)]
	plots_Sgp = [paramplot(res_Sgp, param; ymax = y) for (param, y) in zip(params,ymaxs)]

	title!(plots_Muc[1],"Munich")
	title!(plots_Sgp[1],"Singapore")
	xlabel!(plots_Muc[end],"T_cool in deg C")
	xlabel!(plots_Sgp[end],"T_cool in deg C")
	xticks!.(plots_Muc[1:end-1], :none)
	xticks!.(plots_Sgp[1:end-1], :none)
	ylabel!.(plots_Muc, params)
	yticks!.(plots_Sgp,:none)

	ptot = plot(collect(Iterators.flatten(zip(plots_Muc, plots_Sgp)))...,
							layout=(length(params),2),margin=-5*Plots.mm)

	return ptot
end

function paramplot(res, param; ymin = 0, ymax=Inf)
	#=
	Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
	Created on: 11.01.2022
	Version: Julia 1.7.0

	Description: 
	Generates plot of of a single parameter over the cooling threshold for all thermal designs

	Input: 
	- res (Dict): results from the batch simulation
	- params (String): parameter to be plotted
	- ymin (Float64): lower ylimit
	- ymax (Float64): upper ylimit

	Output:   
	- Plot
	=#

	p = plot(ylims=(ymin,ymax), legend=:none)
	for design in res["DOE"]["Active"]
					plot!(p, res["DOE"]["Cooling_thresholds"],
									res["Data"]["Active"][design][param],
									line = linestyle[design], color=:black,
									label = design)
	end

	#Add annotations of BTMS
	xs =   [res["Data"]["Passive"]["base"]["Tmax"]
					res["Data"]["Passive"]["fins"]["Tmax"]
					res["Data"]["Active"]["base"]["opt"]["Tcool"]
					res["Data"]["Active"]["fins"]["opt"]["Tcool"]
					res["Data"]["Active"]["ins"]["opt"]["Tcool"]
					]
	ys =   [res["Data"]["Passive"]["base"][param]
					res["Data"]["Passive"]["fins"][param]
					res["Data"]["Active"]["base"]["opt"][param]
					res["Data"]["Active"]["fins"]["opt"][param]
					res["Data"]["Active"]["ins"]["opt"][param]
					]
	annotations = [(n, :bottom) for n in ["I", "II", "III", "IV", "V"]]

	scatter!(p,xs, ys, series_annotations = annotations,
					marker = :black, label=:none)
	return p
end

function cost_comparison(res1, res2)
	#=
	Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
	Created on: 11.01.2022
	Version: Julia 1.7.0

	Description: 
	Generates stacked barchart comparing the cost of BTMS at their optimal cooling
	threshold and required heating threshold

	Input: 
	- res1 (Dict): results from the batch simulation for the first climate
	- res1 (Dict): results from the batch simulation for the second climate

	Output:   
	- Figure 8 in the research article entitled "Techno-economic Design of Battery Thermal Management Systems"
	=#

	#Define cost components and extract data from res struct
	costcomponents = ["Cinv" "Cene_loss" "Cene_heat" "Cene_cool" "Cbat" "Ctot"]
	c_pas1 = [[res1["Data"]["Passive"][design][c]
					for c in costcomponents] for design in res1["DOE"]["Passive"]]
	c_act1 = [[res1["Data"]["Active"][design]["opt"][c]
					for c in costcomponents] for design in res1["DOE"]["Active"]]
	c_pas2 = [[res2["Data"]["Passive"][design][c]
					for c in costcomponents] for design in res2["DOE"]["Passive"]]
	c_act2 = [[res2["Data"]["Active"][design]["opt"][c]
					for c in costcomponents] for design in res2["DOE"]["Active"]]
	costs = vcat(c_pas1, c_act1, [zeros(1, length(costcomponents))], c_pas2, c_act2)
	costs = vcat(costs...) #Convert to required plotting format

	p = groupedbar(["I","II","III","IV","V","","a","b","c","d","e"],
							costs[:,1:end-1],
							bar_position = :stack,
							labels = hcat(costcomponents[1:end-1]...),
							ylabel= "EUR/year",
							extra_kwargs = :subplot,
							ylims = (0,12000),
							legend=:top,
							legend_columns=-1
							)

	signsymbol = Dict(1=>"+",-1=>"")
	for i in 2:(length(c_pas1)+length(c_act1))
					change = 100*(costs[i,end]-costs[1,end])/costs[1,end]
					annotate!(p, i-0.5, 1.1costs[i,end],
					text("$(signsymbol[sign(change)])$(round(change, digits=1))%",6),
					label=:none)
	end
	for i in 8:6+(length(c_pas1)+length(c_act1))
					change = 100*(costs[i,end]-costs[7,end])/costs[7,end]
					annotate!(p, i-0.5, 1.1costs[i,end],
					text("$(signsymbol[sign(change)])$(round(change, digits=1))%",6),
					label=:none)
	end

	savefig(p, string(resultpath, "cost_comparison.tex"))
	display(p)
	return nothing
end

function restable(res1, res2)
	#=
	Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
	Created on: 11.01.2022
	Version: Julia 1.7.0

	Description: 
	Generates table containing the key properties of the BTMS at their optimal cooling
	threshold and required heating threshold

	Input: 
	- res1 (Dict): results from the batch simulation for the first climate
	- res1 (Dict): results from the batch simulation for the second climate

	Output:   
	- Table 3 in the research article entitled "Techno-economic Design of Battery Thermal Management Systems"
	=#

	param = ["Tcool","Theat_req", "Ctot", "Cbat", "Cinv", "Cene_cool", "Cene_heat", "Cene_loss", "teol", "Rinc_tot_end", "Tmin", "Tavg","Tmax"]
	designs = [
					res1["Data"]["Passive"]["base"],
					res1["Data"]["Passive"]["fins"],
					res1["Data"]["Active"]["base"]["opt"],
					res1["Data"]["Active"]["fins"]["opt"],
					res1["Data"]["Active"]["ins"]["opt"],
					res2["Data"]["Passive"]["base"],
					res2["Data"]["Passive"]["fins"],
					res2["Data"]["Active"]["base"]["opt"],
					res2["Data"]["Active"]["fins"]["opt"],
					res2["Data"]["Active"]["ins"]["opt"]
					]
	paramtab = [[r[p] for p in param] for r in designs]

	myTable = NamedArray(hcat(paramtab...))
	setnames!(myTable, param, 1)
	setnames!(myTable, ["I","II","III","IV","V","a","b","c","d","e"], 2)
	@show myTable

	return nothing
end

function Ta_plot(Climate_Muc, Climate_Sgp)
	#=
	Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
	Created on: 11.01.2022
	Version: Julia 1.7.0

	Description: 
	Generates plot of ambient temperature in both considered climates

	Input: 
	- Climate_Muc (Dict): Temperature profile in Munich containing: 
		- name (String): powerprofile name
		- Ta (Vector of Float64): ambient temperature in °C
	- Climate_Sgp (Dict): Temperature profile in Singapore containing: 
		- name (String): powerprofile name
		- Ta (Vector of Float64): ambient temperature in °C

	Output:   
	- Fig. S1 in the supplemental material to the research article entitled 
	"Techno-economic Design of Battery Thermal Management Systems"
	=#

	p = plot(grid=:y,
					ylabel = "Ta in degC"
					)
	vline!(p, [1:11], color = "gray", lw=0.25, label=:none)

	t = range(0,12,length=365*2)
	for Climate in [Climate_Muc Climate_Sgp]
					mins = [minimum(Climate["Ta"][(360*24*(n-1)+1):360*24*n]) for n in 1:365]
					maxs = [maximum(Climate["Ta"][(360*24*(n-1)+1):360*24*n]) for n in 1:365]
					minmaxs = collect(Iterators.flatten(zip(maxs, mins)))
					plot!(p,t,
									minmaxs,
									grid=:y,
									xticks=:none,
									tick_direction = :none,
									label = Climate["name"])
	end
	plot!(p, xlims=(0,12),
						xticks=(0.5:11.5,
									["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
									),
						tick_direction = :none,
						extra_kwargs=:subplot,
						legend=:outertop,
						legend_columns=-1
						)

	savefig(p, string(resultpath, "Ta.tex"))
	display(p)
	return nothing
end

function Pbus_plot(Usecase)
	#=
	Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
	Created on: 11.01.2022
	Version: Julia 1.7.0

	Description: 
	Generates plot of power demand resulting from charging buses

	Input: 
	- Usecase (Dict): containing the following parameters: 
        - name (String): powerprofile name
        - tbus (Vector of Float64): time of power profile in seconds
        - Pbus (Vector of Float64): power demand from charging buses in W 
        - dt (Float64): timestep of powerprofile in seconds

	Output:   
	- Fig. S2 in the supplemental material to the research article entitled 
	"Techno-economic Design of Battery Thermal Management Systems"
	=#

	default(lw=0.5)

	#Define start and end time of plot
	hstart = 6
	plotinterval = 6

	#Crop to plotting range
	Pbus = vcat(Usecase["Pbus"][360*hstart:end], Usecase["Pbus"][1:360*hstart-1])         #Slice bus power profile
	Pbus = Pbus[1:plotinterval:end]

	#Plot power
	p = plot(Pbus./1000, label=:none,
					ylabel = "Bus charging power in kW",
					xticks = (0:12*360/plotinterval:length(Pbus),
									Dates.format.(DateTime(1,1,1,6):Second(12*3600) :DateTime(1,1,6,6),"Ip"))
									)

	#Save
	savefig(p, string(resultpath, "Pbus.tex"))
	display(p)
	return nothing
end

function BatOperation_long(res, params, ylimss)
	#=
	Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
	Created on: 11.01.2022
	Version: Julia 1.7.0

	Description: 
	Generates plot of parameters over battery lifetime 

	Input: 
	- res (Dict): results from a single simulation
	- params (Vector of Strings): to be plotted parameters
	- ylimss (Vector of Tuples): ylims of each parameter

	Output:   
	- Fig. S4 in the supplemental material to the research article entitled 
	"Techno-economic Design of Battery Thermal Management Systems"
	=#

	nstep = 5
	t = 0:nstep/365/2:(res["teol"]-nstep/365)
	ps = []

	for (param, param_ylim) in zip(params,ylimss)
					mins = [minimum(res[param][(360*24*(n-nstep)+1):360*24*n]) for n in nstep:nstep:Int(res["teol"]*365)]
					maxs = [maximum(res[param][(360*24*(n-nstep)+1):360*24*n]) for n in nstep:nstep:Int(res["teol"]*365)]
					minmaxs = collect(Iterators.flatten(zip(maxs, mins)))

					p = plot()
					vline!(p, [1:Int(ceil(res["teol"]))], color = "gray", lw=0.05,
									label=:none, ylims=param_ylim)
					plot!(p,t,minmaxs,label=:none,grid=:y,
									xticks=:none,
									tick_direction = :none,
									xlabel="years",
									ylabel=param)
					if param == params[end]
									xticks!(p, 0.5:Int(ceil(res["teol"])), string.(1:Int(ceil(res["teol"]))))
					end
					push!(ps,p)
	end

	ptot = plot(ps...,layout=(length(params),1))

	savefig(ptot, string(resultpath, "Batopslong.tex"))
	display(ptot)
	return nothing

end

function supp_plots(res_Muc, res_Sgp)
	#=
	Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
	Created on: 11.01.2022
	Version: Julia 1.7.0

	Description: 
	Generates plots of parameters over the cooling threshold for all thermal designs

	Input: 
	- res_Muc (Dict): results from the batch simulation in the Munich climate
	- res_Sgp (Dict): results from the batch simulation in the Singapore climate

	Output:   
	- Fig. S5 in the supplemental material to the research article entitled 
	"Techno-economic Design of Battery Thermal Management Systems"
	=#

	params = ["teol", "Rinc_tot_end", "Tmin", "Tavg", "Tmax"]
	ymaxs =  [Inf, Inf, Inf, Inf, Inf]

	ptot = paramplots(res_Muc, res_Sgp, params, ymaxs)

	savefig(ptot, string(resultpath, "supp_plots.tex"))
	display(ptot)
	return nothing
end
end
