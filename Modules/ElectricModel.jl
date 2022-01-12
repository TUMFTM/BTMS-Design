module ElectricModel
using Interpolations
export electricmodel, powerlim

## Constants
const Qnom = 3.0 # Cell nominal capacity in Ah [Naumann et al.]
const Unom = 3.2 # Cell nominal voltage in V [Naumann et al.]
const Umin = 2 # Cell minimum voltage in V [Naumann et al.]
const Umax = 3.6 # Cell maximum voltage in V [Naumann et al.]
const Imax_cont = 3 # Maximum charging current [Naumann et al.]
const Imin_cont = -20 # Maximum discharge current [Naumann et al.]
const SOC_max = 0.95 # Upper SOC limit [Rodrigo et al.]
const SOC_min = 0.05 # Lower SOC limit [Rodrigo et al.]

## Functions
function electricmodel(Pcell, Uocv, Ri_ch, Ri_dch, Q, SOC, dt)
  #=
  Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
  Created on: 11.01.2022
  Version: Julia 1.7.0

  Description: 
  This function calculates the Crate, the ohmic losses and the SOC at the end of 
	the timestep for an equivalent circuit model with a single DC resistance

  Input: 
	- Pcell (Float64): power drawn from or supplied to cell in W
	- Uocv (Float64): open circuit voltage in V
	- Ri_ch (Float64): charging resistance in Ω
	- Ri_dch (Float64): discharging resistance in Ω
	- Q (Float64): cell capacity in Ah
	- SOC (Float64): SOC at beginning of timestep in pu
	- dt (Float64): timestep duration in seconds

	Output: 
	- Crate (Float64): C-rate
	- SOC_new (Float64): SOC after timestep
	- Ploss (Float64): Ohmic losses in W
	- Ri (Float64): Cell internal resistance in Ω
  =#

  Ri = sign(Pcell) > 0 ? Ri_ch : Ri_dch
  Ibat = (-Uocv+sqrt(Uocv^2+4Ri*Pcell))/2Ri #Cell current in A
  Crate = Ibat/Q #Cell C-rate in 1/h
  SOC_new = SOC + Crate*dt/3600 #New cell SOC in pu
  Ploss = Ibat^2*Ri #Generated ohmic losses in W

  return Crate, SOC_new, Ploss, Ri
end

function powerlim(Uocv, Ri_ch, Ri_dch, Q, SOC, T, dt)
	#=
  Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
  Created on: 11.01.2022
  Version: Julia 1.7.0

  Description: 
  This function calculates the charging and discharging power limits

  Input: 
	- Uocv (Float64): open circuit voltage in V
	- Ri_ch (Float64): charging resistance in Ω
	- Ri_dch (Float64): discharging resistance in Ω
	- Q (Float64): cell capacity in Ah
	- SOC (Float64): SOC at beginning of timestep in pu
	- T (Float64): temperature in °C
	- dt (Float64): timestep duration in seconds

	Output: 
	- Pmin (Float64): maximum discharging power in W
	- Pmax (Float64): maximum charging power in W 
  =#

	#Discharging power limit
	Imin_SOC = min(0,(SOC_min-SOC)*3600*Q/dt) #Current limit based on minimum SOC
	Imin_V = (Umin-Uocv)/Ri_dch #Current limit based on minimum voltage
	Imin = max(Imin_cont, Imin_V, Imin_SOC) #Discharge current limit
	Pmin = (Uocv+Ri_dch*Imin)*Imin #Discharge power limit

	#Charging power limit
	Imax_SOC = max(0,(SOC_max-SOC)*3600*Q/dt) #Current limit based on maximum SOC
	Imax_V = (Umax-Uocv)/Ri_ch #Current limit based on maximum voltage
	Imax = min(Tderating(T)*Imax_cont, Imax_V, Imax_SOC) #Charge current limit
	Pmax = (Uocv+Ri_dch*Imax)*Imax #Charge power limit

	return Pmin, Pmax
end

function Tderating(T)
	#=
  Designed by: Olaf Teichert (Institute of Automotive Technology, Technical University of Munich)
  Created on: 11.01.2022
  Version: Julia 1.7.0

  Description: 
  This function calculates the maximum non-harming charging C-rate for a given 
	cell temperature according to Remmlinger et al. DOI: 10.1016/j.jpowsour.2013.12.101
	The function corresponds to a curve fit of the relationship shown in figure 7

  Input: 
	- T (Float64): temperature in °C

	Output: 
	- kder (Float64): derating factor in pu
	=#

	Crate_max = 5.838e-07T^4 + 2.641e-05T^3 + 0.0009767T^2 + 0.02488T + 0.2645;
	kder = clamp(Crate_max, 0.08,1.0)
	return kder
end

function ocv(soc)

	#=
	[Obtained from SimSES, Naumann et al., SimSES: Software for techno-economic Simulation of Stationary Energy Storage Systems]
	#https://ieeexplore.ieee.org/abstract/document/8278770

	BSD 3-Clause License

	Copyright (c) 2020, Team SES | Institute for Electrical Energy Storage
	Technology | Technical University of Munich
	All rights reserved.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions are met:

	* Redistributions of source code must retain the above copyright notice, this
		list of conditions and the following disclaimer.

	* Redistributions in binary form must reproduce the above copyright notice,
		this list of conditions and the following disclaimer in the documentation
		and/or other materials provided with the distribution.

	* Neither the name of the copyright holder nor the names of its
		contributors may be used to endorse or promote products derived from
		this software without specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
	DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
	FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
	DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
	SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
	CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
	OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
	OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
	=#

		a1 = -9.787
		a2 = 566
		a3 = 314.3
		a4 = 43.24
		b1 = -0.01471
		b2 = 0
		k0 = 3.4
		k1 = 0.2919
		k2 = -1.012
		k3 = -0.4422
		k4 = -0.913
		k5 = 0.09669

		ocv = k0 +
		k1 / (1 + exp(a1 * (soc - b1))) +
		k2 / (1 + exp(a2 * (soc - b2))) +
		k3 / (1 + exp(a3 * (soc - 1))) +
		k4 / (1 + exp(a4 * soc)) +
		k5 * soc
end

function Ri_interpolants()

	#=
	[Obtained from SimSES, Naumann et al., SimSES: Software for techno-economic Simulation of Stationary Energy Storage Systems]
	#https://ieeexplore.ieee.org/abstract/document/8278770

	BSD 3-Clause License

	Copyright (c) 2020, Team SES | Institute for Electrical Energy Storage
	Technology | Technical University of Munich
	All rights reserved.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions are met:

	* Redistributions of source code must retain the above copyright notice, this
		list of conditions and the following disclaimer.

	* Redistributions in binary form must reproduce the above copyright notice,
		this list of conditions and the following disclaimer in the documentation
		and/or other materials provided with the distribution.

	* Neither the name of the copyright holder nor the names of its
		contributors may be used to endorse or promote products derived from
		this software without specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
	DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
	FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
	DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
	SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
	CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
	OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
	OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
	=#


	SOC = 0:0.1:1
	T = [10, 25, 40, 60]
	Ri_ch =    [0.053553148 0.037583674 0.030476853 0.030476853;
							0.053553148 0.037583674 0.030476853 0.030476853;
							0.059217307 0.041399276 0.033085519 0.033085519;
							0.058450526 0.041399276 0.033149145 0.033149145;
							0.059280913 0.041779084 0.032640137 0.032640137;
							0.062648279 0.044767041 0.035246691 0.035246691;
							0.067545657 0.047695018 0.037666590 0.03766659;
							0.070789375 0.049790616 0.038939110 0.03893911;
							0.069199318 0.046804711 0.035566933 0.035566933;
							0.080652573 0.053421621 0.039954732 0.039954732;
							0.080652573 0.053421621 0.039954732 0.039954732]

	Ri_dch =   [0.194952789 0.061685557 0.044090174 0.044090174;
							0.194952789 0.061685557 0.044090174 0.044090174;
							0.080774939 0.050374992 0.038430102 0.038430102;
							0.075055233 0.048845551 0.037541588 0.037541588;
							0.074809788 0.051762822 0.040018594 0.040018594;
							0.069390124 0.047827935 0.038375676 0.038375676;
							0.065120954 0.044266280 0.035121551 0.035121551;
							0.061503438 0.041785842 0.032894641 0.032894641;
							0.064822427 0.046555916 0.037723435 0.037723435;
							0.061432462 0.042546505 0.033978319 0.033978319;
							0.061432462 0.042546505 0.033978319 0.033978319]

	Ri_ch_intp  = extrapolate(interpolate((SOC,T,),Ri_ch, Gridded(Linear())), Line())
	Ri_dch_intp = extrapolate(interpolate((SOC,T,),Ri_dch, Gridded(Linear())), Line())

	return Ri_ch_intp, Ri_dch_intp
end

end
