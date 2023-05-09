#########################################################################
##  
##  R code for simulating the operation of an alkaline hydrogen electrolyser (AEL)
##  
##  Written by Kate Ward, Imperial College London, 2019
##  
##  It needs 3 input files:
##    1.  EffCurveFile - Contains the variability of efficiency with power - effEL.csv
##    2.  WindElecFile - A wind power file - in this case half-hourly metered output from individual UK onshore windfarms for Q4 2016 - UKwindQ42016_norm.rds
##    3.  ThresholdFile - A file that gives the variability that can be expected sub-timescale - threshold_30_10_HB.csv
##  
##  The electrolyser is composed of n_stacks stacks, each of which may operate if supplied with power over the fraction EL_thresh of its rated power (P_EL).  
##  If the power falls below this then it is required to remain offline for a number of periods (NOff).  The auxilliary power is required whenever a stack 
##  of the electrolyser is operating - if it is not operational then no power is used.  The type of electrolyser is specified by ELtype.  Two sets of results 
##  are produced: one assumes that the efficiency of the electrolysis is independent of the electrolyser temperature.  The other assumes that the components 
##  are subject to Newtonian heating and cooling with an ambient temperature of T_amb and an electrolyser operating temperature of T_op - which is likely to 
##  be in the range 80 - 120 C.  The cooling/heating coefficient is k_Q.  Published output from field studies suggests that an AEL takes about 8 hours to go 
##  from ambient to optimal operational temperature, so a k_Q value of 1 is correct.  TRes is the temporal resolution of the input power data.
##  
##  



#########################################################################
#
#		INPUTS
#########################################################################

# Files
EffCurveFile = 'effEL.csv'
WindElecFile = 'UKwindQ42016_norm.rds'
ThresholdFile = 'threshold_30_0_HB.csv'		# use 60_0 if you have 1-hour resolution data


# Constants definition
p_aux = 0.075				# Fraction of rated power consumed by auxilliary processes
P_EL = 0.15					# Power of single electrolyser stack
n_stacks = 6				# number of stacks that comprise the electrolyser
EL_thresh = 0.2				# threshold below which the electrolyser cannot operate
NOff = 2					# the number of periods the electrolyser must stay off for
ELtype = 'AEL_Max'			# choice of electrolyser efficiency curve
T_amb = 14					# ambient temperature (degrees C)
T_op = 80					# stack operating temperature
k_Q = 0.5					# heating/cooling rate coefficient (per timestep)
TRes = 0.5					# time resolution of the power data (in hours)

set.seed(23)				# this is a stochastic model, so can enforce reproducible results



#########################################################################
#
#		FUNCTIONS
#########################################################################

right = function(text, num_char)
{
	substr(text, nchar(text) - (num_char-1), nchar(text))
}

Operability = function(Data, NOff, threshold) 
{
	# Determine the periods in which the electrolyser can operate
	# Inputs
	#	The raw minimum power data (Data)
	#	The number of off periods if threshold is exceeded (NOff)
	#	The threshold (threshold)
	#
	# Outputs
	#	Power that is available to the electrolyser (Power_in)

	PowerMin_temp = Data

	PowerAvail = array(0, dim = c(length(Data)))

	N = which(PowerMin_temp < threshold)[1]

	while (! is.na(N) ){
		PowerMin_temp[N:(N+NOff-1)] = 100
		N = which(PowerMin_temp < threshold)[1]
	}

	PowerAvail[which(PowerMin_temp <= 10)] = 1
	PowerAvail[which(PowerMin_temp == 0)] = 0
	PowerAvail[which(PowerMin_temp == 100)] = -1

	PowerAvail = PowerAvail[1:length(Data)]

	return(PowerAvail)
}

MinPower = function(Data, Threshold_lookup)
{
	# Takes the time averaged power input and creates a time series for whether 
	# the threshold was exceeded by the minimum value - stochastic

	N = length(Data)

	PowerOut = array(0, dim=c(N))
	RandVals = runif(N, min=0, max=1)
	PowerBin = 1 + floor(100*Data)


	for (i in 1:length(Data))
	{
		PowerOut[i] = (which(Threshold_lookup[(1+PowerBin[i]),] >= RandVals[i])[2] - 2) / 100
	}

	PowerOut[which(is.na(PowerOut))] = 1

	return(PowerOut)
}



#ThresholdFile = as.character(args[12])	
#WindElecFile = as.character(args[9])	
#EffCurveFile = as.character(args[11])	
#P_EL = as.numeric(args[1])									# Max power of each stack (MW)
#n_stacks = as.numeric(args[2])								# number of stacks
#P_aux = as.numeric(args[3])								# Fraction of max power that is used for auxiliary processes
#Threshold = as.numeric(args[4])							# operating threshold for stack as frac of max power (AEL only)
#T_amb = as.numeric(args[5])								# ambient temperature
#T_op = as.numeric(args[6])									# stack operating temperature
#k_Q = as.numeric(args[7])									# heating/cooling rate coefficient (per half hour)
#NOff = as.numeric(args[8])									# heating/cooling rate coefficient (per half hour)
#L_stack_output = FALSE										# do we want to produce output showing usage of individual stacks?
#TRes = as.numeric(args[13])
#ELtype = as.character(args[10])




#########################################################################
#
#		CALCULATIONS
#########################################################################

P_HG = (1+p_aux) * P_EL * n_stacks			# Hydrogen generator rated power
threshold = EL_thresh*P_EL

#  Read in data
Threshold_lookup_0 = read.csv(ThresholdFile)
Threshold_lookup = as.matrix(Threshold_lookup_0[ , -1])

# Read in efficiency curves
curve_names = read.csv(EffCurveFile, nrows=1)
curve = read.csv(EffCurveFile, skip=2)
colnames(curve) = colnames(curve_names)
UI = curve[ , ELtype]



# Read in the wind power file
PowerIn = readRDS(WindElecFile)

# dimension output array
OUT = data.frame(
	WindFarm  = character(),
	FarmPower = double(),
	P2EL      = double(),
	P4EL      = double(),
	P2H       = double(),
	P2Ht      = double()
)

P2H_Out = NULL

I = 1:n_stacks
J = 1:nrow(PowerIn)

for (i in I)
{
	name = paste0("Stack", i)
	OUT[[name]] = double()
}




# for our visualisations
layout( rbind(1:2, 3:4) )
par(mar = c(5, 4, 1.5, 1.5))
par(font.lab=2)

# run through each farm in turn
for (hh in 2:ncol(PowerIn))
{
	clear('Calculating farm', hh-1, 'of', ncol(PowerIn)-1, '...')

	name = colnames(PowerIn)[hh]
	Power = PowerIn[[name]]

	#Power = seq(0,1,0.05)		#Uncomment for short test sequence

	PowerMin = array(0, dim = c(length(Power)))
	PowerMin = MinPower(Power, Threshold_lookup)
	#PowerMin = Power 			#Uncomment to not use the variability of the winds


	# deduct the power that needs to be used for auxilliary processes
	Power2EL = Power - (p_aux*P_HG)
	Power2EL[ Power2EL < 0 ] = 0


	# deduct the power that needs to be used for auxilliary processes
	Power2ELMin = PowerMin - (p_aux*P_HG)
	Power2ELMin[ Power2ELMin < 0 ] = 0	

	Power2aux = array((p_aux*P_HG),dim = c(length(Power)))
	Power2aux[ Power2ELMin == 0 ] = 0


	# now need to loop over the stacks, allowing each to run independently of the others
	#
	# Power_stack 			Power available to each electrolyser stack
	# Power_stack_Min		Minimum power within each time period that is available to each stack
	# EL_power				Power that actually can be used by the electrolysers
	d = c(length(Power), n_stacks)

	Power_stack = array(0, dim=d)
	Power_stack_Min = array(0, dim=d)
	EL_power = array(0, dim=d)
	Efficiency = array(0, dim=d)
	P2H = array(0, dim=d)
	temperature = array(T_op, dim=d)

	N_stacks_temp = array(0, dim=length(Power))

	Power_temp = Power2EL
	PowerMin_frac = Power2ELMin/Power2EL
	PowerMin_frac[which(Power2EL == 0)] = 0
	availability = array(0,n_stacks)


	for (j in J)
	{
		Power_stack_temp = 0

		N_stacks = min(1 + floor(PowerMin_frac[j] * Power2EL[j] / threshold), sum(availability == 0))

		if (N_stacks != 0)
			Power_stack_temp = min(Power2EL[j] / N_stacks,P_EL)

		if (Power_stack_temp <= threshold)
			Power_stack_temp = 0

		# allocate power to stacks ...
		for (i in 1:N_stacks)
		{
			k = which(availability == 0)[i]
			Power_stack[j,k] = Power_stack_temp
		}

		# if we can't run electrolysers at all then no auxilliary power expended
		if (Power_stack_temp == 0)
			Power2aux[j] = 0     

		# ... and sort out availability & temperature
		if (j != 1)
		{
			for (i in I)
			{
				dT = ifelse(Power_stack[j,i] != 0,
						(T_op - temperature[j-1,i]) * exp(-k_Q),
						(T_amb - temperature[j-1,i]) * exp(-k_Q)
				)

				A = ifelse(Power_stack[j-1,i] != 0 & Power_stack[j,i] == 0,
						NOff - 1,
						max(0, availability[i] - 1)
				)

				temperature[j,i] = temperature[j-1,i] + dT
				availability[i] = A
			}
		}

	}




	EL_power = Power_stack

	i_EL_power = 1 + floor(100 * EL_power / P_EL)

	Efficiency[J, I] = UI[i_EL_power[J,I]]

	EfficiencyTemp = (0.0005 / 1.26) * (T_op - temperature)
	Efficiency2 = 1 / ((1/Efficiency) + EfficiencyTemp)

	P2H = EL_power * Efficiency / 100
	P2H_t = EL_power * Efficiency2 / 100

	P2losses = EL_power * (1 - Efficiency / 100)
	P2losses_t = EL_power * (1 - Efficiency2 / 100)

	P2Wind = (Power - Power2aux - rowSums(EL_power))

	util = apply(EL_power/P_EL, 2, mean)

	out = data.frame(
		WindFarm   = as.character(name),												# Windfarm name
		FarmEnergy = as.double(sum(TRes * Power)),										# Energy output from farm (MWh)
		P2aux      = as.double(sum(TRes * Power2aux)),									# Energy to auxiliary processes (MWh)
		P2wind     = as.double(sum(TRes * (Power - Power2aux - rowSums(EL_power)))),	# Energy supplied to grid from wind (MWh)
		P2H        = as.double(TRes * sum(P2H)),										# Energy stored as hydrogen (MWh) - neglecting temperature sensitivity
		P2losses   = as.double(TRes * sum(P2losses)),									# Energy losses from electrolysis (MWh) - neglecting temperature sensitivity
		P2Ht       = as.double(TRes * sum(P2H_t)),										# Energy stored as hydrogen (MWh) - temperature sensitivity included
		P2losses_t = as.double(TRes * sum(P2losses_t))									# Energy losses from electrolysis (MWh) - temperature sensitivity included
	)

	# Utilisation of each of the stacks (as a % of running flat out all the time)
	for (i in I)
	{
		name1 = paste0("StackUtil", i)
		out[[name1]] = as.double(util[i])
	}

	OUT = rbind(OUT, out)

	# visualise results
	x = OUT$FarmEnergy/nrow(PowerIn)/TRes
	y1 = OUT$P2Ht
	y2 = OUT$P2Ht/OUT$FarmEnergy
	y3 = OUT$P2wind/OUT$FarmEnergy
	hist(y1, col='grey90', main='', xlab='H2 production')
	plot(x, y1, pch=21, bg='grey90', xlim=c(0, max(x)*1.01), ylim=c(0, max(y1)*1.01), xlab='Farm CF', ylab='H2 production')
	plot(x, y2, pch=21, bg='grey90', xlim=c(0, max(x)*1.01), ylim=c(0, max(y2)*1.01), xlab='Farm CF', ylab='H2 conversion rate')
	plot(x, y3, pch=21, bg='grey90', xlim=c(0, max(x)*1.01), ylim=c(0, max(y3)*1.01), xlab='Farm CF', ylab='Electricity spilled')

}

#write.csv(OUT, 'OUT_AEC_ISKWIS.csv')
