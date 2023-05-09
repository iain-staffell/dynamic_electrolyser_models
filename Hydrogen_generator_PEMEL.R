#########################################################################
##  
##  R code for simulating the operation of an polymer electrolyte membrane hydrogen electrolyser (PEMEL)
##  
##  Written by Kate Ward, Imperial College London, 2019
##  
##  This is a deterministic model, but does account for the variability of power outputs over a 10 minute period through analysing the Haute-Borne wind farm data.  
##  This analysis is performed within the model as the choice of parameters affects what data is needed.
##  
##  It needs 3 input files:
##    1.  EffCurveFile - the electrolyser efficiency properties - same as for AEL - effEL.csv
##    2.  WindElecFile - UKwindQ42016_norm.rds - normalised half-hourly wind power outputs from UK onshore wind farms for Q4 2016
##    3.  ThresholdFile - the input Haute-Borne wind power output data - la-haute-borne-data-2012-2020_QC.rds.
##  
##  This model works a little differently to the AEL one.  The analysis of the data from La Haute Borne (or could be any other 10 minute resolution dataset) gives 
##  the likelihood that the power will slip below that required for operation of each of the stacks that comprises the electrolyser, given that the power supplied 
##  was at a particular level during a half-hour/hour (whatever the resolution is of the wind power output).  The PEM electrolyser doesn't require a period of time 
##  offline after stopping, therefore each stack runs with the mean power output for the half hour/hour period for a fraction of the time f, where f is the 
##  probability of the lower threshold being exceeded.  The exception is where the mean power is below the threshold: in this case there is a chance that there is 
##  sufficient power for some of the time, and the device will run with minimum power levels.
##  
##  The output goes to OUT_PEMEL.csv and has the same variables as for the AEL case.
##  
##  



#########################################################################
#
#		INPUTS
#########################################################################

# Files
EffCurveFile = 'effEL.csv'
WindElecFile = 'UKwindQ42016_norm.rds'
ThresholdFile = 'la-haute-borne-data-2012-2020_QC.rds'


# Constants definition
p_aux = 0.075				# Fraction of rated power consumed by auxilliary processes
P_EL = 0.2					# Power of single electrolyser stack
n_stacks = 5				# number of stacks that comprise the electrolyser
EL_thresh = 0.05			# threshold below which the electrolyser cannot operate
ELtype = 'PEMEL_Max'		# choice of electrolyser efficiency curve
T_amb = 14					# ambient temperature (degrees C)
T_op = 80					# stack operating temperature
k_Q = 0.5					# heating/cooling rate coefficient (per timestep)
TRes = 0.5					# time resolution of the power data (in hours)

nbins = 1000				# the desired resolution of the data analysis for the sub-timescale variability




#########################################################################
#
#		FUNCTIONS
#########################################################################

Threshold = function(thresh, p_aux, n_stacks, t_res, nbins)
{
	turbine = c('R80711', 'R80721', 'R80736', 'R80790')

	data = read_data(ThresholdFile)
	

	binsize = 1/nbins
	PowerBins = seq(-binsize, 1, binsize) + binsize/2
	PowerBins[1] = 0
	PowerBins[nbins+2] = 1

	Out_Avg= array(0, dim=c(n_stacks, nbins+2))
	Tot = array(0, nbins+2)

	for (k in 1:length(turbine))
	{

		MaxP = 2050
		name = paste0(turbine[k], "_power_avg")

		#	Normalise the power output and make sure all values are > 0

		PowerAvg10 = data[[name]]/MaxP
		PowerAvg10[which(PowerAvg10 < 0)] = 0

		n = length(PowerAvg10)
		Nstacks = array(0, n)

		#	Find the 30 minute mean power output by reshaping the 1D Power array and then taking the column means
		Power10min30 = matrix(PowerAvg10[1:(n-n%%t_res)], nrow = t_res)
		PowerAvg30 = as.vector(colMeans(Power10min30))

		#	Remove the auxilliary power from the 10 minute powers
		Power = PowerAvg10 - p_aux
		Power[which(Power < 0)] = 0


		for (i in 1:length(Power)) {

			Nstacks[i] = min(floor(Power[i]/thresh),n_stacks)
		}

		Nstacks30 = matrix(Nstacks[1:(n-n%%t_res)], nrow = t_res)

		# Find which bin the 30 minute mean is in so we know where to put the distribution

		# Sort into bins 
		iPowerAvg30 = (PowerAvg30 %/% binsize) + 2
		iPowerAvg30[which(PowerAvg30 == 0)] = 1
		iPowerAvg30[which(PowerAvg30 == 1)] = nbins+2

		for (i in 1:length(PowerAvg30))
		{
			ii = iPowerAvg30[i]
			Tot[ii]	= Tot[ii] + t_res

			for (jj in 1:t_res)
			{
				if (Nstacks30[jj,i] != 0 & !is.na(Nstacks30[jj,i]))
				{
					for (j in 1:Nstacks30[jj,i])
					{
						Out_Avg[j,ii]= Out_Avg[j,ii]+1
					}
				}
			}
		}

	}
	
	Out_Avg = sweep(Out_Avg, 2, FUN="/", STATS=Tot)

	return(Out_Avg)
}





#########################################################################
#
#		CALCULATIONS
#########################################################################

#File0 = as.character(args[12])	
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
#L_stack_output=FALSE										# do we want to produce output showing usage of individual stacks?
#TRes = as.numeric(args[13])
#ELtype = as.character(args[10])


binsize = 1 / nbins

P_HG = (1+p_aux) * P_EL * n_stacks			# Hydrogen generator rated power
threshold = EL_thresh * P_EL

#  Create data for the likelihood that power falls below the threshold for operation
clear('Reading thresholds...')
Threshold_lookup = Threshold(threshold, (P_HG*p_aux), n_stacks, TRes*6, nbins)
clear()


# Read in efficiency curves
curve_names = read.csv(EffCurveFile, nrows=1)
curve = read.csv(EffCurveFile, skip=2)
colnames(curve) = colnames(curve_names)
UI = curve[ , ELtype]



# Read in the wind power file
PowerIn = readRDS(WindElecFile)

# dimension output array
OUT = data.frame(WindFarm = character(),
				FarmEnergy = double(),
				P2aux = double(),
				P2wind = double(),
				P2H = double(),
				P2losses = double(),
				P2Ht = double(),
				P2losses_t = double()
)

for (i in 1:n_stacks)
{
	name2 = paste0("StackUtil", i)
	OUT[[name2]] = double()
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

	d = c(length(Power), n_stacks)
	Power_stack = array(1, dim=d)
	i_EL_power = array(0, dim=d)
	efficiency = array(0, dim=d)
	temperature = array(T_op, dim=d)
	P2H = array(0, dim=d)
	P2H_t = array(0, dim=d)
	EL_power = array(0, dim=d)

	# deduct the power that needs to be used for auxilliary processes
	Power2EL = Power - (p_aux*P_HG)
	Power2EL[which(Power2EL < 0)] = 0

	iPower = (Power %/% binsize) + 2
	iPower[which(Power == 0)] = 1
	iPower[which(Power == 1)] = nbins+2


	Power_stack = sweep(Power_stack, 1, FUN="*", Power2EL/n_stacks)

	Power_stack[which(Power_stack > P_EL)] = P_EL
	# set power values less than threshold to threshold - this will be corrected by Threshold_lookup
	Power_stack[which(Power_stack < threshold)]= threshold


	i_EL_power = 1 + floor(100 * Power_stack / P_EL)

	for (j in 1:length(Power))
	{
		for (i in 1:n_stacks)
		{
			efficiency[j,i] = UI[i_EL_power[j]]

			if (j != 1)
			{
				f = Threshold_lookup[i, iPower[j]]

				#  Runs first for f of the time period then is off for the remaining fraction of the time period
				temperature1_temp = T_op - (T_op - temperature[j-1, 1]) * exp(-f * k_Q)
				temperature1_temp = T_amb - (T_amb - temperature1_temp) * exp(-(1-f) * k_Q)

				# Off for the first fraction of the time period and then runs
				temperature2_temp = T_amb - (T_amb - temperature[j-1, i]) * exp(-(1-f) * k_Q)
				temperature2_temp = T_op - (T_op - temperature2_temp) * exp(-f * k_Q)

				# Say actual temperature is the average of the two
				temperature[j,i] = mean(temperature1_temp, temperature2_temp)
			}
		}
	}

	efficiencyTemp = (0.0005 / 1.26) * (T_op - temperature)
	efficiency2 = 1 / ((1/efficiency) + efficiencyTemp)


	if (n_stacks > 1)
	{
		Power2aux = p_aux * P_HG * apply(Threshold_lookup[ , iPower], 2, max)
		for (j in 1: length(Power))
		{
			EL_power[j, ] = Threshold_lookup[ , iPower[j]] * Power_stack[j]
			P2H[j, ] = EL_power[j, ] * t(efficiency[j, ]) / 100
			P2H_t[j, ] = EL_power[j, ] * t(efficiency2[j, ]) / 100
		}

	} else {

		Power2aux = p_aux * P_HG * Threshold_lookup[iPower]
		for (j in 1:length(Power))
		{
			EL_power[j, ] = Threshold_lookup[iPower[j]] * Power_stack[j]
			P2H[j,] = EL_power[j, ] * t(efficiency[j, ]) / 100
			P2H_t[j,] = EL_power[j, ] * t(efficiency2[j, ]) / 100
		}
	}

	P2losses = EL_power - P2H
	P2losses_t = EL_power - P2H_t

	EL_power_tot = rowSums(EL_power)	
	WindOut = Power - Power2aux - EL_power_tot
	Util = apply(EL_power / P_EL, 2, mean)

	out = data.frame(
		WindFarm   = as.character(name),
		FarmEnergy = as.double(sum(TRes * Power)),
		P2aux      = as.double(sum(TRes * Power2aux)),
		P2wind     = as.double(sum(TRes * WindOut)),
		P2H        = as.double(TRes * sum(P2H)),
		P2losses   = as.double(TRes * sum(P2losses)),
		P2Ht       = as.double(TRes * sum(P2H_t)),
		P2losses_t = as.double(TRes * sum(P2losses_t))
	)

	for (i in 1:n_stacks)
	{
		name1 = paste0("StackUtil", i)
		out[[name1]] = as.double(Util[i])
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


flush('\nWriting results...\n')
write.csv(OUT, 'OUT_PEMEL.csv')
