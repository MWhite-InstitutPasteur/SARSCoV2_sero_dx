
load("C:\\U\\CoronaVirus\\SeroSigPaper\\ANALYSIS\\SeroData\\Data_Process.RData")

part_ID <- as.vector(unique(BICHAT$patient_id))

tt_plot <- seq(from=-10, to=365, by=1)
log2 <- log(2)

N_tt_plot <- length(tt_plot)


Spike_cut <- 0.0001242311
RBD_cut   <- 0.0002260824
S1_cut    <- 0.0007174702
S2_cut    <- 0.0004809369

#######################################################################
#######################################################################
##                                                                   ##
##  #### #   ## ####   #### ##   ## #### ####   ##  ##  ####  ##     ##
##   ##  ##  ## ## ##   ##  ##   ##  ##  ## ##  ##  ## ##  ## ##     ##
##   ##  ### ## ##  ##  ##   ## ##   ##  ##  ## ##  ## ###### ##     ##
##   ##  ## ### ## ##   ##    ###    ##  ## ##  ##  ## ##  ## ##     ## 
##  #### ##  ## ####   ####    #    #### ####    ####  ##  ## #####  ##
##                                                                   ## 
#######################################################################
#######################################################################

#######################################################################
## Spike data

ant_names <- "Spike_IPP"

Spike_data = read.table( paste("C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/1_Data/", ant_names, "_IgG_dil.txt", sep="") )

N_tt <- 15

N_part <- nrow(Spike_data)


AB_max = 0.2
tt_max = max( Spike_data[,3:(2+N_tt)] )


AB_min = 1.95e-5




######################################
######################################
####
####   MEDIUM


Spike_file <-  paste( "C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/2_IgG_mod/Output/", ant_names, "_IgG_dil_local_med.txt", sep="" )




MCMC_ind <- read.table( Spike_file )

MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]

N_posterior <- nrow(MCMC_ind)

Spike_med_mod <- array(NA, dim=c(N_part, N_posterior, N_tt_plot) )
Spike_med_quant <- array(NA, dim=c(N_part, 3, N_tt_plot) )

for(n in 1:N_part)
{
	###################################
	## Model prediction for participant i
	## Posterior projections
	
	for(i in 1:N_posterior)
	{
		Ab_0    = MCMC_ind[i,(n-1)*9+1]
		beta    = MCMC_ind[i,(n-1)*9+2]
		tau     = MCMC_ind[i,(n-1)*9+3]		
		t_delta = MCMC_ind[i,(n-1)*9+4]
		t_short = MCMC_ind[i,(n-1)*9+5]
		t_long  = MCMC_ind[i,(n-1)*9+6]
		t_IgG   = MCMC_ind[i,(n-1)*9+7]
		rho     = MCMC_ind[i,(n-1)*9+8]

		r_delta = log(2)/t_delta
		r_cs = log(2)/t_short
		r_cl = log(2)/t_long
		r_a  = log(2)/t_IgG

		part_1 = rho / ((r_cs - r_a)*(r_cs - r_delta))
		part_2 = (1 - rho) / ((r_cl - r_a)*(r_cl - r_delta))
		part_3 = (r_a + r_cs*(rho-1) - r_cl*rho) / ((r_cs - r_a)*(r_cl - r_a)*(r_a - r_delta))
		part_4 = - (r_delta + r_cs*(rho - 1) - r_cl*rho) / ((r_cs - r_delta)*(r_cl - r_delta)*(r_a - r_delta))

		##AB_tt = Ab_0*exp(-r_cl*tt_plot)
		AB_tt =  rep( Ab_0, length(tt_plot) )

		tt_temp = (tt_plot - tau)[which(tt_plot >= tau )]
		AB_tt[which(tt_plot >= tau)] = AB_tt[which(tt_plot >= tau)] +
                                                  beta*( part_1*exp(-r_cs*tt_temp) + part_2*exp(-r_cl*tt_temp) +  
                                                         part_3*exp(-r_a*tt_temp) + part_4*exp(-r_delta*tt_temp) ) 

		Spike_med_mod[n,i,] = AB_tt
	}

	for(j in 1:length(tt_plot))
	{
		Spike_med_quant[n,,j] <- quantile( Spike_med_mod[n,,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
	}
}



Spike_med_sens <- matrix(NA, nrow=N_posterior, ncol=N_tt_plot)

for(i in 1:N_posterior)
{
	Spike_med_sens[i,] <- colSums(Spike_med_mod[,i,] > Spike_cut)/N_part
}

Spike_med_sens_quant <- matrix(NA, nrow=3, ncol=N_tt_plot)
for(j in 1:N_tt_plot)
{
	Spike_med_sens_quant[,j] <- quantile( Spike_med_sens[,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
}


######################################
######################################
####
####   SHORT


Spike_file <-  paste( "C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/2_IgG_mod/Output/", ant_names, "_IgG_dil_local_short.txt", sep="" )




MCMC_ind <- read.table( Spike_file )

MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]

N_posterior <- nrow(MCMC_ind)

Spike_short_mod <- array(NA, dim=c(N_part, N_posterior, N_tt_plot) )
Spike_short_quant <- array(NA, dim=c(N_part, 3, N_tt_plot) )

for(n in 1:N_part)
{
	###################################
	## Model prediction for participant i
	## Posterior projections
	
	for(i in 1:N_posterior)
	{
		Ab_0    = MCMC_ind[i,(n-1)*9+1]
		beta    = MCMC_ind[i,(n-1)*9+2]
		tau     = MCMC_ind[i,(n-1)*9+3]		
		t_delta = MCMC_ind[i,(n-1)*9+4]
		t_short = MCMC_ind[i,(n-1)*9+5]
		t_long  = MCMC_ind[i,(n-1)*9+6]
		t_IgG   = MCMC_ind[i,(n-1)*9+7]
		rho     = MCMC_ind[i,(n-1)*9+8]

		r_delta = log(2)/t_delta
		r_cs = log(2)/t_short
		r_cl = log(2)/t_long
		r_a  = log(2)/t_IgG

		part_1 = rho / ((r_cs - r_a)*(r_cs - r_delta))
		part_2 = (1 - rho) / ((r_cl - r_a)*(r_cl - r_delta))
		part_3 = (r_a + r_cs*(rho-1) - r_cl*rho) / ((r_cs - r_a)*(r_cl - r_a)*(r_a - r_delta))
		part_4 = - (r_delta + r_cs*(rho - 1) - r_cl*rho) / ((r_cs - r_delta)*(r_cl - r_delta)*(r_a - r_delta))

		##AB_tt = Ab_0*exp(-r_cl*tt_plot)
		AB_tt =  rep( Ab_0, length(tt_plot) )

		tt_temp = (tt_plot - tau)[which(tt_plot >= tau )]
		AB_tt[which(tt_plot >= tau)] = AB_tt[which(tt_plot >= tau)] +
                                                  beta*( part_1*exp(-r_cs*tt_temp) + part_2*exp(-r_cl*tt_temp) +  
                                                         part_3*exp(-r_a*tt_temp) + part_4*exp(-r_delta*tt_temp) ) 

		Spike_short_mod[n,i,] = AB_tt
	}

	for(j in 1:length(tt_plot))
	{
		Spike_short_quant[n,,j] <- quantile( Spike_short_mod[n,,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
	}
}



Spike_short_sens <- matrix(NA, nrow=N_posterior, ncol=N_tt_plot)

for(i in 1:N_posterior)
{
	Spike_short_sens[i,] <- colSums(Spike_short_mod[,i,] > Spike_cut)/N_part
}

Spike_short_sens_quant <- matrix(NA, nrow=3, ncol=N_tt_plot)
for(j in 1:N_tt_plot)
{
	Spike_short_sens_quant[,j] <- quantile( Spike_short_sens[,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
}



######################################
######################################
####
####   LONG


Spike_file <-  paste( "C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/2_IgG_mod/Output/", ant_names, "_IgG_dil_local_long.txt", sep="" )




MCMC_ind <- read.table( Spike_file )

MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]

N_posterior <- nrow(MCMC_ind)

Spike_long_mod <- array(NA, dim=c(N_part, N_posterior, N_tt_plot) )
Spike_long_quant <- array(NA, dim=c(N_part, 3, N_tt_plot) )

for(n in 1:N_part)
{
	###################################
	## Model prediction for participant i
	## Posterior projections
	
	for(i in 1:N_posterior)
	{
		Ab_0    = MCMC_ind[i,(n-1)*9+1]
		beta    = MCMC_ind[i,(n-1)*9+2]
		tau     = MCMC_ind[i,(n-1)*9+3]		
		t_delta = MCMC_ind[i,(n-1)*9+4]
		t_short = MCMC_ind[i,(n-1)*9+5]
		t_long  = MCMC_ind[i,(n-1)*9+6]
		t_IgG   = MCMC_ind[i,(n-1)*9+7]
		rho     = MCMC_ind[i,(n-1)*9+8]

		r_delta = log(2)/t_delta
		r_cs = log(2)/t_short
		r_cl = log(2)/t_long
		r_a  = log(2)/t_IgG

		part_1 = rho / ((r_cs - r_a)*(r_cs - r_delta))
		part_2 = (1 - rho) / ((r_cl - r_a)*(r_cl - r_delta))
		part_3 = (r_a + r_cs*(rho-1) - r_cl*rho) / ((r_cs - r_a)*(r_cl - r_a)*(r_a - r_delta))
		part_4 = - (r_delta + r_cs*(rho - 1) - r_cl*rho) / ((r_cs - r_delta)*(r_cl - r_delta)*(r_a - r_delta))

		##AB_tt = Ab_0*exp(-r_cl*tt_plot)
		AB_tt =  rep( Ab_0, length(tt_plot) )

		tt_temp = (tt_plot - tau)[which(tt_plot >= tau )]
		AB_tt[which(tt_plot >= tau)] = AB_tt[which(tt_plot >= tau)] +
                                                  beta*( part_1*exp(-r_cs*tt_temp) + part_2*exp(-r_cl*tt_temp) +  
                                                         part_3*exp(-r_a*tt_temp) + part_4*exp(-r_delta*tt_temp) ) 

		Spike_long_mod[n,i,] = AB_tt
	}

	for(j in 1:length(tt_plot))
	{
		Spike_long_quant[n,,j] <- quantile( Spike_long_mod[n,,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
	}
}



Spike_long_sens <- matrix(NA, nrow=N_posterior, ncol=N_tt_plot)

for(i in 1:N_posterior)
{
	Spike_long_sens[i,] <- colSums(Spike_long_mod[,i,] > Spike_cut)/N_part
}

Spike_long_sens_quant <- matrix(NA, nrow=3, ncol=N_tt_plot)
for(j in 1:N_tt_plot)
{
	Spike_long_sens_quant[,j] <- quantile( Spike_long_sens[,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
}






#######################################################################
## RBD data

ant_names <- "RBD_IPP"

RBD_data = read.table( paste("C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/1_Data/", ant_names, "_IgG_dil.txt", sep="") )

N_tt <- 15

N_part <- nrow(RBD_data)


AB_max = 0.2
tt_max = max( RBD_data[,3:(2+N_tt)] )


AB_min = 1.95e-5




######################################
######################################
####
####   MEDIUM


RBD_file <- paste( "C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/2_IgG_mod/Output/", ant_names, "_IgG_dil_local_med.txt", sep="" )




MCMC_ind <- read.table( RBD_file )

MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]

N_posterior <- nrow(MCMC_ind)

RBD_med_mod <- array(NA, dim=c(N_part, N_posterior, N_tt_plot) )
RBD_med_quant <- array(NA, dim=c(N_part, 3, N_tt_plot) )

for(n in 1:N_part)
{
	###################################
	## Model prediction for participant i
	## Posterior projections
	
	for(i in 1:N_posterior)
	{
		Ab_0    = MCMC_ind[i,(n-1)*9+1]
		beta    = MCMC_ind[i,(n-1)*9+2]
		tau     = MCMC_ind[i,(n-1)*9+3]		
		t_delta = MCMC_ind[i,(n-1)*9+4]
		t_short = MCMC_ind[i,(n-1)*9+5]
		t_long  = MCMC_ind[i,(n-1)*9+6]
		t_IgG   = MCMC_ind[i,(n-1)*9+7]
		rho     = MCMC_ind[i,(n-1)*9+8]

		r_delta = log(2)/t_delta
		r_cs = log(2)/t_short
		r_cl = log(2)/t_long
		r_a  = log(2)/t_IgG

		part_1 = rho / ((r_cs - r_a)*(r_cs - r_delta))
		part_2 = (1 - rho) / ((r_cl - r_a)*(r_cl - r_delta))
		part_3 = (r_a + r_cs*(rho-1) - r_cl*rho) / ((r_cs - r_a)*(r_cl - r_a)*(r_a - r_delta))
		part_4 = - (r_delta + r_cs*(rho - 1) - r_cl*rho) / ((r_cs - r_delta)*(r_cl - r_delta)*(r_a - r_delta))

		##AB_tt = Ab_0*exp(-r_cl*tt_plot)
		AB_tt =  rep( Ab_0, length(tt_plot) )

		tt_temp = (tt_plot - tau)[which(tt_plot >= tau )]
		AB_tt[which(tt_plot >= tau)] = AB_tt[which(tt_plot >= tau)] +
                                                  beta*( part_1*exp(-r_cs*tt_temp) + part_2*exp(-r_cl*tt_temp) +  
                                                         part_3*exp(-r_a*tt_temp) + part_4*exp(-r_delta*tt_temp) ) 

		RBD_med_mod[n,i,] = AB_tt
	}

	for(j in 1:length(tt_plot))
	{
		RBD_med_quant[n,,j] <- quantile( RBD_med_mod[n,,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
	}
}



RBD_med_sens <- matrix(NA, nrow=N_posterior, ncol=N_tt_plot)

for(i in 1:N_posterior)
{
	RBD_med_sens[i,] <- colSums(RBD_med_mod[,i,] > RBD_cut)/N_part
}

RBD_med_sens_quant <- matrix(NA, nrow=3, ncol=N_tt_plot)
for(j in 1:N_tt_plot)
{
	RBD_med_sens_quant[,j] <- quantile( RBD_med_sens[,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
}


######################################
######################################
####
####   SHORT


RBD_file <-  paste( "C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/2_IgG_mod/Output/", ant_names, "_IgG_dil_local_short.txt", sep="" )




MCMC_ind <- read.table( RBD_file )

MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]

N_posterior <- nrow(MCMC_ind)

RBD_short_mod <- array(NA, dim=c(N_part, N_posterior, N_tt_plot) )
RBD_short_quant <- array(NA, dim=c(N_part, 3, N_tt_plot) )

for(n in 1:N_part)
{
	###################################
	## Model prediction for participant i
	## Posterior projections
	
	for(i in 1:N_posterior)
	{
		Ab_0    = MCMC_ind[i,(n-1)*9+1]
		beta    = MCMC_ind[i,(n-1)*9+2]
		tau     = MCMC_ind[i,(n-1)*9+3]		
		t_delta = MCMC_ind[i,(n-1)*9+4]
		t_short = MCMC_ind[i,(n-1)*9+5]
		t_long  = MCMC_ind[i,(n-1)*9+6]
		t_IgG   = MCMC_ind[i,(n-1)*9+7]
		rho     = MCMC_ind[i,(n-1)*9+8]

		r_delta = log(2)/t_delta
		r_cs = log(2)/t_short
		r_cl = log(2)/t_long
		r_a  = log(2)/t_IgG

		part_1 = rho / ((r_cs - r_a)*(r_cs - r_delta))
		part_2 = (1 - rho) / ((r_cl - r_a)*(r_cl - r_delta))
		part_3 = (r_a + r_cs*(rho-1) - r_cl*rho) / ((r_cs - r_a)*(r_cl - r_a)*(r_a - r_delta))
		part_4 = - (r_delta + r_cs*(rho - 1) - r_cl*rho) / ((r_cs - r_delta)*(r_cl - r_delta)*(r_a - r_delta))

		##AB_tt = Ab_0*exp(-r_cl*tt_plot)
		AB_tt =  rep( Ab_0, length(tt_plot) )

		tt_temp = (tt_plot - tau)[which(tt_plot >= tau )]
		AB_tt[which(tt_plot >= tau)] = AB_tt[which(tt_plot >= tau)] +
                                                  beta*( part_1*exp(-r_cs*tt_temp) + part_2*exp(-r_cl*tt_temp) +  
                                                         part_3*exp(-r_a*tt_temp) + part_4*exp(-r_delta*tt_temp) ) 

		RBD_short_mod[n,i,] = AB_tt
	}

	for(j in 1:length(tt_plot))
	{
		RBD_short_quant[n,,j] <- quantile( RBD_short_mod[n,,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
	}
}



RBD_short_sens <- matrix(NA, nrow=N_posterior, ncol=N_tt_plot)

for(i in 1:N_posterior)
{
	RBD_short_sens[i,] <- colSums(RBD_short_mod[,i,] > RBD_cut)/N_part
}

RBD_short_sens_quant <- matrix(NA, nrow=3, ncol=N_tt_plot)
for(j in 1:N_tt_plot)
{
	RBD_short_sens_quant[,j] <- quantile( RBD_short_sens[,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
}



######################################
######################################
####
####   LONG


RBD_file <- paste( "C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/2_IgG_mod/Output/", ant_names, "_IgG_dil_local_long.txt", sep="" )




MCMC_ind <- read.table( RBD_file )

MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]

N_posterior <- nrow(MCMC_ind)

RBD_long_mod <- array(NA, dim=c(N_part, N_posterior, N_tt_plot) )
RBD_long_quant <- array(NA, dim=c(N_part, 3, N_tt_plot) )

for(n in 1:N_part)
{
	###################################
	## Model prediction for participant i
	## Posterior projections
	
	for(i in 1:N_posterior)
	{
		Ab_0    = MCMC_ind[i,(n-1)*9+1]
		beta    = MCMC_ind[i,(n-1)*9+2]
		tau     = MCMC_ind[i,(n-1)*9+3]		
		t_delta = MCMC_ind[i,(n-1)*9+4]
		t_short = MCMC_ind[i,(n-1)*9+5]
		t_long  = MCMC_ind[i,(n-1)*9+6]
		t_IgG   = MCMC_ind[i,(n-1)*9+7]
		rho     = MCMC_ind[i,(n-1)*9+8]

		r_delta = log(2)/t_delta
		r_cs = log(2)/t_short
		r_cl = log(2)/t_long
		r_a  = log(2)/t_IgG

		part_1 = rho / ((r_cs - r_a)*(r_cs - r_delta))
		part_2 = (1 - rho) / ((r_cl - r_a)*(r_cl - r_delta))
		part_3 = (r_a + r_cs*(rho-1) - r_cl*rho) / ((r_cs - r_a)*(r_cl - r_a)*(r_a - r_delta))
		part_4 = - (r_delta + r_cs*(rho - 1) - r_cl*rho) / ((r_cs - r_delta)*(r_cl - r_delta)*(r_a - r_delta))

		##AB_tt = Ab_0*exp(-r_cl*tt_plot)
		AB_tt =  rep( Ab_0, length(tt_plot) )

		tt_temp = (tt_plot - tau)[which(tt_plot >= tau )]
		AB_tt[which(tt_plot >= tau)] = AB_tt[which(tt_plot >= tau)] +
                                                  beta*( part_1*exp(-r_cs*tt_temp) + part_2*exp(-r_cl*tt_temp) +  
                                                         part_3*exp(-r_a*tt_temp) + part_4*exp(-r_delta*tt_temp) ) 

		RBD_long_mod[n,i,] = AB_tt
	}

	for(j in 1:length(tt_plot))
	{
		RBD_long_quant[n,,j] <- quantile( RBD_long_mod[n,,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
	}
}



RBD_long_sens <- matrix(NA, nrow=N_posterior, ncol=N_tt_plot)

for(i in 1:N_posterior)
{
	RBD_long_sens[i,] <- colSums(RBD_long_mod[,i,] > RBD_cut)/N_part
}

RBD_long_sens_quant <- matrix(NA, nrow=3, ncol=N_tt_plot)
for(j in 1:N_tt_plot)
{
	RBD_long_sens_quant[,j] <- quantile( RBD_long_sens[,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
}




#######################################################################
## S1 data

ant_names <- "S1_NA"

S1_data = read.table( paste("C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/1_Data/", ant_names, "_IgG_dil.txt", sep="") )

N_tt <- 15

N_part <- nrow(S1_data)


AB_max = 0.2
tt_max = max( S1_data[,3:(2+N_tt)] )


AB_min = 1.95e-5




######################################
######################################
####
####   MEDIUM


S1_file <-  paste( "C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/2_IgG_mod/Output/", ant_names, "_IgG_dil_local_med.txt", sep="" )




MCMC_ind <- read.table( S1_file )

MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]

N_posterior <- nrow(MCMC_ind)

S1_med_mod <- array(NA, dim=c(N_part, N_posterior, N_tt_plot) )
S1_med_quant <- array(NA, dim=c(N_part, 3, N_tt_plot) )

for(n in 1:N_part)
{
	###################################
	## Model prediction for participant i
	## Posterior projections
	
	for(i in 1:N_posterior)
	{
		Ab_0    = MCMC_ind[i,(n-1)*9+1]
		beta    = MCMC_ind[i,(n-1)*9+2]
		tau     = MCMC_ind[i,(n-1)*9+3]		
		t_delta = MCMC_ind[i,(n-1)*9+4]
		t_short = MCMC_ind[i,(n-1)*9+5]
		t_long  = MCMC_ind[i,(n-1)*9+6]
		t_IgG   = MCMC_ind[i,(n-1)*9+7]
		rho     = MCMC_ind[i,(n-1)*9+8]

		r_delta = log(2)/t_delta
		r_cs = log(2)/t_short
		r_cl = log(2)/t_long
		r_a  = log(2)/t_IgG

		part_1 = rho / ((r_cs - r_a)*(r_cs - r_delta))
		part_2 = (1 - rho) / ((r_cl - r_a)*(r_cl - r_delta))
		part_3 = (r_a + r_cs*(rho-1) - r_cl*rho) / ((r_cs - r_a)*(r_cl - r_a)*(r_a - r_delta))
		part_4 = - (r_delta + r_cs*(rho - 1) - r_cl*rho) / ((r_cs - r_delta)*(r_cl - r_delta)*(r_a - r_delta))

		##AB_tt = Ab_0*exp(-r_cl*tt_plot)
		AB_tt =  rep( Ab_0, length(tt_plot) )

		tt_temp = (tt_plot - tau)[which(tt_plot >= tau )]
		AB_tt[which(tt_plot >= tau)] = AB_tt[which(tt_plot >= tau)] +
                                                  beta*( part_1*exp(-r_cs*tt_temp) + part_2*exp(-r_cl*tt_temp) +  
                                                         part_3*exp(-r_a*tt_temp) + part_4*exp(-r_delta*tt_temp) ) 

		S1_med_mod[n,i,] = AB_tt
	}

	for(j in 1:length(tt_plot))
	{
		S1_med_quant[n,,j] <- quantile( S1_med_mod[n,,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
	}
}



S1_med_sens <- matrix(NA, nrow=N_posterior, ncol=N_tt_plot)

for(i in 1:N_posterior)
{
	S1_med_sens[i,] <- colSums(S1_med_mod[,i,] > S1_cut)/N_part
}

S1_med_sens_quant <- matrix(NA, nrow=3, ncol=N_tt_plot)
for(j in 1:N_tt_plot)
{
	S1_med_sens_quant[,j] <- quantile( S1_med_sens[,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
}


######################################
######################################
####
####   SHORT


S1_file <-  paste( "C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/2_IgG_mod/Output/", ant_names, "_IgG_dil_local_short.txt", sep="" )




MCMC_ind <- read.table( S1_file )

MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]

N_posterior <- nrow(MCMC_ind)

S1_short_mod <- array(NA, dim=c(N_part, N_posterior, N_tt_plot) )
S1_short_quant <- array(NA, dim=c(N_part, 3, N_tt_plot) )

for(n in 1:N_part)
{
	###################################
	## Model prediction for participant i
	## Posterior projections
	
	for(i in 1:N_posterior)
	{
		Ab_0    = MCMC_ind[i,(n-1)*9+1]
		beta    = MCMC_ind[i,(n-1)*9+2]
		tau     = MCMC_ind[i,(n-1)*9+3]		
		t_delta = MCMC_ind[i,(n-1)*9+4]
		t_short = MCMC_ind[i,(n-1)*9+5]
		t_long  = MCMC_ind[i,(n-1)*9+6]
		t_IgG   = MCMC_ind[i,(n-1)*9+7]
		rho     = MCMC_ind[i,(n-1)*9+8]

		r_delta = log(2)/t_delta
		r_cs = log(2)/t_short
		r_cl = log(2)/t_long
		r_a  = log(2)/t_IgG

		part_1 = rho / ((r_cs - r_a)*(r_cs - r_delta))
		part_2 = (1 - rho) / ((r_cl - r_a)*(r_cl - r_delta))
		part_3 = (r_a + r_cs*(rho-1) - r_cl*rho) / ((r_cs - r_a)*(r_cl - r_a)*(r_a - r_delta))
		part_4 = - (r_delta + r_cs*(rho - 1) - r_cl*rho) / ((r_cs - r_delta)*(r_cl - r_delta)*(r_a - r_delta))

		##AB_tt = Ab_0*exp(-r_cl*tt_plot)
		AB_tt =  rep( Ab_0, length(tt_plot) )

		tt_temp = (tt_plot - tau)[which(tt_plot >= tau )]
		AB_tt[which(tt_plot >= tau)] = AB_tt[which(tt_plot >= tau)] +
                                                  beta*( part_1*exp(-r_cs*tt_temp) + part_2*exp(-r_cl*tt_temp) +  
                                                         part_3*exp(-r_a*tt_temp) + part_4*exp(-r_delta*tt_temp) ) 

		S1_short_mod[n,i,] = AB_tt
	}

	for(j in 1:length(tt_plot))
	{
		S1_short_quant[n,,j] <- quantile( S1_short_mod[n,,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
	}
}



S1_short_sens <- matrix(NA, nrow=N_posterior, ncol=N_tt_plot)

for(i in 1:N_posterior)
{
	S1_short_sens[i,] <- colSums(S1_short_mod[,i,] > S1_cut)/N_part
}

S1_short_sens_quant <- matrix(NA, nrow=3, ncol=N_tt_plot)
for(j in 1:N_tt_plot)
{
	S1_short_sens_quant[,j] <- quantile( S1_short_sens[,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
}



######################################
######################################
####
####   LONG


S1_file <-  paste( "C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/2_IgG_mod/Output/", ant_names, "_IgG_dil_local_long.txt", sep="" )




MCMC_ind <- read.table( S1_file )

MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]

N_posterior <- nrow(MCMC_ind)

S1_long_mod <- array(NA, dim=c(N_part, N_posterior, N_tt_plot) )
S1_long_quant <- array(NA, dim=c(N_part, 3, N_tt_plot) )

for(n in 1:N_part)
{
	###################################
	## Model prediction for participant i
	## Posterior projections
	
	for(i in 1:N_posterior)
	{
		Ab_0    = MCMC_ind[i,(n-1)*9+1]
		beta    = MCMC_ind[i,(n-1)*9+2]
		tau     = MCMC_ind[i,(n-1)*9+3]		
		t_delta = MCMC_ind[i,(n-1)*9+4]
		t_short = MCMC_ind[i,(n-1)*9+5]
		t_long  = MCMC_ind[i,(n-1)*9+6]
		t_IgG   = MCMC_ind[i,(n-1)*9+7]
		rho     = MCMC_ind[i,(n-1)*9+8]

		r_delta = log(2)/t_delta
		r_cs = log(2)/t_short
		r_cl = log(2)/t_long
		r_a  = log(2)/t_IgG

		part_1 = rho / ((r_cs - r_a)*(r_cs - r_delta))
		part_2 = (1 - rho) / ((r_cl - r_a)*(r_cl - r_delta))
		part_3 = (r_a + r_cs*(rho-1) - r_cl*rho) / ((r_cs - r_a)*(r_cl - r_a)*(r_a - r_delta))
		part_4 = - (r_delta + r_cs*(rho - 1) - r_cl*rho) / ((r_cs - r_delta)*(r_cl - r_delta)*(r_a - r_delta))

		##AB_tt = Ab_0*exp(-r_cl*tt_plot)
		AB_tt =  rep( Ab_0, length(tt_plot) )

		tt_temp = (tt_plot - tau)[which(tt_plot >= tau )]
		AB_tt[which(tt_plot >= tau)] = AB_tt[which(tt_plot >= tau)] +
                                                  beta*( part_1*exp(-r_cs*tt_temp) + part_2*exp(-r_cl*tt_temp) +  
                                                         part_3*exp(-r_a*tt_temp) + part_4*exp(-r_delta*tt_temp) ) 

		S1_long_mod[n,i,] = AB_tt
	}

	for(j in 1:length(tt_plot))
	{
		S1_long_quant[n,,j] <- quantile( S1_long_mod[n,,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
	}
}



S1_long_sens <- matrix(NA, nrow=N_posterior, ncol=N_tt_plot)

for(i in 1:N_posterior)
{
	S1_long_sens[i,] <- colSums(S1_long_mod[,i,] > S1_cut)/N_part
}

S1_long_sens_quant <- matrix(NA, nrow=3, ncol=N_tt_plot)
for(j in 1:N_tt_plot)
{
	S1_long_sens_quant[,j] <- quantile( S1_long_sens[,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
}




#######################################################################
## S2 data

ant_names <- "S2_NA"

S2_data = read.table( paste("C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/1_Data/", ant_names, "_IgG_dil.txt", sep="") )

N_tt <- 15

N_part <- nrow(S2_data)


AB_max = 0.2
tt_max = max( S2_data[,3:(2+N_tt)] )


AB_min = 1.95e-5




######################################
######################################
####
####   MEDIUM


S2_file <-  paste( "C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/2_IgG_mod/Output/", ant_names, "_IgG_dil_local_med.txt", sep="" )




MCMC_ind <- read.table( S2_file )

MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]

N_posterior <- nrow(MCMC_ind)

S2_med_mod <- array(NA, dim=c(N_part, N_posterior, N_tt_plot) )
S2_med_quant <- array(NA, dim=c(N_part, 3, N_tt_plot) )

for(n in 1:N_part)
{
	###################################
	## Model prediction for participant i
	## Posterior projections
	
	for(i in 1:N_posterior)
	{
		Ab_0    = MCMC_ind[i,(n-1)*9+1]
		beta    = MCMC_ind[i,(n-1)*9+2]
		tau     = MCMC_ind[i,(n-1)*9+3]		
		t_delta = MCMC_ind[i,(n-1)*9+4]
		t_short = MCMC_ind[i,(n-1)*9+5]
		t_long  = MCMC_ind[i,(n-1)*9+6]
		t_IgG   = MCMC_ind[i,(n-1)*9+7]
		rho     = MCMC_ind[i,(n-1)*9+8]

		r_delta = log(2)/t_delta
		r_cs = log(2)/t_short
		r_cl = log(2)/t_long
		r_a  = log(2)/t_IgG

		part_1 = rho / ((r_cs - r_a)*(r_cs - r_delta))
		part_2 = (1 - rho) / ((r_cl - r_a)*(r_cl - r_delta))
		part_3 = (r_a + r_cs*(rho-1) - r_cl*rho) / ((r_cs - r_a)*(r_cl - r_a)*(r_a - r_delta))
		part_4 = - (r_delta + r_cs*(rho - 1) - r_cl*rho) / ((r_cs - r_delta)*(r_cl - r_delta)*(r_a - r_delta))

		##AB_tt = Ab_0*exp(-r_cl*tt_plot)
		AB_tt =  rep( Ab_0, length(tt_plot) )

		tt_temp = (tt_plot - tau)[which(tt_plot >= tau )]
		AB_tt[which(tt_plot >= tau)] = AB_tt[which(tt_plot >= tau)] +
                                                  beta*( part_1*exp(-r_cs*tt_temp) + part_2*exp(-r_cl*tt_temp) +  
                                                         part_3*exp(-r_a*tt_temp) + part_4*exp(-r_delta*tt_temp) ) 

		S2_med_mod[n,i,] = AB_tt
	}

	for(j in 1:length(tt_plot))
	{
		S2_med_quant[n,,j] <- quantile( S2_med_mod[n,,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
	}
}



S2_med_sens <- matrix(NA, nrow=N_posterior, ncol=N_tt_plot)

for(i in 1:N_posterior)
{
	S2_med_sens[i,] <- colSums(S2_med_mod[,i,] > S2_cut)/N_part
}

S2_med_sens_quant <- matrix(NA, nrow=3, ncol=N_tt_plot)
for(j in 1:N_tt_plot)
{
	S2_med_sens_quant[,j] <- quantile( S2_med_sens[,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
}


######################################
######################################
####
####   SHORT


S2_file <-  paste( "C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/2_IgG_mod/Output/", ant_names, "_IgG_dil_local_short.txt", sep="" )




MCMC_ind <- read.table( S2_file )

MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]

N_posterior <- nrow(MCMC_ind)

S2_short_mod <- array(NA, dim=c(N_part, N_posterior, N_tt_plot) )
S2_short_quant <- array(NA, dim=c(N_part, 3, N_tt_plot) )

for(n in 1:N_part)
{
	###################################
	## Model prediction for participant i
	## Posterior projections
	
	for(i in 1:N_posterior)
	{
		Ab_0    = MCMC_ind[i,(n-1)*9+1]
		beta    = MCMC_ind[i,(n-1)*9+2]
		tau     = MCMC_ind[i,(n-1)*9+3]		
		t_delta = MCMC_ind[i,(n-1)*9+4]
		t_short = MCMC_ind[i,(n-1)*9+5]
		t_long  = MCMC_ind[i,(n-1)*9+6]
		t_IgG   = MCMC_ind[i,(n-1)*9+7]
		rho     = MCMC_ind[i,(n-1)*9+8]

		r_delta = log(2)/t_delta
		r_cs = log(2)/t_short
		r_cl = log(2)/t_long
		r_a  = log(2)/t_IgG

		part_1 = rho / ((r_cs - r_a)*(r_cs - r_delta))
		part_2 = (1 - rho) / ((r_cl - r_a)*(r_cl - r_delta))
		part_3 = (r_a + r_cs*(rho-1) - r_cl*rho) / ((r_cs - r_a)*(r_cl - r_a)*(r_a - r_delta))
		part_4 = - (r_delta + r_cs*(rho - 1) - r_cl*rho) / ((r_cs - r_delta)*(r_cl - r_delta)*(r_a - r_delta))

		##AB_tt = Ab_0*exp(-r_cl*tt_plot)
		AB_tt =  rep( Ab_0, length(tt_plot) )

		tt_temp = (tt_plot - tau)[which(tt_plot >= tau )]
		AB_tt[which(tt_plot >= tau)] = AB_tt[which(tt_plot >= tau)] +
                                                  beta*( part_1*exp(-r_cs*tt_temp) + part_2*exp(-r_cl*tt_temp) +  
                                                         part_3*exp(-r_a*tt_temp) + part_4*exp(-r_delta*tt_temp) ) 

		S2_short_mod[n,i,] = AB_tt
	}

	for(j in 1:length(tt_plot))
	{
		S2_short_quant[n,,j] <- quantile( S2_short_mod[n,,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
	}
}



S2_short_sens <- matrix(NA, nrow=N_posterior, ncol=N_tt_plot)

for(i in 1:N_posterior)
{
	S2_short_sens[i,] <- colSums(S2_short_mod[,i,] > S2_cut)/N_part
}

S2_short_sens_quant <- matrix(NA, nrow=3, ncol=N_tt_plot)
for(j in 1:N_tt_plot)
{
	S2_short_sens_quant[,j] <- quantile( S2_short_sens[,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
}



######################################
######################################
####
####   LONG


S2_file <-  paste( "C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/2_IgG_mod/Output/", ant_names, "_IgG_dil_local_long.txt", sep="" )




MCMC_ind <- read.table( S2_file )

MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]

N_posterior <- nrow(MCMC_ind)

S2_long_mod <- array(NA, dim=c(N_part, N_posterior, N_tt_plot) )
S2_long_quant <- array(NA, dim=c(N_part, 3, N_tt_plot) )

for(n in 1:N_part)
{
	###################################
	## Model prediction for participant i
	## Posterior projections
	
	for(i in 1:N_posterior)
	{
		Ab_0    = MCMC_ind[i,(n-1)*9+1]
		beta    = MCMC_ind[i,(n-1)*9+2]
		tau     = MCMC_ind[i,(n-1)*9+3]		
		t_delta = MCMC_ind[i,(n-1)*9+4]
		t_short = MCMC_ind[i,(n-1)*9+5]
		t_long  = MCMC_ind[i,(n-1)*9+6]
		t_IgG   = MCMC_ind[i,(n-1)*9+7]
		rho     = MCMC_ind[i,(n-1)*9+8]

		r_delta = log(2)/t_delta
		r_cs = log(2)/t_short
		r_cl = log(2)/t_long
		r_a  = log(2)/t_IgG

		part_1 = rho / ((r_cs - r_a)*(r_cs - r_delta))
		part_2 = (1 - rho) / ((r_cl - r_a)*(r_cl - r_delta))
		part_3 = (r_a + r_cs*(rho-1) - r_cl*rho) / ((r_cs - r_a)*(r_cl - r_a)*(r_a - r_delta))
		part_4 = - (r_delta + r_cs*(rho - 1) - r_cl*rho) / ((r_cs - r_delta)*(r_cl - r_delta)*(r_a - r_delta))

		##AB_tt = Ab_0*exp(-r_cl*tt_plot)
		AB_tt =  rep( Ab_0, length(tt_plot) )

		tt_temp = (tt_plot - tau)[which(tt_plot >= tau )]
		AB_tt[which(tt_plot >= tau)] = AB_tt[which(tt_plot >= tau)] +
                                                  beta*( part_1*exp(-r_cs*tt_temp) + part_2*exp(-r_cl*tt_temp) +  
                                                         part_3*exp(-r_a*tt_temp) + part_4*exp(-r_delta*tt_temp) ) 

		S2_long_mod[n,i,] = AB_tt
	}

	for(j in 1:length(tt_plot))
	{
		S2_long_quant[n,,j] <- quantile( S2_long_mod[n,,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
	}
}



S2_long_sens <- matrix(NA, nrow=N_posterior, ncol=N_tt_plot)

for(i in 1:N_posterior)
{
	S2_long_sens[i,] <- colSums(S2_long_mod[,i,] > S2_cut)/N_part
}

S2_long_sens_quant <- matrix(NA, nrow=3, ncol=N_tt_plot)
for(j in 1:N_tt_plot)
{
	S2_long_sens_quant[,j] <- quantile( S2_long_sens[,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
}



























##################################
##################################
##                              ##  
##  #####  ##     ####  ######  ##
##  ##  ## ##    ##  ##   ##    ##
##  #####  ##    ##  ##   ##    ##
##  ##     ##    ##  ##   ##    ## 
##  ##     #####  ####    ##    ##
##                              ##
##################################
##################################


tiff(file="Fig4_long_kinetics_sensitivity.tif", width=40, height=25, units="cm", res=500)

lay.mat <- rbind( c( 1, 2, 3, 4 ),
                  c( 5, 6, 7, 8 ) )
layout(lay.mat, heights=c(10,10), widths=c(1,1,1,1))
layout.show(8)

par(mar=c(4,6.5,3,1))
par(mgp=c(2.5, 1, 0))

point.size = 0.75
lab.size   = 1.5
axis.size  = 1.5
main.size  = 2
line.size  = 2




#############
## Panel 1: Spike

line_seq_x <- c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03)
line_seq_y <- 30*( (-1):12 )

plot(x = 100, y = 100,
     col="orangered", pch=19,
     xlim=c(-30,360), ylim=c(1e-5,0.02), log="y",
     yaxt='n', xaxt='n', bty='n',
     ylab="", xlab="months after symptom onset",
     main=expression(paste( "(A) modelled anti-S"^"tri", " IgG kinetics", sep="" )),
     cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_x))
{
	points(x=c(-100,1000), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
}
for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
}

points(x=c(-100,1000), y=rep(1.95e-5,2), type='l', col="black", lty="dashed")
points(x=c(-100,1000), y=rep(0.02,2), type='l', col="black", lty="dashed")

for(n in 1:N_part)
{
	points(x=tt_plot, y=Spike_med_quant[n,2,], type='l', col="grey")
}

points(x = BICHAT$days_post, y = BICHAT$Spike_IPP_IgG_dil,
     col="orangered", pch=19)

points(x = runif(n=length(COCHIN$Spike_IPP_IgG_dil), min=8, max=12), 
       y = COCHIN$Spike_IPP_IgG_dil,
       col="orangered", pch=19)

points(x = runif(n=length(STRAS$Spike_IPP_IgG_dil), min=8, max=12), 
       y = STRAS$Spike_IPP_IgG_dil,
       col="orangered", pch=19)


points(x = runif(n=length(TRC$Spike_IPP_IgG_dil), min=-25, max=-1), 
       y = TRC$Spike_IPP_IgG_dil,
       col="orangered", pch=4)


points(x = runif(n=length(PNC$Spike_IPP_IgG_dil), min=-25, max=-1), 
       y = PNC$Spike_IPP_IgG_dil,
       col="orangered", pch=4)

points(x = runif(n=length(EFS$Spike_IPP_IgG_dil), min=-25, max=-1), 
       y = EFS$Spike_IPP_IgG_dil,
       col="orangered", pch=4)

points(x=c(-100,1000), y=rep(Spike_cut,2), type='l', col="orangered", lty="dashed", lwd=2)


mtext(side = 2, line = 4, 
cex=0.8*axis.size, 
text="antibody dilution")



axis(1, at=30*c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), 
        label=c("", "0", "", "", "3", "", "", "6", "", "", "9", "", "", "12"), 
        cex.axis=0.8*axis.size )

axis(2, at = c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03), 
          label = c("0.00001", "", "0.0001", "", "0.001", "", "0.01", ""),
          las=2, cex.axis=0.8*axis.size )




#############
## Panel 2: RBD

line_seq_x <- c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03)
line_seq_y <- 30*( (-1):12 )

plot(x = 100, y = 100,
     col="gold", pch=19,
     xlim=c(-30,360), ylim=c(1e-5,0.02), log="y",
     yaxt='n', xaxt='n', bty='n',
     ylab="", xlab="months after symptom onset",
     main=expression(paste( "(B) modelled anti-RBD IgG kinetics", sep="" )),
     cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_x))
{
	points(x=c(-100,1000), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
}
for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
}

points(x=c(-100,1000), y=rep(1.95e-5,2), type='l', col="black", lty="dashed")
points(x=c(-100,1000), y=rep(0.02,2), type='l', col="black", lty="dashed")

for(n in 1:N_part)
{
	points(x=tt_plot, y=RBD_med_quant[n,2,], type='l', col="grey")
}

points(x = BICHAT$days_post, y = BICHAT$RBD_IPP_IgG_dil,
     col="gold", pch=19)

points(x = runif(n=length(COCHIN$RBD_IPP_IgG_dil), min=8, max=12), 
       y = COCHIN$Spike_IPP_IgG_dil,
       col="gold", pch=19)

points(x = runif(n=length(STRAS$RBD_IPP_IgG_dil), min=8, max=12), 
       y = STRAS$Spike_IPP_IgG_dil,
       col="gold", pch=19)


points(x = runif(n=length(TRC$RBD_IPP_IgG_dil), min=-25, max=-1), 
       y = TRC$RBD_IPP_IgG_dil,
       col="gold", pch=4)


points(x = runif(n=length(PNC$RBD_IPP_IgG_dil), min=-25, max=-1), 
       y = PNC$RBD_IPP_IgG_dil,
       col="gold", pch=4)

points(x = runif(n=length(EFS$RBD_IPP_IgG_dil), min=-25, max=-1), 
       y = EFS$RBD_IPP_IgG_dil,
       col="gold", pch=4)


points(x=c(-100,1000), y=rep(RBD_cut,2), type='l', col="gold", lty="dashed", lwd=2)


mtext(side = 2, line = 4, 
cex=0.8*axis.size, 
text="antibody dilution")



axis(1, at=30*c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), 
        label=c("", "0", "", "", "3", "", "", "6", "", "", "9", "", "", "12"), 
        cex.axis=0.8*axis.size )

axis(2, at = c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03), 
          label = c("0.00001", "", "0.0001", "", "0.001", "", "0.01", ""),
          las=2, cex.axis=0.8*axis.size )






#############
## Panel 3: S1

line_seq_x <- c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03)
line_seq_y <- 30*( (-1):12 )

plot(x = 100, y = 100,
     col="gold", pch=19,
     xlim=c(-30,360), ylim=c(1e-5,0.02), log="y",
     yaxt='n', xaxt='n', bty='n',
     ylab="", xlab="months after symptom onset",
     main=expression(paste( "(C) modelled anti-S1 IgG kinetics", sep="" )),
     cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_x))
{
	points(x=c(-100,1000), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
}
for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
}

points(x=c(-100,1000), y=rep(1.95e-5,2), type='l', col="black", lty="dashed")
points(x=c(-100,1000), y=rep(0.02,2), type='l', col="black", lty="dashed")

for(n in 1:N_part)
{
	points(x=tt_plot, y=S1_med_quant[n,2,], type='l', col="grey")
}

points(x = BICHAT$days_post, y = BICHAT$S1_NA_IgG_dil,
     col="limegreen", pch=19)

points(x = runif(n=length(COCHIN$S1_NA_IgG_dil), min=8, max=12), 
       y = COCHIN$S1_NA_IgG_dil,
       col="limegreen", pch=19)

points(x = runif(n=length(STRAS$S1_NA_IgG_dil), min=8, max=12), 
       y = STRAS$S1_NA_IgG_dil,
       col="limegreen", pch=19)


points(x = runif(n=length(TRC$S1_NA_IgG_dil), min=-25, max=-1), 
       y = TRC$S1_NA_IgG_dil,
       col="limegreen", pch=4)


points(x = runif(n=length(PNC$S1_NA_IgG_dil), min=-25, max=-1), 
       y = PNC$S1_NA_IgG_dil,
       col="limegreen", pch=4)


points(x = runif(n=length(EFS$S1_NA_IgG_dil), min=-25, max=-1), 
       y = EFS$S1_NA_IgG_dil,
       col="limegreen", pch=4)

points(x=c(-100,1000), y=rep(S1_cut,2), type='l', col="limegreen", lty="dashed", lwd=2)



mtext(side = 2, line = 4, 
cex=0.8*axis.size, 
text="antibody dilution")



axis(1, at=30*c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), 
        label=c("", "0", "", "", "3", "", "", "6", "", "", "9", "", "", "12"), 
        cex.axis=0.8*axis.size )

axis(2, at = c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03), 
          label = c("0.00001", "", "0.0001", "", "0.001", "", "0.01", ""),
          las=2, cex.axis=0.8*axis.size )




#############
## Panel 4: S2

line_seq_x <- c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03)
line_seq_y <- 30*( (-1):12 )

plot(x = 100, y = 100,
     col="gold", pch=19,
     xlim=c(-30,360), ylim=c(1e-5,0.02), log="y",
     yaxt='n', xaxt='n', bty='n',
     ylab="", xlab="months after symptom onset",
     main=expression(paste( "(D) modelled anti-S2 IgG kinetics", sep="" )),
     cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_x))
{
	points(x=c(-100,1000), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
}
for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
}

points(x=c(-100,1000), y=rep(1.95e-5,2), type='l', col="black", lty="dashed")
points(x=c(-100,1000), y=rep(0.02,2), type='l', col="black", lty="dashed")

for(n in 1:N_part)
{
	points(x=tt_plot, y=S2_med_quant[n,2,], type='l', col="grey")
}

points(x = BICHAT$days_post, y = BICHAT$S2_NA_IgG_dil,
     col="mediumblue", pch=19)

points(x = runif(n=length(COCHIN$S2_NA_IgG_dil), min=8, max=12), 
       y = COCHIN$S2_NA_IgG_dil,
       col="mediumblue", pch=19)

points(x = runif(n=length(STRAS$S2_NA_IgG_dil), min=8, max=12), 
       y = STRAS$S2_NA_IgG_dil,
       col="mediumblue", pch=19)


points(x = runif(n=length(TRC$S2_NA_IgG_dil), min=-25, max=-1), 
       y = TRC$S2_NA_IgG_dil,
       col="mediumblue", pch=4)


points(x = runif(n=length(PNC$S2_NA_IgG_dil), min=-25, max=-1), 
       y = PNC$S2_NA_IgG_dil,
       col="mediumblue", pch=4)

points(x = runif(n=length(EFS$S2_NA_IgG_dil), min=-25, max=-1), 
       y = EFS$S2_NA_IgG_dil,
       col="mediumblue", pch=4)

points(x=c(-100,1000), y=rep(S2_cut,2), type='l', col="mediumblue", lty="dashed", lwd=2)



mtext(side = 2, line = 4, 
cex=0.8*axis.size, 
text="antibody dilution")


axis(1, at=30*c(-1, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), 
        label=c("", "0", "", "", "3", "", "", "6", "", "", "9", "", "", "12"), 
        cex.axis=0.8*axis.size )

axis(2, at = c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03), 
          label = c("0.00001", "", "0.0001", "", "0.001", "", "0.01", ""),
          las=2, cex.axis=0.8*axis.size )





##################################
##                              ##
##  Panel 5: Spike sensitivity  ##
##                              ##
##################################

line_seq_x <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1)

plot( x=100, y=100, pch=19, cex=point.size, col="firebrick1", 
      xlim=c(0,365), ylim=c(0,1), 
	xlab="months after symptom onset", ylab="",
      main=expression(paste( "(E) anti-S"^"tri", " IgG sensitivity", sep="" )),
      xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
      cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size   )


axis(1, at=30*c( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), 
        labels=c("0", "", "", "3", "", "", "6", "", "", "9", "", "", "12"), 
        cex.axis=axis.size)

axis(1, at=30*c(0, 3,  6, 9,  12), 
        labels=c("0", "3", "6",  "9",  "12"), 
        cex.axis=axis.size)

axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        las=2, cex.axis=axis.size)


mtext(side = 2, line = 4, 
	cex=0.8*axis.size, 
	text="  sensitivity")

for(i in 1:length(line_seq_x))
{
	points(x=c(-60,1e4), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}





polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( Spike_med_sens_quant[1,], rev(Spike_med_sens_quant[3,]) ),
        col=rgb(190/256,190/256,190/256,0.4), border=NA)

points(x=tt_plot, y=Spike_med_sens_quant[2,], 
type='l', col="orangered", lwd=line.size )


points(x=tt_plot, y=Spike_long_sens_quant[2,], 
type='l', lty="longdash", col="orangered", lwd=line.size )

points(x=tt_plot, y=Spike_short_sens_quant[2,], 
type='l', lty="dotted", col="orangered", lwd=line.size )


##################################
##                              ##
##  Panel 6: RBD sensitivity    ##
##                              ##
##################################

plot( x=100, y=100, pch=19, cex=point.size, col="firebrick1", 
      xlim=c(0,365), ylim=c(0,1), 
	xlab="months after symptom onset", ylab="",
      main=expression(paste( "(F) anti-RBD IgG sensitivity", sep="" )),
      xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
      cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size   )


axis(1, at=30*c( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), 
        labels=c("0", "", "", "3", "", "", "6", "", "", "9", "", "", "12"), 
        cex.axis=axis.size)

axis(1, at=30*c(0, 3,  6, 9,  12), 
        labels=c("0", "3", "6",  "9",  "12"), 
        cex.axis=axis.size)

axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        las=2, cex.axis=axis.size)


mtext(side = 2, line = 4, 
	cex=0.8*axis.size, 
	text="  sensitivity")

for(i in 1:length(line_seq_x))
{
	points(x=c(-60,1e4), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}





polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( RBD_med_sens_quant[1,], rev(RBD_med_sens_quant[3,]) ),
        col=rgb(190/256,190/256,190/256,0.4), border=NA)

points(x=tt_plot, y=RBD_med_sens_quant[2,], 
type='l', col="gold", lwd=line.size )

points(x=tt_plot, y=RBD_long_sens_quant[2,], 
type='l', lty="longdash", col="gold", lwd=line.size )

points(x=tt_plot, y=RBD_short_sens_quant[2,], 
type='l', lty="dotted", col="gold", lwd=line.size )


##################################
##                              ##
##  Panel 7: S1 sensitivity     ##
##                              ##
##################################

plot( x=100, y=100, pch=19, cex=point.size, col="firebrick1", 
      xlim=c(0,365), ylim=c(0,1), 
	xlab="months after symptom onset", ylab="",
      main=expression(paste( "(G) anti-S1 IgG sensitivity", sep="" )),
      xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
      cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size   )


axis(1, at=30*c( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), 
        labels=c("0", "", "", "3", "", "", "6", "", "", "9", "", "", "12"), 
        cex.axis=axis.size)

axis(1, at=30*c(0, 3,  6, 9,  12), 
        labels=c("0", "3", "6",  "9",  "12"), 
        cex.axis=axis.size)

axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        las=2, cex.axis=axis.size)


mtext(side = 2, line = 4, 
	cex=0.8*axis.size, 
	text="  sensitivity")

for(i in 1:length(line_seq_x))
{
	points(x=c(-60,1e4), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}





polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( S1_med_sens_quant[1,], rev(S1_med_sens_quant[3,]) ),
        col=rgb(190/256,190/256,190/256,0.4), border=NA)

points(x=tt_plot, y=S1_med_sens_quant[2,], 
type='l', col="limegreen", lwd=line.size )


points(x=tt_plot, y=S1_long_sens_quant[2,], 
type='l', lty="longdash", col="limegreen", lwd=line.size )


points(x=tt_plot, y=S1_short_sens_quant[2,], 
type='l', lty="dotted", col="limegreen", lwd=line.size )



##################################
##                              ##
##  Panel 8: S2 sensitivity     ##
##                              ##
##################################

plot( x=100, y=100, pch=19, cex=point.size, col="firebrick1", 
      xlim=c(0,365), ylim=c(0,1), 
	xlab="months after symptom onset", ylab="",
      main=expression(paste( "(H) anti-S2 IgG sensitivity", sep="" )),
      xaxs="i", yaxs="i", xaxt="n", yaxt="n", bty='n',
      cex.main=main.size, cex.axis=axis.size, cex.lab=lab.size   )


axis(1, at=30*c( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), 
        labels=c("0", "", "", "3", "", "", "6", "", "", "9", "", "", "12"), 
        cex.axis=axis.size)

axis(1, at=30*c(0, 3,  6, 9,  12), 
        labels=c("0", "3", "6",  "9",  "12"), 
        cex.axis=axis.size)

axis(2, at=c(0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        las=2, cex.axis=axis.size)


mtext(side = 2, line = 4, 
	cex=0.8*axis.size, 
	text="  sensitivity")

for(i in 1:length(line_seq_x))
{
	points(x=c(-60,1e4), y=rep(line_seq_x[i],2), type='l', lwd=0.5, col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(10e-10,10e10), type='l', lwd=0.5, col="grey", lty="dashed")
}





polygon(x=c(tt_plot, rev(tt_plot)), 
	  y=c( S2_med_sens_quant[1,], rev(S2_med_sens_quant[3,]) ),
        col=rgb(190/256,190/256,190/256,0.4), border=NA)

points(x=tt_plot, y=S2_med_sens_quant[2,], 
type='l', col="mediumblue", lwd=line.size )


points(x=tt_plot, y=S2_long_sens_quant[2,], 
type='l', lty="longdash", col="mediumblue", lwd=line.size )


points(x=tt_plot, y=S2_short_sens_quant[2,], 
type='l', lty="dotted", col="mediumblue", lwd=line.size )






dev.off()






Spike_med_sens_quant[,180]

RBD_med_sens_quant[,180]

S1_med_sens_quant[,180]

S2_med_sens_quant[,180]






Spike_med_sens_quant[,365]

RBD_med_sens_quant[,365]

S1_med_sens_quant[,365]

S2_med_sens_quant[,365]




Spike_1yr <- matrix(NA, nrow=215, ncol=N_posterior)

for(i in 1:215)
{
	for(j in 1:N_posterior)
	{
		Spike_1yr[i,j] <-	1 - Spike_med_mod[i,j,365]/max(Spike_med_mod[i,j,])
	}
}

quantile( Spike_1yr, prob=c(0.5, 0.025, 0.975), na.rm=TRUE )



RBD_1yr <- matrix(NA, nrow=215, ncol=N_posterior)

for(i in 1:215)
{
	for(j in 1:N_posterior)
	{
		RBD_1yr[i,j] <-	1 - RBD_med_mod[i,j,365]/max(RBD_med_mod[i,j,])
	}
}

quantile( RBD_1yr, prob=c(0.5, 0.025, 0.975), na.rm=TRUE )




S1_1yr <- matrix(NA, nrow=215, ncol=N_posterior)

for(i in 1:215)
{
	for(j in 1:N_posterior)
	{
		S1_1yr[i,j] <-	1 - S1_med_mod[i,j,365]/max(S1_med_mod[i,j,])
	}
}

quantile( S1_1yr, prob=c(0.5, 0.025, 0.975), na.rm=TRUE )




S2_1yr <- matrix(NA, nrow=215, ncol=N_posterior)

for(i in 1:215)
{
	for(j in 1:N_posterior)
	{
		S2_1yr[i,j] <-	1 - S2_med_mod[i,j,365]/max(S2_med_mod[i,j,])
	}
}

quantile( S2_1yr, prob=c(0.5, 0.025, 0.975), na.rm=TRUE )







quantile( colSums(Spike_1yr)/215, prob=c(0.5, 0.025, 0.975) )

quantile( colSums(RBD_1yr)/215, prob=c(0.5, 0.025, 0.975), na.rm=TRUE )

quantile( colSums(S1_1yr)/215, prob=c(0.5, 0.025, 0.975) )

quantile( colSums(S2_1yr)/215, prob=c(0.5, 0.025, 0.975) )



