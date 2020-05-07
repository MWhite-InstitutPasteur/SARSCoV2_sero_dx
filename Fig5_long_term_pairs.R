
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



Spike_file <-  paste( "C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/2_IgG_mod/Output/", ant_names, "_IgG_dil_local_med.txt", sep="" )





MCMC_ind <- read.table( Spike_file )

MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]

N_posterior <- nrow(MCMC_ind)

Spike_mod <- array(NA, dim=c(N_part, N_posterior, N_tt_plot) )
Spike_quant <- array(NA, dim=c(N_part, 3, N_tt_plot) )

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

		Spike_mod[n,i,] = AB_tt
	}

	for(j in 1:length(tt_plot))
	{
		Spike_quant[n,,j] <- quantile( Spike_mod[n,,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
	}
}



Spike_sens <- matrix(NA, nrow=N_posterior, ncol=N_tt_plot)

for(i in 1:N_posterior)
{
	Spike_sens[i,] <- colSums(Spike_mod[,i,] > Spike_cut)/N_part
}

Spike_sens_quant <- matrix(NA, nrow=3, ncol=N_tt_plot)
for(j in 1:N_tt_plot)
{
	Spike_sens_quant[,j] <- quantile( Spike_sens[,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
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



RBD_file <-  paste( "C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/2_IgG_mod/Output/", ant_names, "_IgG_dil_local_med.txt", sep="" )






MCMC_ind <- read.table( RBD_file )

MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]

N_posterior <- nrow(MCMC_ind)

RBD_mod <- array(NA, dim=c(N_part, N_posterior, N_tt_plot) )
RBD_quant <- array(NA, dim=c(N_part, 3, N_tt_plot) )

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

		RBD_mod[n,i,] = AB_tt
	}

	for(j in 1:length(tt_plot))
	{
		RBD_quant[n,,j] <- quantile( RBD_mod[n,,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
	}
}



RBD_sens <- matrix(NA, nrow=N_posterior, ncol=N_tt_plot)

for(i in 1:N_posterior)
{
	RBD_sens[i,] <- colSums(RBD_mod[,i,] > RBD_cut)/N_part
}

RBD_sens_quant <- matrix(NA, nrow=3, ncol=N_tt_plot)
for(j in 1:N_tt_plot)
{
	RBD_sens_quant[,j] <- quantile( RBD_sens[,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
}




#######################################################################
## S1 data

ant_names <- "S1_NA"

S1_data = read.table( paste("C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/1_Data/", ant_names, "_IgG_dil.txt", sep="") )


N_tt <- 15

N_part <- nrow(S1_data)


AB_max = 0.2
tt_max = max( RBD_data[,3:(2+N_tt)] )


AB_min = 1.95e-5




S1_file <-  paste( "C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/2_IgG_mod/Output/", ant_names, "_IgG_dil_local_med.txt", sep="" )





MCMC_ind <- read.table( Spike_file )

MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]

N_posterior <- nrow(MCMC_ind)

S1_mod <- array(NA, dim=c(N_part, N_posterior, N_tt_plot) )
S1_quant <- array(NA, dim=c(N_part, 3, N_tt_plot) )

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

		S1_mod[n,i,] = AB_tt
	}

	for(j in 1:length(tt_plot))
	{
		S1_quant[n,,j] <- quantile( S1_mod[n,,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
	}
}



S1_sens <- matrix(NA, nrow=N_posterior, ncol=N_tt_plot)

for(i in 1:N_posterior)
{
	S1_sens[i,] <- colSums(S1_mod[,i,] > S1_cut)/N_part
}

S1_sens_quant <- matrix(NA, nrow=3, ncol=N_tt_plot)
for(j in 1:N_tt_plot)
{
	S1_sens_quant[,j] <- quantile( S1_sens[,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
}



#######################################################################
## S2 data

ant_names <- "S2_NA"

S2_data = read.table( paste("C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/1_Data/", ant_names, "_IgG_dil.txt", sep="") )


N_tt <- 15

N_part <- nrow(S1_data)


AB_max = 0.2
tt_max = max( RBD_data[,3:(2+N_tt)] )


AB_min = 1.95e-5






S2_file <-  paste( "C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/2_IgG_mod/Output/", ant_names, "_IgG_dil_local_med.txt", sep="" )



MCMC_ind <- read.table( S2_file )

MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]

N_posterior <- nrow(MCMC_ind)

S2_mod <- array(NA, dim=c(N_part, N_posterior, N_tt_plot) )
S2_quant <- array(NA, dim=c(N_part, 3, N_tt_plot) )

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

		S2_mod[n,i,] = AB_tt
	}

	for(j in 1:length(tt_plot))
	{
		S2_quant[n,,j] <- quantile( S2_mod[n,,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
	}
}



S2_sens <- matrix(NA, nrow=N_posterior, ncol=N_tt_plot)

for(i in 1:N_posterior)
{
	S2_sens[i,] <- colSums(S2_mod[,i,] > S2_cut)/N_part
}

S2_sens_quant <- matrix(NA, nrow=3, ncol=N_tt_plot)
for(j in 1:N_tt_plot)
{
	S2_sens_quant[,j] <- quantile( S2_sens[,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
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

line_seq_x <- c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01)
line_seq_y <- c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01)



tiff(file="Fig5_long_term_pairs.tif", width=30, height=22, units="cm", res=500)

lay.mat <- rbind( c( 1, 2 ),
                  c( 3, 4 ),
                  c( 5, 5 ) )
layout(lay.mat, heights=c(10,10,1.5), widths=c(1,1))
layout.show(5)

par(mar=c(5,8,3,1))
par(mgp=c(3, 1, 0))

point.size = 0.75
lab.size   = 2
axis.size  = 1.5
main.size  = 1.75


###################################
## Panel 1: Spike_IPP 

plot( x = BICHAT$RBD_IPP_IgG_dil, y = BICHAT$Spike_IPP_IgG_dil,
      xlim=c(1e-5, 0.03), ylim=c(1e-5, 0.03), log="xy",
      pch=19, col="yellowgreen",
	yaxt='n', xaxt='n', bty='n', 
	xlab=expression(paste( "anti-RBD IgG dilution", sep="" )),
      ylab="", 
	main="(A) Measured antibody responses in first month",
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
}

points(x=rep(1.95e-5,2), y=c(1e-10,1e10), type='l', col="black", lty="dashed")
points(x=rep(0.02,2), y=c(1e-10,1e10), type='l', col="black", lty="dashed")

points(x=c(1e-10,1e10), y=rep(1.95e-5,2),  type='l', col="black", lty="dashed")
points(x=c(1e-10,1e10), y=rep(0.02,2), type='l', col="black", lty="dashed")

mtext(side = 2, line = 5, 
cex=axis.size, 
text=expression(paste( "anti-S"^"tri", " IgG dilution", sep="" )))

points( x = STRAS$RBD_IPP_IgG_dil, y = STRAS$Spike_IPP_IgG_dil,
        pch=19, col="sienna1" )

points( x = COCHIN$RBD_IPP_IgG_dil, y = COCHIN$Spike_IPP_IgG_dil,
        pch=19, col="darkmagenta" )


points( x = EFS$RBD_IPP_IgG_dil, y = EFS$Spike_IPP_IgG_dil,
        pch=19, col="dodgerblue" )

points( x = TRC$RBD_IPP_IgG_dil, y = TRC$Spike_IPP_IgG_dil,
        pch=19, col="royalblue" )

points( x = PNC$RBD_IPP_IgG_dil, y = PNC$Spike_IPP_IgG_dil,
        pch=19, col="cornflowerblue" )





axis(1, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
        cex.axis=axis.size )

axis(2, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
        las=2, cex.axis=axis.size )


###################################
## Panel 2: 3 months 


plot( x = 100, y = 100,
      xlim=c(1e-5, 0.03), ylim=c(1e-5, 0.03), log="xy",
      pch=19, col="yellowgreen",
	yaxt='n', xaxt='n', bty='n', 
	xlab=expression(paste( "anti-RBD IgG dilution", sep="" )), ylab="",
	main="(B) Modelled antibody responses at 3 months",
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
}

points(x=rep(1.95e-5,2), y=c(1e-10,1e10), type='l', col="black", lty="dashed")
points(x=rep(0.02,2), y=c(1e-10,1e10), type='l', col="black", lty="dashed")

points(x=c(1e-10,1e10), y=rep(1.95e-5,2),  type='l', col="black", lty="dashed")
points(x=c(1e-10,1e10), y=rep(0.02,2), type='l', col="black", lty="dashed")

mtext(side = 2, line = 5, 
cex=axis.size, 
text=expression(paste( "anti-S"^"tri", " IgG dilution", sep="" )))

points( x = RBD_quant[,2,90], y = Spike_quant[,2,90],
        pch=19, col=c( rep("yellowgreen",4), rep("sienna1",162), rep("darkmagenta", 49)) ) 
        
points( x = EFS$RBD_IPP_IgG_dil, y = EFS$Spike_IPP_IgG_dil,
        pch=19, col="dodgerblue" )

points( x = TRC$RBD_IPP_IgG_dil, y = TRC$Spike_IPP_IgG_dil,
        pch=19, col="royalblue" )

points( x = PNC$RBD_IPP_IgG_dil, y = PNC$Spike_IPP_IgG_dil,
        pch=19, col="cornflowerblue" )



points( x = RBD_quant[30,c(1,3),90], 
        y = Spike_quant[30,c(2,2),90],
        type='l', col="black" ) 

points( x = RBD_quant[30,c(2,2),90], 
        y = Spike_quant[30,c(1,3),90],
        type='l', col="black" ) 

axis(1, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
        cex.axis=axis.size )

axis(2, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
        las=2, cex.axis=axis.size )




###################################
## Panel 3: 6 months 


plot( x = 100, y = 100,
      xlim=c(1e-5, 0.03), ylim=c(1e-5, 0.03), log="xy",
      pch=19, col="yellowgreen",
	yaxt='n', xaxt='n', bty='n', 
	xlab=expression(paste( "anti-RBD IgG dilution", sep="" )), ylab="",
	main="(C) Modelled antibody responses at 6 months",
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
}

points(x=rep(1.95e-5,2), y=c(1e-10,1e10), type='l', col="black", lty="dashed")
points(x=rep(0.02,2), y=c(1e-10,1e10), type='l', col="black", lty="dashed")

points(x=c(1e-10,1e10), y=rep(1.95e-5,2),  type='l', col="black", lty="dashed")
points(x=c(1e-10,1e10), y=rep(0.02,2), type='l', col="black", lty="dashed")

mtext(side = 2, line = 5, 
cex=axis.size, 
text=expression(paste( "anti-S"^"tri", " IgG dilution", sep="" )))

points( x = RBD_quant[,2,180], y = Spike_quant[,2,180],
        pch=19, col=c( rep("yellowgreen",4), rep("sienna1",162), rep("darkmagenta", 49)) )  
        
points( x = EFS$RBD_IPP_IgG_dil, y = EFS$Spike_IPP_IgG_dil,
        pch=19, col="dodgerblue" )

points( x = TRC$RBD_IPP_IgG_dil, y = TRC$Spike_IPP_IgG_dil,
        pch=19, col="royalblue" )

points( x = PNC$RBD_IPP_IgG_dil, y = PNC$Spike_IPP_IgG_dil,
        pch=19, col="cornflowerblue" )

points( x = RBD_quant[30,c(1,3),180], 
        y = Spike_quant[30,c(2,2),180],
        type='l', col="black" ) 

points( x = RBD_quant[30,c(2,2),180], 
        y = Spike_quant[30,c(1,3),180],
        type='l', col="black" ) 

axis(1, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
        cex.axis=axis.size )

axis(2, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
        las=2, cex.axis=axis.size )



###################################
## Panel 4: RBD_IPP 


plot( x = 100, y = 100,
      xlim=c(1e-5, 0.03), ylim=c(1e-5, 0.03), log="xy",
      pch=19, col="yellowgreen",
	yaxt='n', xaxt='n', bty='n', 
	xlab=expression(paste( "anti-RBD IgG dilution", sep="" )), ylab="",
	main="(D) Modelled antibody responses at 12 months",
	cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

for(i in 1:length(line_seq_x))
{
	points(x=c(1e-10,1e10), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_y))
{
	points(x=rep(line_seq_y[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
}

points(x=rep(1.95e-5,2), y=c(1e-10,1e10), type='l', col="black", lty="dashed")
points(x=rep(0.02,2), y=c(1e-10,1e10), type='l', col="black", lty="dashed")

points(x=c(1e-10,1e10), y=rep(1.95e-5,2),  type='l', col="black", lty="dashed")
points(x=c(1e-10,1e10), y=rep(0.02,2), type='l', col="black", lty="dashed")

mtext(side = 2, line = 5, 
cex=axis.size, 
text=expression(paste( "anti-S"^"tri", " IgG dilution", sep="" )))

points( x = RBD_quant[,2,360], y = Spike_quant[,2,360],
        pch=19, col=c( rep("yellowgreen",4), rep("sienna1",162), rep("darkmagenta", 49)) ) 
        
points( x = EFS$RBD_IPP_IgG_dil, y = EFS$Spike_IPP_IgG_dil,
        pch=19, col="dodgerblue" )

points( x = TRC$RBD_IPP_IgG_dil, y = TRC$Spike_IPP_IgG_dil,
        pch=19, col="royalblue" )

points( x = PNC$RBD_IPP_IgG_dil, y = PNC$Spike_IPP_IgG_dil,
        pch=19, col="cornflowerblue" )

points( x = RBD_quant[30,c(1,3),360], 
        y = Spike_quant[30,c(2,2),360],
        type='l', col="black" ) 

points( x = RBD_quant[30,c(2,2),360], 
        y = Spike_quant[30,c(1,3),360],
        type='l', col="black" ) 

axis(1, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
        cex.axis=axis.size )

axis(2, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
        las=2, cex.axis=axis.size )











###############	
##           ##
##  LEGEND   ##
##           ##
###############

par(mar = c(0,0,0,0))
plot.new()

legend(x='center', 
       legend = c("Bichat",     "France negative",
                  "Strasbourg", "Thai negative", 
                  "Cochin",     "Peru negative"), 
       fill = c("yellowgreen", "dodgerblue",
                "sienna1",     "royalblue",
                "darkmagenta", "cornflowerblue"), 
       border = c("yellowgreen", "dodgerblue",
                "sienna1",     "royalblue",
                "darkmagenta", "cornflowerblue"), 
       ncol=3, cex=2, bty="n" )



dev.off()








