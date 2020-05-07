
load("C:\\U\\CoronaVirus\\SeroSigPaper\\ANALYSIS\\SeroData\\Data_Process.RData")

part_ID <- as.vector(unique(BICHAT$patient_id))

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


tt_plot <- seq(from=-10, to=30, by=1)
log2 <- log(2)


Spike_file <- paste( "C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/2_IgG_mod/Output/", ant_names, "_IgG_dil_local_med.txt", sep="" )




MCMC_ind <- read.table( Spike_file )

MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]


Spike_quant <- array(NA, dim=c(4, 3, length(tt_plot)) )

for(n in 1:4)
{
	AB_mod <- matrix(NA, nrow=nrow(MCMC_ind), ncol=length(tt_plot))

	###################################
	## Model prediction for participant i
	## Posterior projections
	
	for(i in 1:nrow(MCMC_ind))
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

		AB_mod[i,] = AB_tt
	}

	for(j in 1:length(tt_plot))
	{
		Spike_quant[n,,j] <- quantile( AB_mod[,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
	}
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


tt_plot <- seq(from=-10, to=30, by=1)
log2 <- log(2)


RBD_file <-  paste( "C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/2_IgG_mod/Output/", ant_names, "_IgG_dil_local_med.txt", sep="" )




MCMC_ind <- read.table( RBD_file )

MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]


RBD_quant <- array(NA, dim=c(4, 3, length(tt_plot)) )

for(n in 1:4)
{
	AB_mod <- matrix(NA, nrow=nrow(MCMC_ind), ncol=length(tt_plot))

	###################################
	## Model prediction for participant i
	## Posterior projections
	
	for(i in 1:nrow(MCMC_ind))
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

		AB_mod[i,] = AB_tt
	}

	for(j in 1:length(tt_plot))
	{
		RBD_quant[n,,j] <- quantile( AB_mod[,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
	}
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


tt_plot <- seq(from=-10, to=30, by=1)
log2 <- log(2)


S1_file <-  paste( "C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/2_IgG_mod/Output/", ant_names, "_IgG_dil_local_med.txt", sep="" )




MCMC_ind <- read.table( S1_file )

MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]


S1_quant <- array(NA, dim=c(4, 3, length(tt_plot)) )

for(n in 1:4)
{
	AB_mod <- matrix(NA, nrow=nrow(MCMC_ind), ncol=length(tt_plot))

	###################################
	## Model prediction for participant i
	## Posterior projections
	
	for(i in 1:nrow(MCMC_ind))
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

		AB_mod[i,] = AB_tt
	}

	for(j in 1:length(tt_plot))
	{
		S1_quant[n,,j] <- quantile( AB_mod[,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
	}
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


tt_plot <- seq(from=-10, to=30, by=1)
log2 <- log(2)


S2_file <-  paste( "C:/U/CoronaVirus/SeroSigPaper/ANALYSIS/CoV_kin/2_IgG_mod/Output/", ant_names, "_IgG_dil_local_med.txt", sep="" )




MCMC_ind <- read.table( S2_file )

MCMC_ind <- MCMC_ind[(0.25*nrow(MCMC_ind)):nrow(MCMC_ind),]


S2_quant <- array(NA, dim=c(4, 3, length(tt_plot)) )

for(n in 1:4)
{
	AB_mod <- matrix(NA, nrow=nrow(MCMC_ind), ncol=length(tt_plot))

	###################################
	## Model prediction for participant i
	## Posterior projections
	
	for(i in 1:nrow(MCMC_ind))
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

		AB_mod[i,] = AB_tt
	}

	for(j in 1:length(tt_plot))
	{
		S2_quant[n,,j] <- quantile( AB_mod[,j], prob=c(0.025,0.5,0.975), na.rm=TRUE )
	}
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


tiff(file="Fig3_short_kinetics.tif", width=40, height=25, units="cm", res=500)

lay.mat <- rbind( c( 1, 5, 9,13 ),
                  c( 2, 6,10,14 ),
                  c( 3, 7,11,15 ),
                  c( 4, 8,12,16 ),
                  c(17,17,17,17 ) )
layout(lay.mat, heights=c(10,10,10,10,2), widths=c(1,1,1,1))
layout.show(17)

par(mar=c(4,6.5,3,1))
par(mgp=c(2.5, 1, 0))

point.size = 1.5
lab.size   = 2
axis.size  = 1.5
main.size  = 2


for(k in 1:4)
{

	index <- which( BICHAT$patient_id == part_ID[k] )
	index <- index[order(BICHAT$days_post[index])]


	#############
	## Spike

	line_seq_x <- c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3)
	line_seq_y <- c(0, 5, 10, 15, 20, 25)

	plot(x = BICHAT$days_post[index], y = BICHAT$Spike_IPP_IgG_dil[index],
	     col="orangered", pch=19, cex=point.size,
	     xlim=c(0,25), ylim=c(1e-5,0.02), log="y",
	     yaxt='n', xaxt='n', bty='n',
	     ylab="", xlab="days after symptom onset",
	     main=paste( "patient: ", part_ID[k], sep="" ),
	     cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

	for(i in 1:length(line_seq_x))
	{
		points(x=c(-100,100), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
	}
	for(i in 1:length(line_seq_y))
	{
		points(x=rep(line_seq_y[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
	}

	points(x=c(-100,100), y=rep(1.95e-5,2), type='l', col="black", lty="dashed")
	points(x=c(-100,100), y=rep(0.02,2), type='l', col="black", lty="dashed")


	points(x=tt_plot, y=Spike_quant[k,2,], type='l', col="black")

	polygon(x=c(tt_plot, rev(tt_plot)), 
		  y=c( Spike_quant[k,1,], rev(Spike_quant[k,3,]) ),
	        col=rgb(190/256,190/256,190/256,0.4), border=NA)


	points(x = BICHAT$days_post[index], y = BICHAT$Spike_IPP_IgG_dil[index],
	     col="orangered", pch=19, cex=point.size)


	mtext(side = 2, line = 5, 
	cex=axis.size, 
	text="dilution")



	axis(1, at=c(0, 5, 10, 15, 20, 25), label=c(0, 5, 10, 15, 20, 25), cex.axis=0.8*axis.size )

	axis(2, at = c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03), 
           label = c("0.00001", "", "0.0001", "", "0.001", "", "0.01", ""),
           las=2, cex.axis=axis.size )



	#############
	## RBD

	line_seq_x <- c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3)
	line_seq_y <- c(0, 5, 10, 15, 20, 25)

	plot(x = BICHAT$days_post[index], y = BICHAT$RBD_IPP_IgG_dil[index],
	     col="gold", pch=19, cex=point.size,
	     xlim=c(0,25), ylim=c(1e-5,0.03), log="y",
	     yaxt='n', xaxt='n', bty='n',
	     ylab="", xlab="days after symptom onset",
	     main=paste( "patient: ", part_ID[k], sep="" ),
	     cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

	for(i in 1:length(line_seq_x))
	{
		points(x=c(-100,100), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
	}
	for(i in 1:length(line_seq_y))
	{
		points(x=rep(line_seq_y[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
	}

	points(x=c(-100,100), y=rep(1.95e-5,2), type='l', col="black", lty="dashed")
	points(x=c(-100,100), y=rep(0.02,2), type='l', col="black", lty="dashed")


	points(x=tt_plot, y=RBD_quant[k,2,], type='l', col="black")

	polygon(x=c(tt_plot, rev(tt_plot)), 
		  y=c( RBD_quant[k,1,], rev(RBD_quant[k,3,]) ),
	        col=rgb(190/256,190/256,190/256,0.4), border=NA)


	points(x = BICHAT$days_post[index], y = BICHAT$RBD_IPP_IgG_dil[index],
	     col="gold", pch=19, cex=point.size)


	mtext(side = 2, line = 5, 
	cex=axis.size, 
	text="dilution")



	axis(1, at=c(0, 5, 10, 15, 20, 25), label=c(0, 5, 10, 15, 20, 25), cex.axis=0.8*axis.size )

	axis(2, at = c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03), 
           label = c("0.00001", "", "0.0001", "", "0.001", "", "0.01", ""),
           las=2, cex.axis=axis.size )


	#############
	## S1

	line_seq_x <- c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3)
	line_seq_y <- c(0, 5, 10, 15, 20, 25)

	plot(x = BICHAT$days_post[index], y = BICHAT$S1_NA_IgG_dil[index],
	     col="limegreen", pch=19, cex=point.size,
	     xlim=c(0,25), ylim=c(1e-5,0.03), log="y",
	     yaxt='n', xaxt='n', bty='n',
	     ylab="", xlab="days after symptom onset",
	     main=paste( "patient: ", part_ID[k], sep="" ),
	     cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

	for(i in 1:length(line_seq_x))
	{
		points(x=c(-100,100), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
	}
	for(i in 1:length(line_seq_y))
	{
		points(x=rep(line_seq_y[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
	}

	points(x=c(-100,100), y=rep(1.95e-5,2), type='l', col="black", lty="dashed")
	points(x=c(-100,100), y=rep(0.02,2), type='l', col="black", lty="dashed")



	points(x=tt_plot, y=S1_quant[k,2,], type='l', col="black")

	polygon(x=c(tt_plot, rev(tt_plot)), 
		  y=c( S1_quant[k,1,], rev(S1_quant[k,3,]) ),
	        col=rgb(190/256,190/256,190/256,0.4), border=NA)

	points(x = BICHAT$days_post[index], y = BICHAT$S1_NA_IgG_dil[index],
	     col="limegreen", pch=19, cex=point.size)

	mtext(side = 2, line = 5, 
	cex=axis.size, 
	text="dilution")



	axis(1, at=c(0, 5, 10, 15, 20, 25), label=c(0, 5, 10, 15, 20, 25), cex.axis=0.8*axis.size )

	axis(2, at = c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03), 
           label = c("0.00001", "", "0.0001", "", "0.001", "", "0.01", ""),
           las=2, cex.axis=axis.size )

	#############
	## S2

	line_seq_x <- c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3)
	line_seq_y <- c(0, 5, 10, 15, 20, 25)

	plot(x = BICHAT$days_post[index], y = BICHAT$S2_NA_IgG_dil[index],
	     col="mediumblue", pch=19, cex=point.size,
	     xlim=c(0,25), ylim=c(1e-5,0.03), log="y",
	     yaxt='n', xaxt='n', bty='n',
	     ylab="", xlab="days after symptom onset",
	     main=paste( "patient: ", part_ID[k], sep="" ),
	     cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)

	for(i in 1:length(line_seq_x))
	{
		points(x=c(-100,100), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
	}
	for(i in 1:length(line_seq_y))
	{
		points(x=rep(line_seq_y[i],2), y=c(1e-10,1e10), type='l', col="grey", lty="dashed")
	}

	points(x=c(-100,100), y=rep(1.95e-5,2), type='l', col="black", lty="dashed")
	points(x=c(-100,100), y=rep(0.02,2), type='l', col="black", lty="dashed")


	points(x=tt_plot, y=S2_quant[k,2,], type='l', col="black")

	polygon(x=c(tt_plot, rev(tt_plot)), 
		  y=c( S2_quant[k,1,], rev(S2_quant[k,3,]) ),
	        col=rgb(190/256,190/256,190/256,0.4), border=NA)

	points(x = BICHAT$days_post[index], y = BICHAT$S2_NA_IgG_dil[index],
	     col="mediumblue", pch=19, cex=point.size)

	mtext(side = 2, line = 5, 
	cex=axis.size, 
	text="dilution")



	axis(1, at=c(0, 5, 10, 15, 20, 25), label=c(0, 5, 10, 15, 20, 25), cex.axis=0.8*axis.size )

	axis(2, at = c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01, 0.03), 
           label = c("0.00001", "", "0.0001", "", "0.001", "", "0.01", ""),
           las=2, cex.axis=axis.size )
}


###############	
##           ##
##  LEGEND   ##
##           ##
###############

par(mar = c(0,0,0,0))
plot.new()

legend(x='center', 
       legend = c( expression(paste( "S"^"tri", " IgG" )), expression(paste( "RBD IgG" )), expression(paste( "S1 IgG" )), expression(paste( "S2 IgG" )) ), 
       fill = c("orangered", "gold",  "limegreen", "mediumblue"), 
       border = c("orangered", "gold",  "limegreen", "mediumblue"), 
       ncol=4, cex=2, bty="n" )


dev.off()








