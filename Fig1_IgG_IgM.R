library(binom)

load("C:\\U\\CoronaVirus\\SeroSigPaper\\ANALYSIS\\SeroData\\Data_Process.RData")




N_negative <- N_trc + N_pnc + N_efs

N_positive <- N_bichat + N_stras + N_cochin



t.test( x = log(c(BICHAT$Spike_IPP_IgG_dil, STRAS$Spike_IPP_IgG_dil, COCHIN$Spike_IPP_IgG_dil)),
 	  y = log(c(TRC$Spike_IPP_IgG_dil, PNC$Spike_IPP_IgG_dil, EFS$Spike_IPP_IgG_dil)) )


t.test( x = log(c(BICHAT$RBD_IPP_IgG_dil, STRAS$RBD_IPP_IgG_dil, COCHIN$RBD_IPP_IgG_dil)),
 	  y = log(c(TRC$RBD_IPP_IgG_dil, PNC$RBD_IPP_IgG_dil, EFS$RBD_IPP_IgG_dil)) )



t.test( x = log(c(BICHAT$S1_NA_dil, STRAS$S1_NA_IgG_dil, COCHIN$S1_NA_IgG_dil)),
 	  y = log(c(TRC$S1_NA_IgG_dil, PNC$S1_NA_IgG_dil, EFS$S1_NA_IgG_dil)) )



t.test( x = log(c(BICHAT$S2_NA_IgG_dil, STRAS$S2_NA_IgG_dil, COCHIN$S2_NA_IgG_dil)),
 	  y = log(c(TRC$S2_NA_IgG_dil, PNC$S2_NA_IgG_dil, EFS$S2_NA_IgG_dil)) )


t.test( x = log(c(BICHAT$Spike_IPP_IgM_dil, STRAS$Spike_IPP_IgM_dil, COCHIN$Spike_IPP_IgM_dil)),
 	  y = log(c(TRC$Spike_IPP_IgM_dil, PNC$Spike_IPP_IgM_dil)) )


t.test( x = log(c(BICHAT$RBD_IPP_IgM_dil, STRAS$RBD_IPP_IgM_dil, COCHIN$RBD_IPP_IgM_dil)),
 	  y = log(c(TRC$RBD_IPP_IgM_dil, PNC$RBD_IPP_IgM_dil)) )



t.test( x = log(c(BICHAT$S1_NA_IgM_dil, STRAS$S1_NA_IgM_dil, COCHIN$S1_NA_IgM_dil)),
 	  y = log(c(TRC$S1_NA_IgM_dil, PNC$S1_NA_IgM_dil)) )



t.test( x = log(c(BICHAT$S2_NA_IgM_dil, STRAS$S2_NA_IgM_dil, COCHIN$S2_NA_IgM_dil)),
 	  y = log(c(TRC$S2_NA_IgM_dil, PNC$S2_NA_IgM_dil)) )




####################################
####################################
##                                ##
##  ONE AT A TIME                 ##
##                                ##
####################################
####################################


N_ant <- 8
N_SS <- 50000

SS_cut <- exp(seq(from=log(1e-5), to=log(0.3), length=N_SS))

sens_mat <- matrix(NA, ncol=N_ant, nrow=N_SS)
spec_mat <- matrix(NA, ncol=N_ant, nrow=N_SS)

colnames(sens_mat) <- c( "Spike_IPP_IgG", "RBD_IPP_IgG", "S1_NA_IgG", "S2_NA_IgG", "Spike_IPP_IgM", "RBD_IPP_IgM", "S1_NA_IgM", "S2_NA_IgM")
colnames(spec_mat) <- c( "Spike_IPP_IgG", "RBD_IPP_IgG", "S1_NA_IgG", "S2_NA_IgG", "Spike_IPP_IgM", "RBD_IPP_IgM", "S1_NA_IgM", "S2_NA_IgM")

for(i in 1:N_SS)
{
	sens_mat[i,1] <- length(which( c( BICHAT$Spike_IPP_IgG_dil, STRAS$Spike_IPP_IgG_dil, COCHIN$Spike_IPP_IgG_dil) > SS_cut[i] ))/( N_positive )
	spec_mat[i,1] <- length(which( c( TRC$Spike_IPP_IgG_dil, PNC$Spike_IPP_IgG_dil, EFS$Spike_IPP_IgG_dil) < SS_cut[i] ))/( N_negative )

	sens_mat[i,2] <- length(which( c( BICHAT$RBD_IPP_IgG_dil, STRAS$RBD_IPP_IgG_dil, COCHIN$RBD_IPP_IgG_dil) > SS_cut[i] ))/( N_positive )
	spec_mat[i,2] <- length(which( c( TRC$RBD_IPP_IgG_dil, PNC$RBD_IPP_IgG_dil, EFS$RBD_IPP_IgG_dil) < SS_cut[i] ))/( N_negative )

	sens_mat[i,3] <- length(which( c( BICHAT$S1_NA_IgG_dil, STRAS$S1_NA_IgG_dil, COCHIN$S1_NA_IgG_dil) > SS_cut[i] ))/( N_positive )
	spec_mat[i,3] <- length(which( c( TRC$S1_NA_IgG_dil, PNC$S1_NA_IgG_dil, EFS$S1_NA_IgG_dil) < SS_cut[i] ))/( N_negative )

	sens_mat[i,4] <- length(which( c( BICHAT$S2_NA_IgG_dil, STRAS$S2_NA_IgG_dil, COCHIN$S2_NA_IgG_dil) > SS_cut[i] ))/( N_positive )
	spec_mat[i,4] <- length(which( c( TRC$S2_NA_IgG_dil, PNC$S2_NA_IgG_dil, EFS$S2_NA_IgG_dil) < SS_cut[i] ))/( N_negative )
}

for(i in 1:N_SS)
{
	sens_mat[i,5] <- length(which( c( BICHAT$Spike_IPP_IgM_dil, STRAS$Spike_IPP_IgM_dil, COCHIN$Spike_IPP_IgM_dil) > SS_cut[i] ))/( N_positive )
	spec_mat[i,5] <- length(which( c( TRC$Spike_IPP_IgM_dil, PNC$Spike_IPP_IgM_dil ) < SS_cut[i] ))/( N_trc + N_pnc )

	sens_mat[i,6] <- length(which( c( BICHAT$RBD_IPP_IgM_dil, STRAS$RBD_IPP_IgM_dil, COCHIN$RBD_IPP_IgM_dil) > SS_cut[i] ))/( N_positive )
	spec_mat[i,6] <- length(which( c( TRC$RBD_IPP_IgM_dil, PNC$RBD_IPP_IgM_dil ) < SS_cut[i] ))/( N_trc + N_pnc )

	sens_mat[i,7] <- length(which( c( BICHAT$S1_NA_IgM_dil, STRAS$S1_NA_IgM_dil, COCHIN$S1_NA_IgM_dil) > SS_cut[i] ))/( N_positive )
	spec_mat[i,7] <- length(which( c( TRC$S1_NA_IgM_dil, PNC$S1_NA_IgM_dil ) < SS_cut[i] ))/( N_trc + N_pnc )

	sens_mat[i,8] <- length(which( c( BICHAT$S2_NA_IgM_dil, STRAS$S2_NA_IgM_dil, COCHIN$S2_NA_IgM_dil) > SS_cut[i] ))/( N_positive )
	spec_mat[i,8] <- length(which( c( TRC$S2_NA_IgM_dil, PNC$S2_NA_IgM_dil ) < SS_cut[i] ))/( N_trc + N_pnc )
}




######################################
## Calculate Area Under Curve (AUC)


AUC_one_ant <- rep(NA, N_ant)

for(j in 1:N_ant)
{
	AUC_one_ant[j] <- sum( (sens_mat[1:(N_SS-1),j] - sens_mat[2:N_SS,j])*
                             0.5*(spec_mat[1:(N_SS-1),j] + spec_mat[2:N_SS,j]) )
}





sens_target <- matrix(NA, nrow=3, ncol=N_ant)
colnames(sens_target) <- colnames(sens_mat)
rownames(sens_target) <- c("sens_99", "sens_spec", "spec_99")

spec_target <- matrix(NA, nrow=3, ncol=N_ant)
colnames(spec_target) <- colnames(spec_mat)
rownames(spec_target) <- c("sens_99", "sens_spec", "spec_99")

for(j in 1:N_ant)
{
	sens_target[1,j] <- sens_mat[max(which(sens_mat[,j] > 0.99)),j]
	spec_target[1,j] <- spec_mat[max(which(sens_mat[,j] > 0.99)),j]

	sens_target[2,j] <- sens_mat[which.min(abs(sens_mat[,j] - spec_mat[,j])),j]
	spec_target[2,j] <- spec_mat[which.min(abs(sens_mat[,j] - spec_mat[,j])),j]

	sens_target[3,j] <- sens_mat[min(which(spec_mat[,j] > 0.99)),j]
	spec_target[3,j] <- spec_mat[min(which(spec_mat[,j] > 0.99)),j]
}




sens_target_lwr <- sens_target

for(i in 1:nrow(sens_target_lwr))
{
	for(j in 1:ncol(sens_target_lwr))
	{
		sens_target_lwr[i,j] <- binom.confint( sens_target_lwr[i,j]*N_positive, N_positive, method="wilson")[1,5]
		sens_target_lwr[i,j] <- round( 100*sens_target_lwr[i,j], 1)  
	}

}

sens_target_upr <- sens_target

for(i in 1:nrow(sens_target_upr))
{
	for(j in 1:ncol(sens_target_upr))
	{
		sens_target_upr[i,j] <- binom.confint( sens_target_upr[i,j]*N_positive, N_positive, method="wilson")[1,6]
		sens_target_upr[i,j] <- round( 100*sens_target_upr[i,j], 1)  
	}
}


spec_target_lwr <- spec_target

for(i in 1:nrow(spec_target_lwr))
{
	for(j in 1:4)
	{
		spec_target_lwr[i,j] <- binom.confint( spec_target_lwr[i,j]*N_negative, N_negative, method="wilson")[1,5]
		spec_target_lwr[i,j] <- round( 100*spec_target_lwr[i,j], 1)  
	}

	for(j in 5:8)
	{
		spec_target_lwr[i,j] <- binom.confint( spec_target_lwr[i,j]*(N_trc + N_pnc), (N_trc + N_pnc), method="wilson")[1,5]
		spec_target_lwr[i,j] <- round( 100*spec_target_lwr[i,j], 1)  
	}
}

spec_target_upr <- spec_target

for(i in 1:nrow(spec_target_upr))
{
	for(j in 1:4)
	{
		spec_target_upr[i,j] <- binom.confint( spec_target_upr[i,j]*N_negative, N_negative, method="wilson")[1,6]
		spec_target_upr[i,j] <- round( 100*spec_target_upr[i,j], 1)  
	}

	for(j in 5:8)
	{
		spec_target_upr[i,j] <- binom.confint( spec_target_upr[i,j]*(N_trc + N_pnc), (N_trc + N_pnc), method="wilson")[1,6]
		spec_target_upr[i,j] <- round( 100*spec_target_upr[i,j], 1)  
	}
}


sens_target <- round( 100*sens_target, 1)



spec_target <- round( 100*spec_target, 1)
















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




tiff(file="Fig1_IgG_IgM.tif", width=40, height=20, units="cm", res=500)

lay.mat <- rbind( c( 1, 2, 3, 4, 9 ),
                  c( 5, 6, 7, 8, 9 ) )
layout(lay.mat, heights=c(1,1), widths=c(1,1,1,1,2))
layout.show(9)

par(mar=c(3,5,3,1))
par(mgp=c(2.5, 1, 0))

point.size = 0.75
lab.size   = 2
axis.size  = 1.1
main.size  = 1.5


line_seq_x <- c(0.00001, 0.00003, 0.0001, 0.0003, 0.001, 0.003, 0.01)


#####################################
##                                  
##  PANEL 1                     
##  Spike

boxplot( BICHAT$Spike_IPP_IgG_dil,
	   STRAS$Spike_IPP_IgG_dil,
	   COCHIN$Spike_IPP_IgG_dil,
 	   TRC$Spike_IPP_IgG_dil,
	   PNC$Spike_IPP_IgG_dil, 
	   EFS$Spike_IPP_IgG_dil, 
	   log="y",
	   pch=19, yaxt='n', xaxt='n',
	   ylim=c(1e-5, 0.03),
	   col=c("orangered", "orangered", "orangered", "orangered", "orangered", "orangered"),
	   ylab="", 
	   main=expression(paste( "(A) anti-S"^"tri", " IgG antibodies", sep="" )),
	   cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


for(i in 1:length(line_seq_x))
{
	points(x=c(0,100), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
}

mtext(side = 2, line = 3.5, 
cex=axis.size, 
text="     antibody dilution")



axis(1, at=c(1, 3, 5), label=c("Bichat", "Cochin",  "PNC"), cex.axis=0.8*axis.size )

axis(1, at=c(2, 4, 6), label=c( "Sburg", "TRC", "EFS"), cex.axis=0.8*axis.size )


axis(2, at=c(0.00001, 0.0001, 0.001, 0.01), 
        label=c("0.00001", "0.0001", "0.001", "0.01"),
        las=2, cex.axis=axis.size )



#####################################
##                                  
##  PANEL 2                     
##  RBD

par(mar=c(3,4,3,1))


boxplot( BICHAT$RBD_IPP_IgG_dil,
	   STRAS$RBD_IPP_IgG_dil,
	   COCHIN$RBD_IPP_IgG_dil,
 	   TRC$RBD_IPP_IgG_dil,
	   PNC$RBD_IPP_IgG_dil, 
	   EFS$RBD_IPP_IgG_dil, 
	   log="y",
	   pch=19, yaxt='n', xaxt='n',
	   ylim=c(1e-5, 0.03),
	   col=c("gold", "gold", "gold", "gold", "gold", "gold"),
	   ylab="", 
	   main=expression(paste( "(B) anti-RBD IgG antibodies", sep="" )),
	   cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


for(i in 1:length(line_seq_x))
{
	points(x=c(0,100), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
}


axis(1, at=c(1, 3, 5), label=c("Bichat", "Cochin",  "PNC"), cex.axis=0.8*axis.size )

axis(1, at=c(2, 4, 6), label=c( "Sburg", "TRC", "EFS"), cex.axis=0.8*axis.size )


axis(2, at=c(0.00001, 0.0001, 0.001, 0.01), 
        label=c("0.00001", "0.0001", "0.001", "0.01"),
        las=2, cex.axis=axis.size )



#####################################
##                                  
##  PANEL 3                     
##  S1_NA

boxplot( BICHAT$S1_NA_IgG_dil,
	   STRAS$S1_NA_IgG_dil,
	   COCHIN$S1_NA_IgG_dil,
 	   TRC$S1_NA_IgG_dil,
	   PNC$S1_NA_IgG_dil, 
	   EFS$S1_NA_IgG_dil, 
         log="y",
	   pch=19, yaxt='n', xaxt='n',
	   ylim=c(1e-5, 0.03),
	   col=c("limegreen", "limegreen", "limegreen", "limegreen", "limegreen", "limegreen"),
	   ylab="", 
	   main=expression(paste( "(C) anti-S1 IgG antibodies", sep="" )),
	   cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


for(i in 1:length(line_seq_x))
{
	points(x=c(0,100), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
}




axis(1, at=c(1, 3, 5), label=c("Bichat", "Cochin",  "PNC"), cex.axis=0.8*axis.size )

axis(1, at=c(2, 4, 6), label=c( "Sburg", "TRC", "EFS"), cex.axis=0.8*axis.size )


axis(2, at=c(0.00001, 0.0001, 0.001, 0.01), 
        label=c("0.00001", "0.0001", "0.001", "0.01"),
        las=2, cex.axis=axis.size )


#####################################
##                                  
##  PANEL 4                     
##  S2_NA

boxplot( BICHAT$S2_NA_IgG_dil,
	   STRAS$S1_NA_IgG_dil,
	   COCHIN$S1_NA_IgG_dil,
 	   TRC$S1_NA_IgG_dil,
	   PNC$S1_NA_IgG_dil, 
	   EFS$S1_NA_IgG_dil, 
         log="y",
	   pch=19, yaxt='n', xaxt='n',
	   ylim=c(1e-5, 0.03),
	   col=c("mediumblue", "mediumblue", "mediumblue", "mediumblue", "mediumblue", "mediumblue"),
	   ylab="", 
	   main=expression(paste( "(D) anti-S2 IgG antibodies", sep="" )),
	   cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)





for(i in 1:length(line_seq_x))
{
	points(x=c(0,100), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
}


axis(1, at=c(1, 3, 5), label=c("Bichat", "Cochin",  "PNC"), cex.axis=0.8*axis.size )

axis(1, at=c(2, 4, 6), label=c( "Sburg", "TRC", "EFS"), cex.axis=0.8*axis.size )


axis(2, at=c(0.00001, 0.0001, 0.001, 0.01), 
        label=c("0.00001", "0.0001", "0.001", "0.01"),
        las=2, cex.axis=axis.size )


#####################################
##                                  
##  PANEL 5                     
##  Spike IgM

par(mar=c(3,5,3,1))


boxplot( BICHAT$Spike_IPP_IgM_dil,
	   STRAS$Spike_IPP_IgM_dil,
	   COCHIN$Spike_IPP_IgM_dil,
 	   TRC$Spike_IPP_IgM_dil,
	   PNC$Spike_IPP_IgM_dil, 
	   EFS$Spike_IPP_IgM_dil, 
         log="y",
	   pch=19, yaxt='n', xaxt='n',
	   ylim=c(1e-5, 0.03),
	   col=c("orangered", "orangered", "orangered", "orangered", "orangered", "orangered"),
	   ylab="", 
	   main=expression(paste( "(E) anti-S"^"tri", " IgM antibodies", sep="" )),
	   cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


for(i in 1:length(line_seq_x))
{
	points(x=c(0,100), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
}

mtext(side = 2, line = 3.5, 
cex=axis.size, 
text="     antibody dilution")



axis(1, at=c(1, 3, 5), label=c("Bichat", "Cochin",  "PNC"), cex.axis=0.8*axis.size )

axis(1, at=c(2, 4, 6), label=c( "Sburg", "TRC", "EFS"), cex.axis=0.8*axis.size )


axis(2, at=c(0.00001, 0.0001, 0.001, 0.01), 
        label=c("0.00001", "0.0001", "0.001", "0.01"),
        las=2, cex.axis=axis.size )



#####################################
##                                  
##  PANEL 6                     
##  RBD IgM 

par(mar=c(3,4,3,1))


boxplot( BICHAT$RBD_IPP_IgM_dil,
	   STRAS$RBD_IPP_IgM_dil,
	   COCHIN$RBD_IPP_IgM_dil,
 	   TRC$RBD_IPP_IgM_dil,
	   PNC$RBD_IPP_IgM_dil, 
	   EFS$RBD_IPP_IgM_dil, 
         log="y",
	   pch=19, yaxt='n', xaxt='n',
	   ylim=c(1e-5, 0.03),
	   col=c("gold", "gold", "gold", "gold", "gold", "gold"),
	   ylab="", 
	   main=expression(paste( "(F) anti-RBD IgM antibodies", sep="" )),
	   cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


for(i in 1:length(line_seq_x))
{
	points(x=c(0,100), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
}


axis(1, at=c(1, 3, 5), label=c("Bichat", "Cochin",  "PNC"), cex.axis=0.8*axis.size )

axis(1, at=c(2, 4, 6), label=c( "Sburg", "TRC", "EFS"), cex.axis=0.8*axis.size )

axis(2, at=c(0.00001, 0.0001, 0.001, 0.01), 
        label=c("0.00001", "0.0001", "0.001", "0.01"),
        las=2, cex.axis=axis.size )



#####################################
##                                  
##  PANEL 7                     
##  S1_NA IgM

boxplot( BICHAT$S1_NA_IgM_dil,
	   STRAS$S1_NA_IgM_dil,
	   COCHIN$S1_NA_IgM_dil,
 	   TRC$S1_NA_IgM_dil,
	   PNC$S1_NA_IgM_dil, 
	   EFS$S1_NA_IgM_dil, 
         log="y",
	   pch=19, yaxt='n', xaxt='n',
	   ylim=c(1e-5, 0.03),
	   col=c("limegreen", "limegreen", "limegreen", "limegreen", "limegreen", "limegreen"),
	   ylab="", 
	   main=expression(paste( "(G) anti-S1 IgM antibodies", sep="" )),
	   cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)


for(i in 1:length(line_seq_x))
{
	points(x=c(0,100), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
}



axis(1, at=c(1, 3, 5), label=c("Bichat", "Cochin",  "PNC"), cex.axis=0.8*axis.size )

axis(1, at=c(2, 4, 6), label=c( "Sburg", "TRC", "EFS"), cex.axis=0.8*axis.size )


axis(2, at=c(0.00001, 0.0001, 0.001, 0.01), 
        label=c("0.00001", "0.0001", "0.001", "0.01"),
        las=2, cex.axis=axis.size )


#####################################
##                                  
##  PANEL 8                     
##  S2_NA IgM

boxplot( BICHAT$S2_NA_IgM_dil,
	   STRAS$S1_NA_IgM_dil,
	   COCHIN$S1_NA_IgM_dil,
 	   TRC$S1_NA_IgM_dil,
	   PNC$S1_NA_IgM_dil, 
	   EFS$S1_NA_IgM_dil, 
	   log="y",
	   pch=19, yaxt='n', xaxt='n',
	   ylim=c(1e-5, 0.03),
	   col=c("mediumblue", "mediumblue", "mediumblue", "mediumblue", "mediumblue", "mediumblue"),
	   ylab="", 
	   main=expression(paste( "(H) anti-S2 IgM antibodies", sep="" )),
	   cex.lab=lab.size, cex.axis=axis.size, cex.main=main.size)




for(i in 1:length(line_seq_x))
{
	points(x=c(0,100), y=rep(line_seq_x[i],2), type='l', col="grey", lty="dashed")
}


axis(1, at=c(1, 3, 5), label=c("Bichat", "Cochin",  "PNC"), cex.axis=0.8*axis.size )

axis(1, at=c(2, 4, 6), label=c( "Sburg", "TRC", "EFS"), cex.axis=0.8*axis.size )



axis(2, at=c(0.00001, 0.0001, 0.001, 0.01), 
        label=c("0.00001", "0.0001", "0.001", "0.01"),
        las=2, cex.axis=axis.size )







#####################################
#####################################
##                                 ## 
##  PANELS 9                       ##
##  ROC analysis


line_seq <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0)

par(mar=c(5,7,2,1.5))
par(mgp=c(2.5, 1.25,0))


plot(x=c(0,1), y=c(0,1), type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="(I) ROC analysis",
cex.lab=3, cex.axis=2, cex.main=3)

mtext(side = 1, line = 3.4, 
cex=2, 
text="1 - specificity")

mtext(side = 2, line = 5, 
cex=2, 
text="sensitivity")

for(i in 1:length(line_seq))
{
	points(x=c(0,1), y=rep(line_seq[i],2), type='l', col="grey", lty="dashed")
	points(x=rep(line_seq[i],2), y=c(0,1), type='l', col="grey", lty="dashed")
}

points( x=1-spec_mat[,1], y=sens_mat[,1], 
    	  type='S', lwd=2, col="orangered" )


points( x=1-spec_mat[,2], y=sens_mat[,2], 
    	  type='S', lwd=2, col="gold" )


points( x=1-spec_mat[,3], y=sens_mat[,3], 
    	  type='S', lwd=2, col="limegreen" )


points( x=1-spec_mat[,4], y=sens_mat[,4], 
    	  type='S', lwd=2, col="mediumblue" )


points( x=1-spec_mat[,5], y=sens_mat[,5], 
    	  type='S', lty="dotted", lwd=2, col="orangered" )


points( x=1-spec_mat[,6], y=sens_mat[,6], 
    	  type='S', lty="dotted", lwd=2, col="gold" )


points( x=1-spec_mat[,7], y=sens_mat[,7], 
    	  type='S', lty="dotted", lwd=2, col="limegreen" )


points( x=1-spec_mat[,8], y=sens_mat[,8], 
    	  type='S', lty="dotted", lwd=2, col="mediumblue" )


legend(x="bottomright", 
cex=2, 
bg="white", box.col="white",
#fill = c("orangered", "gold", "limegreen", "mediumblue"),
#border = c("orangered", "gold",  "limegreen", "mediumblue"), 
##legend = c("Spike", "RBD", "S1", "S2") )
col = c("orangered", "gold",  "limegreen", "mediumblue", "orangered", "gold",  "limegreen", "mediumblue"),
lty=c(1,1,1,1,2,2,2,2), lwd=2,
legend = c( expression(paste( "S"^"tri", " IgG" )), expression(paste( "RBD IgG" )), expression(paste( "S1 IgG" )), expression(paste( "S2 IgG" )),
            expression(paste( "S"^"tri", " IgM" )), expression(paste( "RBD IgM" )), expression(paste( "S1 IgM" )), expression(paste( "S2 IgM" ))   )   )

axis(1, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        cex.axis=2) 

axis(2, at=c(0.0, 0.2, 0.4, 0.6, 0.8, 1), 
        labels=c("0%", "20%", "40%", "60%", "80%", "100%"), 
        las=2, cex.axis=2 ) 


dev.off()




