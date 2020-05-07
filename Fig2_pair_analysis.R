library(binom)


load("C:\\U\\CoronaVirus\\SeroSigPaper\\ANALYSIS\\SeroData\\Data_Process.RData")


set.seed(1234)

N_tree <- 50000

N_negative <- N_trc + N_pnc + N_efs

N_positive <- N_bichat + N_stras + N_cochin

set.seed(1234)

library(MASS)
library(ROCR)
library(randomForest)
library(pROC)

status <- c( rep("pos", N_bichat), rep("pos", N_stras), rep("pos", N_cochin), rep("neg", N_trc), rep("neg", N_pnc), rep("neg", N_efs) )
status <- as.factor(status)

AB_data <- cbind( c(BICHAT$Spike_IPP_IgG_dil, STRAS$Spike_IPP_IgG_dil, COCHIN$Spike_IPP_IgG_dil, TRC$Spike_IPP_IgG_dil, PNC$Spike_IPP_IgG_dil, EFS$Spike_IPP_IgG_dil),
                  c(BICHAT$RBD_IPP_IgG_dil,   STRAS$RBD_IPP_IgG_dil,   COCHIN$RBD_IPP_IgG_dil,   TRC$RBD_IPP_IgG_dil,   PNC$RBD_IPP_IgG_dil,   EFS$RBD_IPP_IgG_dil),
                  c(BICHAT$S1_NA_IgG_dil,     STRAS$S1_NA_IgG_dil,     COCHIN$S1_NA_IgG_dil,     TRC$S1_NA_IgG_dil,     PNC$S1_NA_IgG_dil,     EFS$S1_NA_IgG_dil),
                  c(BICHAT$S2_NA_IgG_dil,     STRAS$S2_NA_IgG_dil,     COCHIN$S2_NA_IgG_dil,     TRC$S2_NA_IgG_dil,     PNC$S2_NA_IgG_dil,     EFS$S2_NA_IgG_dil),
                  c(BICHAT$Spike_IPP_IgM_dil, STRAS$Spike_IPP_IgM_dil, COCHIN$Spike_IPP_IgM_dil, TRC$Spike_IPP_IgM_dil, PNC$Spike_IPP_IgM_dil, EFS$Spike_IPP_IgM_dil),
                  c(BICHAT$RBD_IPP_IgM_dil,   STRAS$RBD_IPP_IgM_dil,   COCHIN$RBD_IPP_IgM_dil,   TRC$RBD_IPP_IgM_dil,   PNC$RBD_IPP_IgM_dil,   EFS$RBD_IPP_IgM_dil),
                  c(BICHAT$S1_NA_IgM_dil,     STRAS$S1_NA_IgM_dil,     COCHIN$S1_NA_IgM_dil,     TRC$S1_NA_IgM_dil,     PNC$S1_NA_IgM_dil,     EFS$S1_NA_IgM_dil),
                  c(BICHAT$S2_NA_IgM_dil,     STRAS$S2_NA_IgM_dil,     COCHIN$S2_NA_IgM_dil,     TRC$S2_NA_IgM_dil,     PNC$S2_NA_IgM_dil,     EFS$S2_NA_IgM_dil) )  

colnames(AB_data) <- c( "Spike_IPP_IgG","RBD_IPP_IgG", "S1_NA_IgG", "S2_NA_IgG", "Spike_IPP_IgM", "RBD_IPP_IgM", "S1_NA_IgM", "S2_NA_IgM")

AB_data <- as.data.frame(AB_data)







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




####################################
## Variable Importance

RF_all = randomForest( status ~ ., data=AB_data[,c(1,2,3,4)], 
                                      importance=TRUE, ntree=N_tree)

rf.roc_all <- roc(status, RF_all$votes[,2])


varImpPlot( RF_all )



####################################
## 2 antigens: Spike & RBD

RF_Spike_RBD = randomForest( status ~ ., data=AB_data[,c(1,2)], 
                                      importance=TRUE, ntree=N_tree)

rf.roc_Spike_RBD <- roc(status, RF_Spike_RBD$votes[,2])



####################################
## 2 antigens: Spike & S1

RF_Spike_S1 = randomForest( status ~ ., data=AB_data[,c(1,3)], 
                                      importance=TRUE, ntree=N_tree)

rf.roc_Spike_S1 <- roc(status, RF_Spike_S1$votes[,2])



####################################
## 2 antigens: Spike & S2

RF_Spike_S2 = randomForest( status ~ ., data=AB_data[,c(1,4)], 
                                      importance=TRUE, ntree=N_tree)

rf.roc_Spike_S2 <- roc(status, RF_Spike_S2$votes[,2])


####################################
## 2 antigens: RBD & S1

RF_RBD_S1 = randomForest( status ~ ., data=AB_data[,c(2,3)], 
                                      importance=TRUE, ntree=N_tree)

rf.roc_RBD_S1 <- roc(status, RF_RBD_S1$votes[,2])



####################################
## 2 antigens: RBD & S2

RF_RBD_S2 = randomForest( status ~ ., data=AB_data[,c(2,4)], 
                                      importance=TRUE, ntree=N_tree)

rf.roc_RBD_S2 <- roc(status, RF_RBD_S2$votes[,2])


####################################
## 2 antigens: S1 & S2

RF_S1_S2 = randomForest( status ~ ., data=AB_data[,c(3,4)], 
                                      importance=TRUE, ntree=N_tree)

rf.roc_S1_S2 <- roc(status, RF_S1_S2$votes[,2])



####################################
## 3 antigens: Spike, RBD & S1

RF_Spike_RBD_S1 = randomForest( status ~ ., data=AB_data[,c(1,2,3)], 
                                      importance=TRUE, ntree=N_tree)

rf.roc_Spike_RBD_S1 <- roc(status, RF_Spike_RBD_S1$votes[,2])



####################################
## 3 antigens: Spike, RBD & S2

RF_Spike_RBD_S2 = randomForest( status ~ ., data=AB_data[,c(1,2,4)], 
                                      importance=TRUE, ntree=N_tree)

rf.roc_Spike_RBD_S2 <- roc(status, RF_Spike_RBD_S2$votes[,2])




####################################
## 3 antigens: Spike, S1 & S2

RF_Spike_S1_S2 = randomForest( status ~ ., data=AB_data[,c(1,3,4)], 
                                      importance=TRUE, ntree=N_tree)

rf.roc_Spike_S1_S2 <- roc(status, RF_Spike_S1_S2$votes[,2])


####################################
## 3 antigens: RBD, S1 & S2

RF_RBD_S1_S2 = randomForest( status ~ ., data=AB_data[,c(2,3,4)], 
                                      importance=TRUE, ntree=N_tree)

rf.roc_RBD_S1_S2 <- roc(status, RF_RBD_S1_S2$votes[,2])









sens_target <- matrix(NA, nrow=15, ncol=3)
rownames(sens_target) <- c("Spike", "RBD", "S1", "S2",
                           "Spike_RBD", "Spike_S1", "Spike_S2", "RBD_S1", "RBD_S2", "S1_S2",
                           "Spike_RBD_S1", "Spike_RBD_S2", "Spike_S1_S2", "RBD_S1_S2", 
	                     "all")
colnames(sens_target) <- c("spec_99", "sens_spec", "sens_99")

sens_target_lwr <- sens_target
sens_target_upr <- sens_target

spec_target <- sens_target
spec_target_lwr <- sens_target
spec_target_upr <- sens_target



for(i in 5:15)
{
	if( i == 5 )
	{
		RFX <- rf.roc_Spike_RBD
	}

	if( i == 6 )
	{
		RFX <- rf.roc_Spike_S1
	}		

	if( i == 7 )
	{
		RFX <- rf.roc_Spike_S2
	}		

	if( i == 8 )
	{
		RFX <- rf.roc_RBD_S1
	}		

	if( i == 9 )
	{
		RFX <- rf.roc_RBD_S2
	}		

	if( i == 10 )
	{
		RFX <- rf.roc_S1_S2
	}		

	if( i == 11 )
	{
		RFX <- rf.roc_Spike_RBD_S1
	}

	if( i == 12 )
	{
		RFX <- rf.roc_Spike_RBD_S2
	}

	if( i == 13 )
	{
		RFX <- rf.roc_Spike_S1_S2
	}

	if( i == 14 )
	{
		RFX <- rf.roc_RBD_S1_S2
	}

	if( i == 15 )
	{
		RFX <- rf.roc_all
	}

	################################
	## High specificity target

	sens_target[i,1]     <- binom.confint( RFX$sensitivities[min(which(RFX$specificities > 0.99))]*N_positive, N_positive, method="wilson")[1,4]
	sens_target_lwr[i,1] <- binom.confint( RFX$sensitivities[min(which(RFX$specificities > 0.99))]*N_positive, N_positive, method="wilson")[1,5]
	sens_target_upr[i,1] <- binom.confint( RFX$sensitivities[min(which(RFX$specificities > 0.99))]*N_positive, N_positive, method="wilson")[1,6]

	spec_target[i,1]     <- binom.confint( RFX$specificities[min(which(RFX$specificities > 0.99))]*N_negative, N_negative, method="wilson")[1,4]
	spec_target_lwr[i,1] <- binom.confint( RFX$specificities[min(which(RFX$specificities > 0.99))]*N_negative, N_negative, method="wilson")[1,5]
	spec_target_upr[i,1] <- binom.confint( RFX$specificities[min(which(RFX$specificities > 0.99))]*N_negative, N_negative, method="wilson")[1,6]


	################################
	## Balanced sensitivty and specificity target

	sens_target[i,2]     <- binom.confint( RFX$sensitivities[which.min(abs(RFX$sensitivities - RFX$specificities))]*N_positive, N_positive, method="wilson")[1,4]
	sens_target_lwr[i,2] <- binom.confint( RFX$sensitivities[which.min(abs(RFX$sensitivities - RFX$specificities))]*N_positive, N_positive, method="wilson")[1,5]
	sens_target_upr[i,2] <- binom.confint( RFX$sensitivities[which.min(abs(RFX$sensitivities - RFX$specificities))]*N_positive, N_positive, method="wilson")[1,6]

	spec_target[i,2]     <- binom.confint( RFX$specificities[which.min(abs(RFX$sensitivities - RFX$specificities))]*N_negative, N_negative, method="wilson")[1,4]
	spec_target_lwr[i,2] <- binom.confint( RFX$specificities[which.min(abs(RFX$sensitivities - RFX$specificities))]*N_negative, N_negative, method="wilson")[1,5]
	spec_target_upr[i,2] <- binom.confint( RFX$specificities[which.min(abs(RFX$sensitivities - RFX$specificities))]*N_negative, N_negative, method="wilson")[1,6]


	################################
	## High sensitivity target

	sens_target[i,3]     <- binom.confint( RFX$sensitivities[max(which(RFX$sensitivities > 0.99))]*N_positive, N_positive, method="wilson")[1,4]
	sens_target_lwr[i,3] <- binom.confint( RFX$sensitivities[max(which(RFX$sensitivities > 0.99))]*N_positive, N_positive, method="wilson")[1,5]
	sens_target_upr[i,3] <- binom.confint( RFX$sensitivities[max(which(RFX$sensitivities > 0.99))]*N_positive, N_positive, method="wilson")[1,6]

	spec_target[i,3]     <- binom.confint( RFX$specificities[max(which(RFX$sensitivities > 0.99))]*N_negative, N_negative, method="wilson")[1,4]
	spec_target_lwr[i,3] <- binom.confint( RFX$specificities[max(which(RFX$sensitivities > 0.99))]*N_negative, N_negative, method="wilson")[1,5]
	spec_target_upr[i,3] <- binom.confint( RFX$specificities[max(which(RFX$sensitivities > 0.99))]*N_negative, N_negative, method="wilson")[1,6]

}



spec_target     <- round(100*spec_target,1)
spec_target_lwr <- round(100*spec_target_lwr,1)
spec_target_upr <- round(100*spec_target_upr,1)

sens_target     <- round(100*sens_target,1)
sens_target_lwr <- round(100*sens_target_lwr,1)
sens_target_upr <- round(100*sens_target_upr,1)




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



tiff(file="Fig2_pair_analysis.tif", width=40, height=25, units="cm", res=500)

lay.mat <- rbind( c( 1, 2, 6 ),
                  c( 3, 4, 6 ),
                  c( 5, 5, 6 ) )
layout(lay.mat, heights=c(10,10,2), widths=c(1.2,1.2,2))
layout.show(6)

par(mar=c(5,7.5,3,1))
par(mgp=c(3, 1, 0))

point.size = 0.75
lab.size   = 2
axis.size  = 1.5
main.size  = 2


###################################
## Panel 1: Spike_IPP 

plot( x = BICHAT$RBD_IPP_IgG_dil, y = BICHAT$Spike_IPP_IgG_dil,
      xlim=c(1e-5, 0.03), ylim=c(1e-5, 0.03), log="xy",
      pch=19, col="yellowgreen",
	yaxt='n', xaxt='n', bty='n', 
	xlab=expression(paste( "anti-RBD IgG dilution", sep="" )), 
	ylab="", 
	main=expression(paste( "(A) anti-S"^"tri", " IgG vs anti-RBD IgG", sep="" )),
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

mtext(side = 2, line = 4.5, 
cex=axis.size, 
text=expression(paste( "anti-S"^"tri", " IgG dilution", sep="" )))

points( x = STRAS$RBD_IPP_IgG_dil, y = STRAS$Spike_IPP_IgG_dil,
        pch=19, col="sienna1" )

points( x = COCHIN$RBD_IPP_IgG_dil, y = COCHIN$Spike_IPP_IgG_dil,
        pch=19, col="darkmagenta" )

points( x = TRC$RBD_IPP_IgG_dil, y = TRC$Spike_IPP_IgG_dil,
        pch=19, col="royalblue" )

points( x = PNC$RBD_IPP_IgG_dil, y = PNC$Spike_IPP_IgG_dil,
        pch=19, col="cornflowerblue" )

points( x = EFS$RBD_IPP_IgG_dil, y = EFS$Spike_IPP_IgG_dil,
        pch=19, col="dodgerblue" )

axis(1, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
        cex.axis=axis.size )

axis(2, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
        las=2, cex.axis=axis.size )


###################################
## Panel 2: RBD_IPP 

plot( x = BICHAT$S1_NA_IgG_dil, y = BICHAT$Spike_IPP_IgG_dil,
      xlim=c(1e-5, 0.03), ylim=c(1e-5, 0.03), log="xy",
      pch=19, col="yellowgreen",
	yaxt='n', xaxt='n', bty='n', 
	xlab=expression(paste( "anti-S1 IgG dilution", sep="" )), 
	ylab="", 
	main=expression(paste( "(B) anti-S"^"tri", " IgG vs anti-S1 IgG", sep="" )),
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

mtext(side = 2, line = 4.5, 
cex=axis.size, 
text=expression(paste( "anti-S"^"tri", " IgG dilution", sep="" )))


points( x = STRAS$S1_NA_IgG_dil, y = STRAS$Spike_IPP_IgG_dil,
        pch=19, col="sienna1" )

points( x = COCHIN$S1_NA_IgG_dil, y = COCHIN$Spike_IPP_IgG_dil,
        pch=19, col="darkmagenta" )

points( x = TRC$S1_NA_IgG_dil, y = TRC$Spike_IPP_IgG_dil,
        pch=19, col="royalblue" )

points( x = PNC$S1_NA_IgG_dil, y = PNC$Spike_IPP_IgG_dil,
        pch=19, col="cornflowerblue" )

points( x = EFS$S1_NA_IgG_dil, y = EFS$Spike_IPP_IgG_dil,
        pch=19, col="dodgerblue" )


axis(1, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
        cex.axis=axis.size )

axis(2, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
        las=2, cex.axis=axis.size )



###################################
## Panel 3:  S1_NA

plot( x = BICHAT$S2_NA_IgG_dil, y = BICHAT$Spike_IPP_IgG_dil,
      xlim=c(1e-5, 0.03), ylim=c(1e-5, 0.03), log="xy",
      pch=19, col="yellowgreen",
	yaxt='n', xaxt='n', bty='n', 
	xlab=expression(paste( "anti-S2 IgG dilution", sep="" )), 
	ylab="", 
	main=expression(paste( "(C) anti-S"^"tri", " IgG vs anti-S2 IgG", sep="" )),
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

mtext(side = 2, line = 4.5, 
cex=axis.size, 
text=expression(paste( "anti-S"^"tri", " IgG dilution", sep="" )))


points( x = STRAS$S2_NA_IgG_dil, y = STRAS$Spike_IPP_IgG_dil,
        pch=19, col="sienna1" )

points( x = COCHIN$S2_NA_IgG_dil, y = COCHIN$Spike_IPP_IgG_dil,
        pch=19, col="darkmagenta" )

points( x = TRC$S2_NA_IgG_dil, y = TRC$Spike_IPP_IgG_dil,
        pch=19, col="royalblue" )

points( x = PNC$S2_NA_IgG_dil, y = PNC$Spike_IPP_IgG_dil,
        pch=19, col="cornflowerblue" )

points( x = EFS$S2_NA_IgG_dil, y = EFS$Spike_IPP_IgG_dil,
        pch=19, col="dodgerblue" )

axis(1, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
        cex.axis=axis.size )

axis(2, at=c(0.00001, 0.0001, 0.001, 0.01, 0.03), 
        label=c("0.00001", "0.0001", "0.001", "0.01", ""),
        las=2, cex.axis=axis.size )



###################################
## Panel 4: S2_NA

plot( x = BICHAT$Spike_IPP_IgM_dil, y = BICHAT$Spike_IPP_IgG_dil,
      xlim=c(1e-5, 0.03), ylim=c(1e-5, 0.03), log="xy",
      pch=19, col="yellowgreen",
	yaxt='n', xaxt='n', bty='n', 
	xlab=expression(paste( "anti-S"^"tri", " IgM dilution", sep="" )), 
	ylab="", 
	main=expression(paste( "(D) anti-S"^"tri", " IgG vs anti-S"^"tri", " IgM", sep="" )),
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

mtext(side = 2, line = 4.5, 
cex=axis.size, 
text=expression(paste( "anti-S"^"tri", " IgG dilution", sep="" )))


points( x = STRAS$Spike_IPP_IgM_dil, y = STRAS$Spike_IPP_IgG_dil,
        pch=19, col="sienna1" )

points( x = COCHIN$Spike_IPP_IgM_dil, y = COCHIN$Spike_IPP_IgG_dil,
        pch=19, col="darkmagenta" )

points( x = TRC$Spike_IPP_IgM_dil, y = TRC$Spike_IPP_IgG_dil,
        pch=19, col="royalblue" )

points( x = PNC$Spike_IPP_IgM_dil, y = PNC$Spike_IPP_IgG_dil,
        pch=19, col="cornflowerblue" )


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
       ncol=3, cex=2.5, bty="n" )


#####################################
#####################################
##                                 ## 
##  PANELS 9                       ##
##  ROC analysis


PP <- 10


line_seq_x <- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)
line_seq_y <- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)



par(mar=c(5,7,3,1.5))
par(mgp=c(2.5, 1.25,0))


plot(x=100, y=100, type='l', lty="dashed",
xlim=c(0,1.002), ylim=c(0,1.002),
#xaxs='i', yaxs='i', 
xaxt='n', yaxt='n', bty='n',
xlab="", ylab="",
main="(E) ROC analysis",
cex.lab=3, cex.axis=2, cex.main=3)

mtext(side = 1, line = 3.4, 
cex=2, 
text="1 - specificity")

mtext(side = 2, line = 4.5, 
cex=2, 
text="sensitivity")

for(i in 1:length(line_seq_y))
{
	points(x=c(0,1), y=rep(line_seq_y[i]^PP,2), type='l', col="grey", lty="dashed")
}

for(i in 1:length(line_seq_x))
{
	points(x=rep(1 - line_seq_x[i]^PP,2), y=c(0,1), type='l', col="grey", lty="dashed")
}


points( x=1-spec_mat[,1]^PP, y=sens_mat[,1]^PP, 
    	  type='S', lwd=2, col="orangered" )


points( x=1-rf.roc_Spike_RBD$specificities^PP, y=rf.roc_Spike_RBD$sensitivities^PP, 
    	  type='S', lwd=2, col="gold" )

points( x=1-rf.roc_Spike_RBD_S1$specificities^PP, y=rf.roc_Spike_RBD_S1$sensitivities^PP, 
    	  type='S', lwd=2, col="mediumblue" )



points( x=1-rf.roc_all$specificities^PP, y=rf.roc_all$sensitivities^PP, 
    	  type='S', lwd=2, col="green4" )


legend(x="bottomright", 
cex=2, 
bg="white", box.col="white",
fill = c("orangered", "gold", "mediumblue",  "green4"),
border = c("orangered", "gold", "mediumblue",  "green4"), 
legend = c( expression(paste( "S"^"tri", " IgG", sep="" )),
            expression(paste( "S"^"tri", " IgG + RBD IgG", sep="" )),
            expression(paste( "S"^"tri", " IgG + RBD IgG + S1 IgG", sep="" )),  
            expression(paste( "S"^"tri", " IgG + RBD IgG + S1 IgG + S2 IgG", sep="" )) ) )  

axis(1, at=1- c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=rev(c("0%", "1%", "5%", "10%", "20%", "100%")), 
        cex.axis=2) 

axis(2, at=c(0.0, 0.8, 0.9, 0.95, 0.99, 1)^PP, 
        labels=c("0%", "80%", "90%", "95%", "99%", "100%"), 
        las=2, cex.axis=2 ) 


dev.off()






