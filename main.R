rm(list=ls())
seed <- 1234
set.seed(seed)
## go to the directory which store the sample dataset and the codes downloaded from github
setwd("Z:/...")

## load some required R packages and helper functions.
library(expm)
library("survival")
library("MASS")
library("rootSolve")
source("HelperFunctions.R")
library(pracma)

#--------------------------------- Read in data -------------------------------#
load("dat.RData")


#-----------preprocess data to make it fit the following codes.----------------#
n=nrow(DATA1)

no_var <- length(grep("X", colnames(DATA1)))

tau = 30; # last FU time.
#----- recurrent process #1 ------#
IDR_1 = matrix(DATA.1$id1)
DeltaR_1 = matrix(DATA.1$cen1)
nobs_1 = nrow(IDR_1)

XXR_1 = as.matrix(cbind(rep(1,nobs_1),    # match with the model parameters
                        rep(DATA1$X1, times=table(DATA.1$id1)),
                        rep(DATA1$X2, times=table(DATA.1$id1)),
                        matrix(0,nrow=nobs_1,ncol=(1+no_var)))) 
ZZR_1 = matrix(c(rep(1, nobs_1), rep(0, nobs_1)), ncol=2)                       
YR_1=matrix(DATA.1$t1)
NNR_1=length(YR_1)
nnR_1=sum(DeltaR_1==1)    

### generate at risk indicator at each event time; each row correspond to a observations in a process
od <- order(YR_1[DeltaR_1==1])
OYR_1=YR_1[DeltaR_1==1][od]
indYR_1 <- 1*(matrix(rep(YR_1, nnR_1), ncol=nnR_1) >= matrix(rep(OYR_1, each=NNR_1), ncol=nnR_1))

### generate 0/1 version of the subject ID w.r.t to YR; 
### e.g., 1st row is for the 1st subject. the first three entries are 1 means the first three YR/data is for subject 1
###       2nd row us for the 2nd subject, the fourth entry is 1 means the fourth entry in YR/data is for subject 1
indIDR_1=matrix(rep(0,n*NNR_1),ncol=NNR_1)
for (i in 1:n){
  indIDR_1[i, IDR_1==i]=1
}
omega_1=sort(YR_1[DeltaR_1==1])  # ordered event times
#----- recurrent process #2 ------#
IDR_2 = matrix(DATA.2$id2)
DeltaR_2 = matrix(DATA.2$cen2)
nobs_2 = nrow(IDR_2)


XXR_2 = as.matrix(cbind(matrix(0,nrow=nobs_2,ncol=(1+no_var)),  # match with the model parameters
                        rep(1,nobs_2), 
                        rep(DATA1$X1, times=table(DATA.2$id2)),
                        rep(DATA1$X2, times=table(DATA.2$id2)))) 
ZZR_2 = matrix(c(rep(0, nobs_2), rep(1, nobs_2)), ncol=2)
YR_2 =  matrix(DATA.2$t2)


NNR_2=length(YR_2)
nnR_2=sum(DeltaR_2==1)  

od <- order(YR_2[DeltaR_2==1])
OYR_2=YR_2[DeltaR_2==1][od]
indYR_2 <- 1*(matrix(rep(YR_2, nnR_2), ncol=nnR_2) >= matrix(rep(OYR_2, each=NNR_2), ncol=nnR_2))


indIDR_2=matrix(rep(0,n*NNR_2),ncol=NNR_2)
for (i in 1:n){
  indIDR_2[i, IDR_2==i]=1
}

omega_2=sort(YR_2[DeltaR_2==1]) 

designMat_1 <- outer(YR_1[DeltaR_1==0], omega_1, FUN = ">=")
designMat_2 <- outer(YR_2[DeltaR_2==0], omega_2, FUN = ">=")

#---------------------------Set Initial Values ----------------------------#
oldxiR = matrix(c(0))                # coefficients for time-varying random effects; START arbitrarily

oldbetaR = matrix(c(1, rep(0,no_var),   # logLambda01(tau),regression coefficients for process 1,
                    1, rep(0,no_var)))  # logLambda02(tau),regression coefficients for process 2      
                                        # match with XXR_1 and XXR_2

oldalphaR_1 <- 1/nnR_1             # f_m1_1
oldalphaR_2 <- 1/nnR_2             # f_m2_2   
oldfR_1 <- rep(1/nnR_1, nnR_1)     # f_1_1, ..., f_m1_1 
oldfR_2 <- rep(1/nnR_2, nnR_2)     # f_1_2, ..., f_m2_2 
oldSigma <- matrix(c(1, 0, 0, 1), ncol=2) # covariance matrix of the random effects


## Gauss Hermite nodes and weights - two dimensional
bb <- cbind(rep(sqrt(2)*b, each=length(b)), rep(sqrt(2)*b, times=length(b)) )
ww <- as.matrix(as.vector(w%o%w))

mb <- nrow(bb)
mz <- ncol(bb)
nxi <- length(oldxiR)
nbetaR=length(oldbetaR)


### initiate the EM algorithm
epsilon <- 0.001             
maxiter <- 500              
error <- 1
iter <- 0

while(error>epsilon & iter<maxiter){
  
  #------------------------------ E step ---------------------------------# 
  Tbb <- as.matrix(expm::sqrtm(oldSigma)%*%t(bb))
  
  ptm1 <- proc.time()
  tempGamma_1 <- E_step(oldxiR=oldxiR, 
                        oldbetaR=oldbetaR, 
                        oldfR=oldfR_1, 
                        IDR=IDR_1, 
                        YR=YR_1, 
                        DeltaR=DeltaR_1, 
                        XXR=XXR_1,  
                        P_no=1, 
                        omega=omega_1, 
                        NNR=NNR_1, 
                        nnR=nnR_1, 
                        mb, 
                        Tbb,
                        indIDR=indIDR_1, 
                        designMat=designMat_1)
  tempGamma_2 <- E_step(oldxiR=oldxiR, 
                        oldbetaR=oldbetaR, 
                        oldfR=oldfR_2, 
                        IDR=IDR_2, 
                        YR=YR_2, 
                        DeltaR=DeltaR_2, 
                        XXR=XXR_2,  
                        P_no=2, 
                        omega=omega_2, 
                        NNR=NNR_2, 
                        nnR=nnR_2, 
                        mb, 
                        Tbb,
                        indIDR=indIDR_2, 
                        designMat=designMat_2)
  ptm2 <- proc.time()
  
  cat("E step:", round((ptm2-ptm1)["elapsed"]/60), "min \n")

  
  ### after multiply the design matrix indIDR, we get the contribution from each subject
  Gamma <- exp(tempGamma_1 + tempGamma_2) * matrix(rep(t(ww),n), byrow=TRUE,nrow=n)
  
  
  ### The denominator in the formula for the conditional expectation
  Gamma_sum <- rowSums(Gamma) 
  
  Gamma <- Gamma / matrix(rep(Gamma_sum, mb),ncol=mb)
  ### then Gamma can be used to calculate conditional expectation of g(b_i)
  
  
  #------------------------------ M step ---------------------------------# 
  ### calculate newSigma strictly follow formula
  newSigma <- c()
  for (i in 1:n){
    temp <- c()
    for (kkk in 1:mb){
      temp <- cbind(temp, as.numeric(Tbb[,kkk]%*%t(Tbb[,kkk])))
    }
    E_hat_b_tb <- temp%*%Gamma[i,] 
    
    newSigma <- cbind(newSigma,  E_hat_b_tb)
  }
  newSigma <- matrix(rowMeans(newSigma),ncol=2)
  
  ### next, to update other parameters, calculate the equations and their first derivatives
  ### from each recurrent process.
  
  ptm3 <- proc.time()
  newton_P1 <- newton(oldxiR=oldxiR, 
                      oldbetaR=oldbetaR, 
                      oldalphaR=oldalphaR_1, 
                      n=n,
                      IDR=IDR_1, 
                      YR=YR_1, 
                      DeltaR=DeltaR_1,
                      XXR=XXR_1,
                      P_no=1, 
                      omega=omega_1, 
                      NNR=NNR_1, 
                      nnR=nnR_1, 
                      mb=mb, 
                      Tbb=Tbb, 
                      Gamma=Gamma,
                      indIDR=indIDR_1,
                      nxi=nxi,
                      designMat=designMat_1)
  
  
  newton_P2 <- newton(oldxiR=oldxiR, 
                      oldbetaR=oldbetaR, 
                      oldalphaR=oldalphaR_2, 
                      n=n,
                      IDR=IDR_2, 
                      YR=YR_2, 
                      DeltaR=DeltaR_2,
                      XXR=XXR_2,
                      P_no=2, 
                      omega=omega_2, 
                      NNR=NNR_2, 
                      nnR=nnR_2, 
                      mb=mb, 
                      Tbb=Tbb, 
                      Gamma=Gamma,
                      indIDR=indIDR_2,
                      nxi=nxi,
                      designMat=designMat_2)
  
  ptm4 <- proc.time()
  
  cat("M step:", round((ptm4-ptm3)["elapsed"]/60), "min \n")
  
  ### the eqns to solve
  eqnxi <- newton_P1$eqnxi + newton_P2$eqnxi
  
  eqnbeta <- newton_P1$eqnbeta + newton_P2$eqnbeta
  
  eqnalpha <- c(newton_P1$eqnalpha, newton_P2$eqnalpha)
  
  
  ### the first derivatives of the eqns to be used in Newton Raphson
  eqnxi_xi <- newton_P1$eqnxi_xi + newton_P2$eqnxi_xi
  eqnxi_beta <- newton_P1$eqnxi_beta + newton_P2$eqnxi_beta
  eqnxi_alpha <- cbind(newton_P1$eqnxi_alpha, newton_P2$eqnxi_alpha)
  
  eqnbeta_xi <- newton_P1$eqnbeta_xi + newton_P2$eqnbeta_xi
  eqnbeta_beta <- newton_P1$eqnbeta_beta + newton_P2$eqnbeta_beta
  eqnbeta_alpha <- cbind(newton_P1$eqnbeta_alpha, newton_P2$eqnbeta_alpha)
  
  eqnalpha_xi <- rbind(newton_P1$eqnalpha_xi, newton_P2$eqnalpha_xi)
  eqnalpha_beta <- rbind(newton_P1$eqnalpha_beta, newton_P2$eqnalpha_beta)
  eqnalpha_alpha <- rbind(c(newton_P1$eqnalpha_alpha,0), c(0,newton_P2$eqnalpha_alpha))
  
  
  
  HH = rbind(cbind(eqnxi_xi, eqnxi_beta, eqnxi_alpha),
             cbind(eqnbeta_xi, eqnbeta_beta, eqnbeta_alpha), 
             cbind(eqnalpha_xi, eqnalpha_beta, eqnalpha_alpha))
  
  ### update unknown parameters; 
  newpar <- matrix(c(oldxiR, oldbetaR, oldalphaR_1, oldalphaR_2)) - solve(HH) %*% matrix(c(eqnxi, eqnbeta, eqnalpha))
  
  
  newxiR <- matrix(newpar[1:nxi])
  newbetaR <- matrix(newpar[(nxi+1):(nxi+nbetaR)])
  newalphaR_1 <- newpar[(nxi+nbetaR+1)]
  newalphaR_2 <- newpar[(nxi+nbetaR+2)]
  
  ptm5 <- proc.time()
  
  ### update f for baseline hazar
  fs1 <- recursive(newxiR, newbetaR, newalphaR_1, n, IDR_1, YR_1, DeltaR_1, XXR_1,  P_no=1, omega_1, NNR_1, nnR_1, mb, Tbb, Gamma)
  newfR_1 <- fs1$f   # f_j,k

  
  fs2 <- recursive(newxiR, newbetaR, newalphaR_2, n, IDR_2, YR_2, DeltaR_2, XXR_2,  P_no=2, omega_2, NNR_2, nnR_2, mb, Tbb, Gamma)
  newfR_2 <- fs2$f   # f_j,k

  ptm6 <- proc.time()
  cat("Recursive step:", round((ptm6-ptm4)["elapsed"]/60), "min \n")
  
  cat("Total time:", round((ptm6-ptm1)["elapsed"]/60), "min \n")
  
  
  ### some house-keeping
  error <- (c(abs(newxiR - oldxiR),
              abs(newbetaR-oldbetaR), abs(newfR_1-oldfR_1), abs(newfR_2-oldfR_2), 
              abs(diag(newSigma-oldSigma)), abs(newSigma[1,2]-oldSigma[1,2])))
  # print(error)
  
  error <- max(error)
  
  iter <- iter+1
  
  oldxiR=newxiR
  oldbetaR=newbetaR
  oldfR_1=newfR_1
  oldfR_2=newfR_2
  oldSigma=newSigma
  oldalphaR_1=newfR_1[nnR_1]
  oldalphaR_2=newfR_2[nnR_2]
  cat(iter, "/")
  
  ### print to monitor intermediate results
  # if (iter%%20==0){
    print("error,oldxiR,oldbetaR,            oldalphaR_1, oldalphaR_2")
    print(round(c(error, oldxiR, oldbetaR, oldalphaR_1, oldalphaR_2),3))
  # }
  

}


if (iter==maxiter){
  convergenceR = 1    # 1 means the algorithm reached the maximum iterations
}else{
  convergenceR=0
}

### record cumulative baseline hazard on a grid of 0~tau
### NOTE: should multiply the exp(longLambda), and add up lambda(t) to get the correct Lambda.

time_seq <- sort(unique(c(seq(0,tau,length.out = 100), c(0.25, 0.5, 0.75, 1)*tau)))
cumhaz_hat_1 <- stepfun(omega_1, c(0,cumsum(oldfR_1))*exp(oldbetaR[1]))(time_seq) 
cumhaz_hat_2 <- stepfun(omega_2, c(0,cumsum(oldfR_2))*exp(oldbetaR[4]))(time_seq)  
                                                                                    

res <- c(seed, convergenceR,  oldxiR,
         oldbetaR,
         oldalphaR_1,
         oldalphaR_2,
         oldSigma,
         cumhaz_hat_1,
         cumhaz_hat_2)

#------------------------------------------------------------------------------#
### calculate ESE
#------------------------------------------------------------------------------#

source("estimate_R_cov_simple.R")
source("cov_Sigma.R")
pt <- c(0.25, 0.5, 1)*30   # hard coded tau to be 30
npt <- length(pt)


if ((convergenceR==0)){
  estpar <- c(newxiR, newbetaR[-c(1,2+no_var)], 
              newSigma[1,1], newSigma[1,2], newSigma[2,2],
              newfR_1*exp(newbetaR[1]),
              newfR_2*exp(newbetaR[2+no_var]))
  newbb <- Tbb; newweight=Gamma;
  
  #----------------------------------------------------------------------------#
  nbetaonly <- length(newbetaR)-2  
  nSigma <- mz*(mz+1)/2  
  
  #-----------------------------------#

  newXR_1=matrix(XXR_1[,-c(1,2+no_var)], ncol=ncol(XXR_1)-2)  # only those covariates
  
  newlambdaR_1=newfR_1*exp(newbetaR[1])                        # back to the original scale, the jump size at each event time
  
  newLambdaR_1=indYR_1%*%newlambdaR_1                          # this way, the unordered cumhaz is computed, with redundant values.
  # this newLambdaR is the cumulative hazard correspond
  # to the observation times in YR  
  
  
  newXR_2=matrix(XXR_2[,-c(1,2+no_var)], ncol=ncol(XXR_2)-2)
  
  newlambdaR_2=newfR_2*exp(newbetaR[2+no_var])  
  
  newLambdaR_2=indYR_2%*%newlambdaR_2 
  #-----------------------------------#
  
  npar <- length(as.numeric(newxiR)) + nbetaonly + nSigma   
  
  
  #--------------#
  ### define the second derivative matrix; store the info from each subject computed from the loop that follows.
  ### There are 16 matrices, but since some are equivalent except for dimensions, we only need to compute 11 matrices.
  Dbeta2R=matrix(0, nbetaonly,nbetaonly);    # sum i # second derivatives wrt beta
  Dbeta2R_1=matrix(0,nbetaonly,nbetaonly);           # the part from the first process
  Dbeta2R_2=matrix(0,nbetaonly,nbetaonly);           # the part from the second process
  
  DbetaLambdaR_1=matrix(0,nbetaonly, nnR_1);  # sum i       # partial beta partial Lambda
  DbetaLambdaR_2=matrix(0,nbetaonly, nnR_2); 
  
  DLambda2R_1=matrix(0,nnR_1,nnR_1);          # sum i       # second derivatives wrt Lambda
  DLambda2R_2=matrix(0,nnR_2,nnR_2);
  
  
  Dxi2R=matrix(0, length(as.numeric(newxiR)), length(as.numeric(newxiR)))   # sum i # second derivative wrt xi
  Dxi2R_1=matrix(0, length(as.numeric(newxiR)), length(as.numeric(newxiR))) # the part from the first process
  Dxi2R_2=matrix(0, length(as.numeric(newxiR)), length(as.numeric(newxiR))) # the part from the second process
  
  DxibetaR=matrix(0, length(as.numeric(newxiR)), nbetaonly)  # sum i # partial xi partial beta
  DxibetaR_1=matrix(0, length(as.numeric(newxiR)), nbetaonly) # the part from the first process
  DxibetaR_2=matrix(0, length(as.numeric(newxiR)), nbetaonly) # the part from the second process
  
  DxiLambdaR_1 = matrix(0, length(as.numeric(newxiR)), nnR_1)  # sum i # partial xi partial Lambda
  DxiLambdaR_2 = matrix(0, length(as.numeric(newxiR)), nnR_2)
  
  #--------------#
  DL=matrix(0,npar+nnR_1+nnR_2,n);                           # E_hat{l'}
  DL2=matrix(0,npar+nnR_1+nnR_2,npar+nnR_1+nnR_2);   # sum i # E_hat{l't(l')}
  Dgamma2=matrix(0,nSigma,nSigma);                   # sum i # second derivatives wrt the parameters in the covariance matrix Sigma
  
  for (i in 1:n){
    cat(paste(i, "/"))
    # print(Dbeta2R)
    est_cov1 <- estimate_R_cov_simple(newbetaR=newbetaR,
                                      newxiR=newxiR,
                                      nxi=length(as.numeric(newxiR)),
                                      no_var=no_var,
                                      omega=omega_1,  # omega_1 and OYR_1 are same.
                                      ZZR=ZZR_1,
                                      newXR=newXR_1,
                                      newLambdaR=newLambdaR_1,
                                      newlambdaR=newlambdaR_1,
                                      YR=YR_1,
                                      DeltaR=DeltaR_1,
                                      nnR=nnR_1,
                                      IDR=IDR_1,
                                      OYR=OYR_1,
                                      
                                      Dbeta2R=Dbeta2R_1,        # Note: only feed the part from the 1st process 
                                      DbetaLambdaR=DbetaLambdaR_1,
                                      DLambda2R=DLambda2R_1,
                                      DxibetaR=DxibetaR_1,      # Note: only feed the part from the 1st process 
                                      DxiLambdaR=DxiLambdaR_1,
                                      Dxi2R = Dxi2R_1,          # Note: only feed the part from the 1st process 
                                      
                                      nbetaonly,mb,newbb,newweight,i)
    
    
    est_cov2 <- estimate_R_cov_simple(newbetaR=newbetaR,
                                      newxiR=newxiR,
                                      nxi=length(as.numeric(newxiR)),
                                      no_var=no_var,
                                      omega=omega_2,
                                      ZZR=ZZR_2,
                                      newXR=newXR_2,
                                      newLambdaR=newLambdaR_2,
                                      newlambdaR=newlambdaR_2,
                                      YR=YR_2,
                                      DeltaR=DeltaR_2,
                                      nnR=nnR_2,
                                      IDR=IDR_2,
                                      OYR=OYR_2,
                                      
                                      Dbeta2R=Dbeta2R_2,
                                      DbetaLambdaR=DbetaLambdaR_2,
                                      DLambda2R=DLambda2R_2,
                                      DxibetaR=DxibetaR_2,
                                      DxiLambdaR=DxiLambdaR_2,
                                      Dxi2R = Dxi2R_2,
                                      
                                      nbetaonly,mb,newbb,newweight,i)
    
    
    
    
    ### assign the intermediate results from the previous 1 to i subjects to Dbeta2R_1 and Dbeta2R_2.
    ### same for DbetaLambda_1&2
    ### the following updates the values and will be used in this for loop.
    DbetaLambdaR_1 <- est_cov1$DbetaLambdaR
    Dbeta2R_1 <- est_cov1$Dbeta2R           
    DLambda2R_1 <- est_cov1$DLambda2R
    DxibetaR_1 <- est_cov1$DxibetaR
    DxiLambdaR_1 <- est_cov1$DxiLambdaR
    Dxi2R_1 <- est_cov1$Dxi2R
    
    DbetaLambdaR_2 <- est_cov2$DbetaLambdaR
    Dbeta2R_2 <- est_cov2$Dbeta2R
    DLambda2R_2 <- est_cov2$DLambda2R
    DxibetaR_2 <- est_cov2$DxibetaR
    DxiLambdaR_2 <- est_cov2$DxiLambdaR
    Dxi2R_2 <- est_cov2$Dxi2R
    
    #----------#
    
    tempDbetaR <- est_cov1$tempDbetaR + est_cov2$tempDbetaR   # first derivative wrt beta, before multiplying the weights
    # this step corresponds to sum over k
    
    tempDxiR <- est_cov1$tempDxiR + est_cov2$tempDxiR         # first derivative wrt xi, before multiplying the weights
    # this step corresponds to sum over k
    
    cov_sig = cov_Sigma(newSigma,newbb,newweight,i)  ## since the part relevant to Sigma is irrelevant to other parameters. 
    ## We can just use the same codes.
    
    tempDgamma = cov_sig$scoreSig_i    # 3*200; first derivative wrt Sigma, before multiplying the weights, contributed by subject i.
    Dgamma2=Dgamma2 + cov_sig$dscoreSig_i   # 3*3; second derivative wrt Sigma
    
    ###
    temp <- rbind(tempDxiR, tempDbetaR, tempDgamma, est_cov1$tempDLambdaR, est_cov2$tempDLambdaR)
    
    DL[,i] <- temp%*%matrix(newweight[i,])           
    ## So, each column in DL stores the first derivatives after multiplying the weights.
    ## Based on Louis's paper, DL is S*(y,theta) should be zero at the last iteration of the EM procedure.Check.
    
    DL2 = DL2 + temp%*%t(t(matrix(rep(newweight[i,], npar+nnR_1+nnR_2), ncol=npar+nnR_1+nnR_2))*temp); 
    ## E(first derivative %*% t(first derivative)), and incorporating the weights
    
    
  }
  
  Dbeta2R <-  Dbeta2R_1 + Dbeta2R_2  
  
  Dxi2R <- Dxi2R_1 + Dxi2R_2
  
  DxibetaR <- DxibetaR_1 + DxibetaR_2
  
  
  ### E_hat{l''}; the order in the matrix is (1) xi (2) beta (3) Sigma (4) Lambda_1 (5) Lambda_2
  tempD2 = rbind(
    cbind(Dxi2R, DxibetaR, matrix(0, length(as.numeric(newxiR)), nSigma), DxiLambdaR_1, DxiLambdaR_2),
    cbind(t(DxibetaR), Dbeta2R, matrix(0, nbetaonly, nSigma), DbetaLambdaR_1, DbetaLambdaR_2),
    cbind(matrix(0, nSigma, length(as.numeric(newxiR))+nbetaonly), Dgamma2, matrix(0, nSigma, nnR_1), matrix(0, nSigma, nnR_2)),
    cbind(t(DxiLambdaR_1), t(DbetaLambdaR_1), matrix(0, nnR_1, nSigma), DLambda2R_1, matrix(0, nnR_1, nnR_2)),
    cbind(t(DxiLambdaR_2), t(DbetaLambdaR_2), matrix(0, nnR_2, nSigma), matrix(0, nnR_2, nnR_1), DLambda2R_2)
  )
  
  cov = -tempD2 - (DL2 - DL%*%t(DL)) 

  cov = pracma::inv(cov)   # Checked with time-constant version of cov -> Correct.
  
  
  #----------------------------------------------------------------------------#
  
  
  
  # Note that we have estimates for the jump sizes Lambda{}, not Lambda(). So, apply multivariate delta method.
  MX <- rbind(cbind(diag(rep(1, nxi + nbetaonly + nSigma)), matrix(0, nxi + nbetaonly + nSigma, nnR_1+nnR_2)),
              cbind(matrix(0, 3, nxi + nbetaonly + nSigma), matrix(rep(pt,nnR_1),ncol=nnR_1) >= t(matrix(rep(omega_1,npt),nrow=nnR_1)),
                    matrix(0, npt, nnR_2)),
              cbind(matrix(0, 3, nxi + nbetaonly + nSigma), matrix(0, npt, nnR_1), 
                    matrix(rep(pt,nnR_2),ncol=nnR_2) >= t(matrix(rep(omega_2,npt),nrow=nnR_2))  )
  )
  VAR <- diag(MX%*%cov%*%t(MX))
  
  
} else {
  
  VAR <- rep(NA,   nxi + nbetaonly + nSigma + 1 + npt*2)
  
}
res <- c(res, VAR)
names(res) <- c("seed",  "convergence",  "beta",
                  "logLambda_tau_1", "alpha1_1", "alpha1_2",
                  "logLambda_tau_2", "alpha2_1", "alpha2_2",
                  "f_m1_1", "f_m2_2", "sigma_sq_1", "Cov", "Cov_rep", "sigma_sq_2",
                  paste0("Lambda01",1:103),
                  paste0("Lambda02",1:103),
                  "SE_beta" ,   
                  "SE_alpha11",  "SE_alpha12",  "SE_alpha21",  "SE_alpha22",    
                  "SE_sigma11", "SE_sigma12", "SE_sigma21",
                  "SE_Lambda11", "SE_Lambda12", "SE_Lambda13", 
                  "SE_Lambda21", "SE_Lambda22", "SE_Lambda23")



