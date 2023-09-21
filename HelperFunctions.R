### Acknowledgement:
### The following codes borrowed heavily from the MATLAB codes from
### (1) Zhu L, Sun J, Tong X et al. Regression analysis of multivariate recurrent 
###     event data with a dependent terminal event. Lifetime Data Analysis 2010; 
###     16(4): 478–490.
### (2) Zeng D and Lin D. Semiparametric transformation models with random 
###     effects for recurrent events. Journal of the American Statistical 
###     Association 2007; 102(477): 167–180

### define the time-varying h(t) function; its functional form is predetermined.
### change it to fit your data; note that the first component in h(t) was fixed at 1
### for identifiability.

poly_fun <- function(xi, t){
  ### xi is a vector
  ### t is a scalar
  ### return a scalar
  bases <- matrix(c(1*(t<=12)))     
  
  return(as.numeric(t(xi)%*%bases))
}

poly_bases <- function(t) {
  bases <- matrix(c(1*(t<=12)))    
  
  return(bases)
}

### GaussHermiteNodeWeight
b=c(.24534070830090124990383653063362,
    -.24534070830090124990383653063362,
    .73747372854539435870560514425210,
    -.73747372854539435870560514425210,
    1.2340762153953230078858183469594,
    -1.2340762153953230078858183469594,
    1.7385377121165862067808656621364,
    -1.7385377121165862067808656621364,
    2.2549740020892755230823333447346,
    -2.2549740020892755230823333447346,
    2.7888060584281304805250337564032,
    -2.7888060584281304805250337564032,
    3.3478545673832163269149245229964,
    -3.3478545673832163269149245229964,
    3.9447640401156252103756288005244,
    -3.9447640401156252103756288005244,
    4.6036824495507442730776752489784,
    -4.6036824495507442730776752489784,
    5.3874808900112328620169004106810,
    -5.3874808900112328620169004106810)

w=c(    .46224366960061003033827572574491,
        .46224366960061003033827572574491,
        .28667550536283409293535337283309,
        .28667550536283409293535337283309,
        .10901720602002330602538760473843,
        .10901720602002330602538760473843,
        .24810520887463607698642874311139e-1,
        .24810520887463607698642874311139e-1,
        .32437733422378614159636885549073e-2,
        .32437733422378614159636885549073e-2,
        .22833863601635393795824113589902e-3,
        .22833863601635393795824113589902e-3,
        .78025564785320626929735083989987e-5,
        .78025564785320626929735083989987e-5,
        .10860693707692815546422802702487e-6,
        .10860693707692815546422802702487e-6,
        .43993409922731799891344886730418e-9,
        .43993409922731799891344886730418e-9,
        .22293936455341510064614411221584e-12,
        .22293936455341510064614411221584e-12)



E_step <- function(oldxiR, oldbetaR, oldfR, IDR, YR, DeltaR, XXR,  P_no, omega, NNR, nnR, mb, Tbb,  indIDR, designMat){

  LHS <- matrix(1+sapply(YR, function(x) poly_fun(oldxiR,x)))
  
  LHS[DeltaR==0] <- 0
  
  LHS <- LHS%*%Tbb[P_no,] 
  
  LHS <- indIDR %*% LHS 
  
  RHS <- matrix(rep(matrix(oldfR), mb), ncol=mb) * exp(matrix(1+sapply(omega, function(x) poly_fun(oldxiR,x)))%*%Tbb[P_no,])
  
  RHS <- designMat %*% RHS 
  
  RHS <- RHS * exp(matrix( rep(XXR[DeltaR==0,]%*%oldbetaR,mb), ncol=mb))
  
  tempGamma <- LHS - RHS
  
  return(tempGamma)
}


recursive <- function(xiR, betaR, alphaR_1, n, IDR_1, YR_1, DeltaR_1, XXR_1,  P_no=1, omega_1, NNR_1, nnR_1, mb, Tbb, Gamma){
  ### omega is the ordered event times
  ### returns f = f_j, j=1,..mk
  ### returns f1 = [partial f_1/partial alpha    (partial f_1/partial beta)^T]
  ###              ...
  ###              [partial f_mk/partial alpha   (partial f_mk/partial beta)^T]
  
  nbetaR <- length(betaR)
  nxiR <- length(xiR)
  
  
  ### initialize
  f <- rep(NA, nnR_1)
  f1 <- matrix(rep(NA, nnR_1*(nxiR+nbetaR+1)), nrow=nnR_1)  # store summands of the partial derivatives of f_m
  # w.r.t xi, beta, alpha_k
  
  
  ### initial values
  f[nnR_1] <- alphaR_1
  f1[nnR_1, 1:nxiR] <- 0             # for partial f_mk,k/partial xi
  f1[nnR_1, (nxiR+1):(nxiR+nbetaR)] <- 0   # for partial f_mk,k/partial beta t
  f1[nnR_1, (nxiR+nbetaR+1)] <- 1       # for partial fm / partial alpha_k
  
  ### first given f_mk,k=1, calculate the other f_j,k according to the recursive formula
  for (k in 2:nnR_1){
    #print(k)              # for debug
    
    j <- nnR_1 - k + 1         # j = nn-1, ..., 1
    
    RHS_i <- c()
    RHS_i_xi <- c()
    RHS_i_beta <- c()
    
    ### Matrix operation to replace the for loop over the n subjects
    RHS_recursive_j <- rowSums(t(matrix(rep(exp((1+poly_fun(xiR, omega_1[j]))*Tbb[P_no,]), n), ncol=n)) * Gamma)
    RHS_recursive_j <- (omega_1[j]<=YR_1[DeltaR_1==0]) %*% ( RHS_recursive_j * exp(XXR_1[DeltaR_1==0,]%*%betaR))
    RHS_recursive_jp1 <- rowSums(t(matrix(rep(exp((1+poly_fun(xiR, omega_1[j+1]))*Tbb[P_no,]), n), ncol=n)) * Gamma)
    RHS_recursive_jp1 <- (omega_1[j+1]<=YR_1[DeltaR_1==0]) %*% (RHS_recursive_jp1 * exp(XXR_1[DeltaR_1==0,]%*%betaR))
    RHS_recursive <- RHS_recursive_j - RHS_recursive_jp1
    
    RHS_xi_j <- rowSums(t(matrix(rep(exp((1+poly_fun(xiR, omega_1[j]))*Tbb[P_no,]) * Tbb[P_no,], n), ncol=n)) * Gamma)
    RHS_xi_j <- (omega_1[j]<=YR_1[DeltaR_1==0]) %*% (RHS_xi_j * exp(XXR_1[DeltaR_1==0,]%*%betaR))
    RHS_xi_j <- RHS_xi_j * poly_bases(omega_1[j])
    RHS_xi_jp1 <- rowSums(t(matrix(rep(exp((1+poly_fun(xiR, omega_1[j+1]))*Tbb[P_no,]) * Tbb[P_no,], n), ncol=n)) * Gamma)
    RHS_xi_jp1 <- (omega_1[j+1]<=YR_1[DeltaR_1==0]) %*% (RHS_xi_jp1 * exp(XXR_1[DeltaR_1==0,]%*%betaR))
    RHS_xi_jp1 <- RHS_xi_jp1 * poly_bases(omega_1[j+1])
    RHS_xi <- RHS_xi_j - RHS_xi_jp1
    
    
    RHS_beta_j <- rowSums(t(matrix(rep(exp((1+poly_fun(xiR, omega_1[j]))*Tbb[P_no,]), n), ncol=n)) * Gamma)
    RHS_beta_j <- matrix(rep(RHS_beta_j * exp(XXR_1[DeltaR_1==0,]%*%betaR),nbetaR), ncol=nbetaR) * XXR_1[DeltaR_1==0,]
    RHS_beta_j <- (omega_1[j]<=YR_1[DeltaR_1==0]) %*% RHS_beta_j
    RHS_beta_jp1 <- rowSums(t(matrix(rep(exp((1+poly_fun(xiR, omega_1[j+1]))*Tbb[P_no,]), n), ncol=n)) * Gamma)
    RHS_beta_jp1 <- matrix(rep(RHS_beta_jp1 * exp(XXR_1[DeltaR_1==0,]%*%betaR),nbetaR), ncol=nbetaR) * XXR_1[DeltaR_1==0,]
    RHS_beta_jp1 <- (omega_1[j+1]<=YR_1[DeltaR_1==0]) %*% RHS_beta_jp1
    RHS_beta <- RHS_beta_j - RHS_beta_jp1
    # 
    ### f_j,k are found
    f[j] <- (1/f[j+1] + sum(RHS_recursive))^(-1)
    
    ### next, use the recursive formula in the calculation of the derivative for eqn3
    ### these temporary results will be used in the calculation of the derivatives for eqn1-2
    
    tempfalpha_xi <- -f[j+1]^(-2) * f1[j+1,1:nxiR] + RHS_xi
    tempfalpha_beta <- -f[j+1]^(-2) * f1[j+1,(nxiR+1):(nxiR+nbetaR)] + RHS_beta
    tempfalpha_alpha <- -f[j+1]^(-2) * f1[j+1,(nxiR+nbetaR+1)]
    
    ### part of partial f_mk,k/partial alpha and partial f_mk,k/partial beta are found
    f1[j,] = -f[j]^2*c(tempfalpha_xi, tempfalpha_beta, tempfalpha_alpha)
    
  }
  
  
  return(list(f=f, f1=f1))
}


newton <- function(oldxiR, oldbetaR, oldalphaR, n, IDR, YR, DeltaR, XXR,P_no=1, omega, NNR, nnR, mb, Tbb, Gamma,
                   indIDR, nxiR, designMat){
  ### To solve the equations, we calculate the values of the equations and 
  ### the 1st derivatives of the eqns
  
  sumLHS <- 0
  sumRHS <- 0
  
  nxiR <- length(oldxiR)
  nbetaR <- length(oldbetaR)
  
  
  ### NOTE: do_not assign 1 or 2 to P_no here
  fs <- recursive(oldxiR, oldbetaR, oldalphaR, n, IDR, YR, DeltaR, XXR,  P_no, omega, NNR, nnR, mb, Tbb, Gamma)
  f <- fs$f   # f_j,k
  f1 <- fs$f1 # partial derivatives
  
  
  #----------------- use matrix operation to speed up calculation -----------------------#
  EXPB <- exp(matrix(1+sapply(omega, function(x) poly_fun(oldxiR,x)))%*%matrix(Tbb[P_no,],nrow=1))  # nnR*no_nodes
  EXPX <- exp(XXR[DeltaR==0,]%*%oldbetaR)
  
  ### (1)
  LHS_xi <- t(matrix(sapply(YR, function(x) poly_bases(x)), nrow = nxiR))
  LHS_xi[DeltaR==0] <- 0
  LHS_xi =  indIDR %*% LHS_xi  
  Ebik <- Gamma %*% Tbb[P_no,]
  LHS_xi <- LHS_xi * matrix(rep(Ebik,nxiR), ncol=nxiR)
  LHS_xi <- matrix(colSums(LHS_xi))
  
  RHS_xi_all <- matrix(NA, nrow=nxiR)
  for(q_xi in 1:nxiR){
    RHS_xi <-  EXPB *
      t(matrix(rep(matrix(Tbb[P_no,]), nnR), ncol=nnR)) *
      matrix(rep(matrix(f), mb), ncol=mb) *
      matrix(rep(sapply(omega, function(x) poly_bases(x)[q_xi]), mb), ncol=mb)
    RHS_xi <- designMat %*% RHS_xi
    RHS_xi <- RHS_xi * exp(matrix( rep(XXR[DeltaR==0,]%*%oldbetaR,mb), ncol=mb))
    RHS_xi <- sum(rowSums(RHS_xi * Gamma))
    RHS_xi_all[q_xi] <- RHS_xi
    
  }
  
  eqn_xi <- LHS_xi - RHS_xi_all
  
  ### (2) 
  LHS_beta <- colSums(XXR[DeltaR==1,])
  RHS_beta <- EXPB * matrix(rep(matrix(f), mb), ncol=mb)
  RHS_beta <- designMat %*% RHS_beta   # each line corresponds to sum_j (those related to j)
  RHS_beta <- rowSums(RHS_beta * Gamma)    
  RHS_beta <- matrix(rep(matrix(RHS_beta),nbetaR),ncol=nbetaR) * XXR[DeltaR==0,] * matrix(rep(EXPX,nbetaR),ncol=nbetaR)
  RHS_beta <- colSums(RHS_beta)
  eqn_beta <- LHS_beta - RHS_beta
  
  ### (3) 
  # see the last section 
  
  ### 1.
  eqn_xi_xi <- matrix(NA, nxiR, nxiR)
  for (p_xi in 1:nxiR){
    for (q_xi in 1:nxiR){
      LHS_xi_xi <- EXPB *
        t(matrix(rep(matrix(Tbb[P_no,]), nnR), ncol=nnR))^2 *
        matrix(rep(sapply(omega, function(x) poly_bases(x)[p_xi]*poly_bases(x)[q_xi]), mb), ncol=mb) *
        matrix(rep(matrix(f), mb), ncol=mb) 
      LHS_xi_xi <- designMat %*% LHS_xi_xi
      LHS_xi_xi <-  rowSums(LHS_xi_xi * Gamma)
      LHS_xi_xi <- LHS_xi_xi %*% EXPX
      
      RHS_xi_xi <- EXPB *
        t(matrix(rep(matrix(Tbb[P_no,]), nnR), ncol=nnR)) *
        matrix(rep(sapply(omega, function(x) poly_bases(x)[p_xi]), mb), ncol=mb) *
        matrix(rep(matrix(f1[,q_xi]),mb), ncol=mb)                
      RHS_xi_xi <- designMat %*% RHS_xi_xi
      RHS_xi_xi <-  rowSums(RHS_xi_xi * Gamma)
      RHS_xi_xi <- RHS_xi_xi %*% EXPX
      
      eqn_xi_xi[p_xi,q_xi] <- - LHS_xi_xi - RHS_xi_xi
    }
  }
  
  ### 2.
  LHS_xi_beta_all <- matrix(NA, nxiR, nbetaR)
  RHS_xi_beta_all <- matrix(NA, nrow=nxiR, ncol=nbetaR)
  
  for (q_xi in 1:nxiR){
    LHS_xi_beta <- EXPB *
      t(matrix(rep(matrix(Tbb[P_no,]), nnR), ncol=nnR)) *
      matrix(rep(sapply(omega, function(x) poly_bases(x)[q_xi]), mb), ncol=mb) *
      matrix(rep(matrix(f), mb), ncol=mb) 
    LHS_xi_beta <- designMat %*% LHS_xi_beta
    LHS_xi_beta <- rowSums(LHS_xi_beta * Gamma)
    LHS_xi_beta <- matrix(rep(LHS_xi_beta*EXPX, nbetaR),ncol=nbetaR) * XXR[DeltaR==0,]
    LHS_xi_beta_all[q_xi,] <- colSums(LHS_xi_beta)
    
    for (q_beta in 1:nbetaR){
      RHS_xi_beta <- EXPB *
        t(matrix(rep(matrix(Tbb[P_no,]), nnR), ncol=nnR)) *
        matrix(rep(sapply(omega, function(x) poly_bases(x)[q_xi]), mb), ncol=mb) *
        matrix(rep(matrix(f1[, nxiR+q_beta]),mb), ncol=mb)               
      RHS_xi_beta <- designMat %*% RHS_xi_beta
      RHS_xi_beta <-  rowSums(RHS_xi_beta * Gamma)
      RHS_xi_beta <- RHS_xi_beta %*% EXPX
      RHS_xi_beta_all[q_xi,q_beta] <- RHS_xi_beta 
    }
  }
  eqn_xi_beta <- -LHS_xi_beta_all - RHS_xi_beta_all  
  
  ### 3.
  eqn_xi_alpha <- matrix(NA, nxiR, 1)
  for (q_xi in 1:nxiR){
    temp <- EXPB *
      t(matrix(rep(matrix(Tbb[P_no,]), nnR), ncol=nnR)) *
      matrix(rep(sapply(omega, function(x) poly_bases(x)[q_xi]), mb), ncol=mb) *
      matrix(rep(matrix(f1[, nxiR+nbetaR+1]),mb), ncol=mb)               
    temp <- designMat %*% temp
    temp <-  rowSums(temp * Gamma)
    eqn_xi_alpha[q_xi,] <- -temp %*% EXPX
  }
  
  
  ### 4.  
  LHS_beta_xi <- t(matrix(LHS_xi_beta_all, nrow=nxiR, ncol=nbetaR))
  RHS_beta_xi_all <- matrix(NA, nrow=nbetaR, ncol=nxiR)
  for (q_xi in 1:nxiR){
    RHS_beta_xi <- EXPB *
      matrix(rep(f1[,q_xi], mb), ncol=mb)
    RHS_beta_xi <- designMat %*% RHS_beta_xi 
    RHS_beta_xi <- rowSums(RHS_beta_xi * Gamma)
    RHS_beta_xi <- colSums(matrix(rep((RHS_beta_xi*EXPX), nbetaR), ncol=nbetaR) * XXR[DeltaR==0,])
    RHS_beta_xi_all[,q_xi] <- RHS_beta_xi
  }
  eqn_beta_xi <- - LHS_beta_xi - RHS_beta_xi_all
  
  ### 5.
  LHS_beta_beta <- rowSums((designMat %*% (EXPB *  matrix(rep(matrix(f), mb), ncol=mb))) * Gamma) 
  LHS_beta_beta <- t((matrix(rep(LHS_beta_beta * EXPX,nbetaR),ncol=nbetaR) * XXR[DeltaR==0,])) %*% (XXR[DeltaR==0,])
  
  RHS_beta_beta_all <- matrix(NA, nbetaR, nbetaR)
  for (q_beta in 1:nbetaR){
    RHS_beta_beta <- rowSums((designMat %*% (EXPB * matrix(rep(f1[,nxiR+q_beta],mb),ncol=mb))) * Gamma)
    RHS_beta_beta <- colSums(matrix(rep(RHS_beta_beta * EXPX, nbetaR),ncol=nbetaR) *  XXR[DeltaR==0,])
    RHS_beta_beta_all[,q_beta] <- RHS_beta_beta
  }  
  eqn_beta_beta <- - LHS_beta_beta - RHS_beta_beta_all
  
  ### 6.
  eqn_beta_alpha <- rowSums((designMat %*% (EXPB * matrix(rep(f1[,(nxiR+nbetaR+1)], mb), ncol=mb)) ) * Gamma) 
  eqn_beta_alpha <- - matrix(t(matrix(eqn_beta_alpha * EXPX)) %*% XXR[DeltaR==0,])
  
  
  eqnalpha <- sum(f)-1   # (3) is here
  
  
  ### 7-9 are here
  eqnalpha_xi <- matrix(colSums(matrix(f1[,(1:nxiR)],ncol=nxiR)), nrow=1, ncol=nxiR)
  
  eqnalpha_beta <- matrix(colSums(f1[,(nxiR+1):(nxiR+nbetaR)]), nrow=1, ncol=nbetaR)
  
  eqnalpha_alpha <- matrix(sum(f1[,nxiR+nbetaR+1]), nrow=1, ncol=1)
  
  
  return(list(eqnxi=eqn_xi, eqnbeta=eqn_beta, eqnalpha=eqnalpha,
              eqnxi_xi=eqn_xi_xi, eqnxi_beta=eqn_xi_beta, eqnxi_alpha=eqn_xi_alpha,
              eqnbeta_xi=eqn_beta_xi, eqnbeta_beta=eqn_beta_beta, eqnbeta_alpha=eqn_beta_alpha,
              eqnalpha_xi=eqnalpha_xi, eqnalpha_beta=eqnalpha_beta, eqnalpha_alpha=eqnalpha_alpha))
  
  
}




estimate_R_cov_simple <- function(newbetaR,newxiR,nxi,no_var,omega,ZZR,newXR,newLambdaR,newlambdaR,YR,DeltaR,nnR,
                                  IDR,OYR,Dbeta2R,DbetaLambdaR,DLambda2R,DxibetaR,DxiLambdaR,Dxi2R,
                                  nbetaonly,mb,newbb,newweight,i){
  
  TXR = matrix(newXR[IDR==i, ], nrow=sum(IDR==i))  # X matrix for subject i
  TLambdaR = newLambdaR[IDR==i]                    # the cumulative baseline intensity at the observation times of subject i
  
  
  
  TYR=YR[IDR==i,] 
  TDeltaR=DeltaR[IDR==i,]
  TZZR=ZZR[IDR==i,]
  
  
  ### compute the summation over j=1,..,m_k, for use in the next for loop.
  
  sumj_ind_EXZ_Lam <- matrix(0, 1, mb)
  sumj_ind_EXZ_Lam_Omegat <- matrix(0, nxi, mb) 
  sumj_ind_EXZ_Lam_Omegat_Omegat <- matrix(0, nxi*nxi, mb)
  
  ### Also compute first derivative wrt LambdaR, j=1,.., m_k. Each row in the matrix correspond to the jth entry.
  tempDLambdaR=matrix(0, nnR, mb) # first derivative wrt LambdaR, at each node
  
  
  for (j in 1:nnR){
    
    omega_j <- omega[j]
    ind <- 1*(omega_j <= TYR[TDeltaR==0])
    
    if (ind==1){
      
      #-------------------#
      EXZ <- exp( (1+poly_fun(newxiR,omega_j))*(ZZR[i,]%*%newbb)  + matrix(rep(TXR[1,]%*%newbetaR[-c(1,2+no_var)],mb),ncol=mb) )
      Lam <- newlambdaR[j]
      Omegat <- poly_bases(omega_j)    # Omegat is a column vector
      
      #-------------------#
      sumj_ind_EXZ_Lam <- sumj_ind_EXZ_Lam + ind*EXZ*Lam
      
      sumj_ind_EXZ_Lam_Omegat <- sumj_ind_EXZ_Lam_Omegat + Omegat %*% (ind*EXZ*Lam)  # Note: use (ind*EXZ*Lam), not sumj_ind_EXZ_Lam
      
      sumj_ind_EXZ_Lam_Omegat_Omegat <- sumj_ind_EXZ_Lam_Omegat_Omegat + matrix(Omegat %*% t(Omegat)) %*% (ind*EXZ*Lam)  # (n_xiR * n_xiR) * n_nodes
      
      #-------------------#
      
      ### only added the second part. The first part will be added outside this loop over j
      tempDLambdaR[j,] <- tempDLambdaR[j,] - ind * EXZ
      
      
      ### partial beta partial Lambda; omega_j, j=1,..,nnR, represented by column
      DbetaLambdaR[,j] = DbetaLambdaR[,j] - ind * as.numeric(EXZ %*% as.matrix(newweight[i,])) * as.matrix(TXR[1,])
      
      ### the resulting Dbeta2R should be the same as the one calculated outside the loop.-> Correct.
      # Dbeta2R = Dbeta2R - ind*as.numeric(EXZ%*%as.matrix(newweight[i,]))*Lam*(as.matrix(TXR[1,])%*%t(as.matrix(TXR[1,])))
      
      
      ### the resulting DxibetaR should be the same as below -> Correct.
      # DxibetaR = DxibetaR - ind * as.numeric((EXZ*as.numeric(ZZR[i,]%*%newbb))%*% newweight[i,]) * Lam *
      #                                           poly_bases(omega_j) %*% t(as.matrix(TXR[1,]))   # check poly_bases
      
      # partial xi partial Lambda
      DxiLambdaR[,j] =  DxiLambdaR[,j] - ind * as.numeric(   (EXZ*as.numeric(ZZR[i,]%*%newbb))   %*% as.matrix(newweight[i,])   ) * poly_bases(omega_j)
      
      
    }
    
    
  }
  
  
  
  ### first derivative wrt beta, at each node
  tempDbetaR=sum(TDeltaR==1)*matrix(rep(as.matrix(TXR[1,]),mb), ncol=mb) - as.matrix(TXR[1,])%*%sumj_ind_EXZ_Lam
  # So, this is info contributed by subject i, to partial logL/partial beta
  
  
  ### second derivative wrt beta; check if this is the same as computation in the loop;
  Dbeta2R = Dbeta2R - as.numeric(sumj_ind_EXZ_Lam%*%as.matrix(newweight[i,]))*(as.matrix(TXR[1,])%*%t(as.matrix(TXR[1,])))
  
  
  #----------#
  # Next, add 1/Lambda{} to first derivatives wrt Lambda{}, tempDLambdaR; 
  # and compute second derivative wrt Lambda{.}, DLambda2R
  # ith subject's contribution, is when the event time is this subject's event time
  
  if(length(which(omega%in%TYR[TDeltaR==1])) != sum(TDeltaR)) {
    stop("omega and TYR have different number decimal points")
  } 
  
  ### first derivative wrt Lambda, at each node
  tempDLambdaR[which(omega%in%TYR[TDeltaR==1]),] <- matrix(rep(as.matrix(1/newlambdaR[which(omega%in%TYR[TDeltaR==1])]),mb),ncol=mb) - 
    tempDLambdaR[which(omega%in%TYR[TDeltaR==1]),] 
  ### second derivative wrt Lambda
  if (sum(TDeltaR)==1) { 
    
    DLambda2R[which(omega%in%TYR[TDeltaR==1]), which(omega%in%TYR[TDeltaR==1])] <- DLambda2R[which(omega%in%TYR[TDeltaR==1]), which(omega%in%TYR[TDeltaR==1])]-
      newlambdaR[which(omega%in%TYR[TDeltaR==1])]^(-2)
    
  } else if (sum(TDeltaR)>1){
    
    diag(DLambda2R[which(omega%in%TYR[TDeltaR==1]), which(omega%in%TYR[TDeltaR==1])]) <- diag(DLambda2R[which(omega%in%TYR[TDeltaR==1]), which(omega%in%TYR[TDeltaR==1])])-
      newlambdaR[which(omega%in%TYR[TDeltaR==1])]^(-2)
    
  }
  
  
  
  ### this should return the same value as the one in the loop; check- > Correct.
  DxibetaR = DxibetaR - (   (sumj_ind_EXZ_Lam_Omegat * as.numeric(ZZR[i,]%*%newbb)) %*% as.matrix(newweight[i,])   ) %*% t(as.matrix(TXR[1,]))
  
  # second derivative wrt xi
  Dxi2R = Dxi2R - matrix((  sumj_ind_EXZ_Lam_Omegat_Omegat * t(matrix(rep((ZZR[i,]%*%newbb)^2, nxi*nxi), ncol=nxi*nxi))  )  %*% as.matrix(newweight[i,]),
                         ncol=nxi)
  
  
  ### first derivative wrt xi, before multiplying the weights; Note: vectorized event time in the first term.
  sum_tikl_tilde <- matrix(0, nrow=nxi)
  if (sum(TDeltaR)>0){
    for (lll in 1:sum(TDeltaR)){
      sum_tikl_tilde <- sum_tikl_tilde + poly_bases(TYR[TDeltaR==1][lll])
    }
  }
  
  tempDxiR= matrix(rep(sum_tikl_tilde,mb), ncol=mb) * t(matrix(rep((ZZR[i,]%*%newbb), nxi), ncol=nxi)) - sumj_ind_EXZ_Lam_Omegat * t(matrix(rep((ZZR[i,]%*%newbb),nxi), ncol=nxi))
  
  
  list(tempDbetaR=tempDbetaR,     # n_beta * n_nodes
       tempDLambdaR=tempDLambdaR, # n_events * n_nodes
       tempDxiR=tempDxiR,         # n_xi * n_nodes
       Dbeta2R=Dbeta2R,           # n_beta * n_beta
       DbetaLambdaR=DbetaLambdaR, # n_beta * n_Lambda
       DLambda2R=DLambda2R,        # n_events * n_events
       DxiLambdaR=DxiLambdaR,
       DxibetaR=DxibetaR,
       Dxi2R=Dxi2R
  )     
  
}


cov_Sigma = function(Sigma,Tbb,weight,i){
  
  mb=ncol(Tbb)

  ISigma=inv(Sigma)
  
  
  nSig=nrow(Sigma)   
  
  scoreSig=matrix(0, nSig*(nSig+1)/2, mb);
  dscoreSig=matrix(0, nSig*(nSig+1)/2, nSig*(nSig+1)/2)
  
  
  ind1=0
  for (k in 1:nSig){
    
    for (m in k:nSig){
      ind1=ind1+1
      dkm=matrix(0,nSig,nSig)
      dkm[k,m]=1
      dkm[m,k]=1
      scoreSig[ind1,]=t(diag(t(Tbb)%*%(ISigma%*%dkm%*%ISigma)%*%Tbb))/2-sum(sum(ISigma*dkm))/2  
      ind2=0                                                                                    
      for (p in 1:nSig){
        for (q in p:nSig){
          ind2=ind2+1;
          dpq=matrix(0, nSig, nSig)
          dpq[p,q]=1
          dpq[q,p]=1
          temp=ISigma%*%dpq%*%ISigma%*%dkm%*%ISigma+ISigma%*%dkm%*%ISigma%*%dpq%*%ISigma
          dscoreSig[ind1,ind2]=-t(diag(t(Tbb)%*%temp%*%Tbb)/2)%*%matrix(weight[i,]) + sum(sum((ISigma%*%dpq%*%ISigma)*dkm))/2
          dscoreSig[ind2,ind1]=dscoreSig[ind1,ind2]
        }
      }                                                                    
      
    }
  }
  
  
  
  scoreSig_i=scoreSig
  dscoreSig_i=dscoreSig
  
  
  return(list(scoreSig_i=scoreSig_i,
              dscoreSig_i=dscoreSig_i))
  
}



