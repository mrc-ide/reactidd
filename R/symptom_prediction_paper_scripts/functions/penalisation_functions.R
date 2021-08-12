GetArgmaxId=function(calib_object=NULL, S=NULL){
  if ((is.null(calib_object))&(is.null(S))){
    stop("Please provide calib_object or S.")
  }
  
  if (is.null(S)){
    argmax_id=matrix(NA, nrow=ncol(calib_object$Lambda), ncol=2)
    for (block_id in 1:ncol(calib_object$Lambda)){
      if (ncol(calib_object$Lambda)==1){
        myS=calib_object$S
      } else {
        myS=calib_object$S_blocks[,block_id,drop=FALSE]
      }
      myS[is.na(myS)]=0
      myid=which.max(myS[,1])
      argmax_id[block_id,]=c(myid, which(calib_object$params$pi_list==calib_object$P[myid,block_id]))
    }
  } else {
    argmax_id=matrix(NA, nrow=1, ncol=2)
    myS=apply(S,1,max,na.rm=TRUE)
    myS[is.na(myS)]=0
    myid=which.max(myS)
    argmax_id[1,]=c(myid, max(which(S[myid,]==myS[myid])))
  }
  colnames(argmax_id)=c("lambda_id","pi_id")
  return(argmax_id)
}


GetArgmax=function(calib_object){
  argmax=matrix(NA, nrow=ncol(calib_object$Lambda), ncol=2)
  for (block_id in 1:ncol(calib_object$Lambda)){
    if (ncol(calib_object$Lambda)==1){
      myS=calib_object$S
    } else {
      myS=calib_object$S_blocks[,block_id,drop=FALSE]
    }
    myS[is.na(myS)]=0
    myid=which.max(myS[,1])
    argmax[block_id,]=c(calib_object$Lambda[myid,block_id], calib_object$P[myid,block_id])
  }
  colnames(argmax)=c("lambda","pi")
  return(argmax)
}


GetLambdaPath=function(lmax, lmin, cardinal=100){
  return(seq(sqrt(lmax),sqrt(lmin),length.out=cardinal)^2)
}


GetSubsample=function(ydata, family="gaussian", tau=0.5, resampling_method="subsampling", resampling_function=NULL){
  # resampling_function has to be a function of tau
  # argument resampling_method is ignored if resampling_function is provided
  if (resampling_method=="subsampling"){
    if (family=="gaussian"){
      if (is.null(resampling_function)){
        s=sample(length(ydata), size=tau*length(ydata))
      } else {
        s=resampling_function(tau=tau)
      }
    }
    if (family=="binomial"){
      if (is.null(resampling_function)){
        s0=sample(which(ydata=="0"), size=tau*sum(ydata=="0"))
        s1=sample(which(ydata=="1"), size=tau*sum(ydata=="1"))
        s=c(s0,s1)
      } else {
        s=resampling_function(tau=tau)
      }
    }
    if (family=="cox"){
      if (is.null(resampling_function)){
        s0=sample(which(ydata[,2]=="0"), size=tau*sum(ydata[,2]=="0"))
        s1=sample(which(ydata[,2]=="1"), size=tau*sum(ydata[,2]=="1"))
        s=c(s0,s1)
      } else {
        s=resampling_function(tau=tau)
      }
    }
  }
  if (resampling_method=="bootstrap"){
    if (family=="gaussian"){
      if (is.null(resampling_function)){
        s=sample(length(ydata), size=length(ydata), replace=TRUE)
      } else {
        s=resampling_function(tau=tau)
      }
    }
    if (family=="binomial"){
      if (is.null(resampling_function)){
        s0=sample(which(ydata=="0"), size=sum(ydata=="0"), replace=TRUE)
        s1=sample(which(ydata=="1"), size=sum(ydata=="1"), replace=TRUE)
        s=c(s0,s1)
      } else {
        s=resampling_function(tau=tau)
      }
    }
    if (family=="cox"){
      if (is.null(resampling_function)){
        s0=sample(which(ydata[,2]=="0"), size=sum(ydata[,2]=="0"), replace=TRUE)
        s1=sample(which(ydata[,2]=="1"), size=sum(ydata[,2]=="1"), replace=TRUE)
        s=c(s0,s1)
      } else {
        s=resampling_function(tau=tau)
      }
    }
  }
  return(s)
}


ComputePFER=function(q, pi, N, K, method="MB"){
  if (!method%in%c("MB", "SS")){
    stop('PFER method must be "MB" or "SS".')
  }
  if (method=="MB"){
    upperbound=1/(2*pi-1)*q^2/N
  } 
  if (method=="SS"){
    ### Adapted code from stabsel package (to avoid "out of bounds" error)
    cutoff=pi
    B=K
    theta=q/N
    if (cutoff <= 3/4) {
      tmp=2 * (2 * cutoff - 1 - 1/(2*B))
    } else {
      tmp=(1 + 1/B)/(4 * (1 - cutoff + 1 / (2*B)))
    }
    upperbound=q^2/N/tmp
    if ((cutoff < 1/2 + min(theta^2, 1 / (2*B) + 3/4 * theta^2))|(cutoff>1)){
      # If out of bounds for SS formula under unimodality assumption, set to Inf
      upperbound=Inf
    }
  }
  return(upperbound)
}


ComputeFDP=function(PFER, selection_proportions, pi){
  # Computing the proportion of false discoveries among discoveries
  if (length(dim(selection_proportions))==2){
    selprops=selection_proportions[upper.tri(selection_proportions)] # vector of selection proportions
  } else {
    selprops=selection_proportions
  }
  S=sum(selprops>=pi) # number of stability-selected edges
  if (S!=0){
    FDP=PFER/S
  } else {
    FDP=0
  }
  return(FDP)
}


GetBinomialProbabilities=function(q, N, pi, K){
  p_1=pbinom(round(K*pi), size=K, prob=q/N, log.p=TRUE) # proportion <= (1-pi)
  p_2=log(pbinom(round(K*(1-pi))-1, size=K, prob=q/N)-pbinom(round(K*pi), size=K, prob=q/N)) # 1-pi < proportion < pi
  if (is.infinite(p_2)|is.na(p_2)){
    p_2=0
    for (i in seq(round(K*(1-pi))-1, round(K*pi)+1)){
      p_2=p_2+dbinom(i, size=K, prob=q/N)
    }
    p_2=log(p_2)
  }
  p_3=pbinom(round(K*(1-pi))-1, size=K, prob=q/N, lower.tail=FALSE, log.p=TRUE) # proportion >= pi 
  
  if (abs(exp(p_1)+exp(p_2)+exp(p_3)-1)>1e-3){
    print(paste("N:", N))
    print(paste("q:", q))
    print(paste("K:", K))
    print(paste("pi:", pi))
    stop(paste0("Probabilities do not sum to 1 (Binomial distribution) \n p_1+p_2+p_3=", exp(p_1)+exp(p_2)+exp(p_3)))
  }
  return(list(p_1=p_1, p_2=p_2, p_3=p_3))
}


GetBinomialScore=function(stab_iter, q=NULL, N=NULL, pi, K){
  if (is.matrix(stab_iter)){
    stab_iter=stab_iter[upper.tri(stab_iter)]
  }
  N=length(stab_iter)
  
  if (is.null(q)){
    q=round(sum(stab_iter))
  }
  
  p_vect=GetBinomialProbabilities(q, N, 1-pi, K)
  
  S_0=sum(stab_iter<=(1-pi)) # stable-out
  S_1=sum(stab_iter>=pi) # stable-in
  U=sum((stab_iter<pi)&(stab_iter>(1-pi))) # unstable
  
  if (S_0+S_1+U!=N){
    stop(paste0("Inconsistency in number of edges \n S_0+S_1+U=", S_0+S_1+U, " instead of ", N))
  }
  
  l=S_0*p_vect$p_1+U*p_vect$p_2+S_1*p_vect$p_3
  
  if (is.infinite(l)){
    l=NA
  }
  
  return(l)
}


StabilityCalibRegression=function(xdata, ydata, pk=NULL, Lambda, pi_list=seq(0.6,0.9,by=0.01), K=100, tau=0.5, seed=1, 
                                  family="gaussian", penalty.factor=NULL,
                                  resampling_method="subsampling", resampling_function=NULL, PFER_method="MB", PFER_thr=+Inf, FDP_thr=+Inf, 
                                  verbose=TRUE, debug=FALSE, ...){
  # Browser mode for debugging
  if (debug){
    browser()
  }
  
  # Prepare dimensions
  if (is.null(pk)){
    pk=ncol(xdata)
  }
  
  # Prepare penalty.factor
  if (is.null(penalty.factor)){
    penalty.factor=rep(1, sum(pk))
  }
  
  # Checking consistency between arguments 
  if (!is.null(pk)){
    if (sum(pk)!=ncol(xdata)){
      stop("Argument pk is not consistent with the number of variables in X Please make sure that sum(pk) is equal to ncol(X).")
    }
  }
  
  # Refine K if using complementary pairs (SS)
  if (PFER_method=="SS"){
    K=ceiling(K/2)*2
  }
  
  # Name columns of xdata
  if (is.null(colnames(xdata))){
    colnames(xdata)=paste0("var", 1:ncol(xdata))
  }
  
  # Create matrix with block indices
  bigblocks_vect=diag(GetBlockMatrix(pk))
  N_blocks=unname(table(bigblocks_vect))
  blocks=unique(bigblocks_vect)
  names(N_blocks)=blocks
  nblocks=max(blocks)
  
  # Re-format Lambda
  if (is.vector(Lambda)){
    Lambda=matrix(rep(Lambda, nblocks), ncol=nblocks)
  }
  rownames(Lambda)=paste0("s",seq(0,nrow(Lambda)-1))
  
  # Re-format ydata
  if (is.vector(ydata)|is.factor(ydata)){
    ydata=matrix(ydata, ncol=1)
  }
  
  # Check consistency between X and Y
  if (nrow(xdata)!=nrow(ydata)){
    stop("Different numbers of observations in xdata and ydata.")
  }
  
  # Prepare the PFER and FDP thresholds
  if (length(PFER_thr)==1){
    PFER_thr_blocks=ceiling(prop.table(N_blocks)*PFER_thr)
  } else {
    if (length(PFER_thr)==nblocks){
      PFER_thr_blocks=PFER_thr
    } else {
      stop(paste0("Please provide a single number or as many values as there are blocks in PFER_thr (currently ",length(PFER_thr)," values provided)."))
    }
  }
  if (length(FDP_thr)==1){
    FDP_thr_blocks=rep(FDP_thr, nblocks)
  } else {
    if (length(FDP_thr)==nblocks){
      FDP_thr_blocks=FDP_thr
    } else {
      stop(paste0("Please provide a single number or as many values as there are blocks in FDP_thr (currently ",length(FDP_thr)," values provided)."))
    }
  }
  
  # Initialise objects to be filled
  Q=P=Q_s=matrix(NA, nrow=nrow(Lambda), ncol=1) # average number of included variable per lambda
  loglik=PFER=FDP=matrix(NA, nrow=nrow(Lambda), ncol=length(pi_list))
  N=N_block=ncol(xdata)
  Beta=array(0, dim=c(nrow(Lambda), ncol(xdata), K))
  rownames(Beta)=rownames(Lambda)
  colnames(Beta)=colnames(xdata)
  best_loglik=best_PFER=best_FDP=matrix(NA, nrow=nrow(Lambda), ncol=nblocks)
  
  # Printing message 
  if (verbose){
    if (all(!is.infinite(PFER_thr_blocks))){
      print("Threshold(s) in PFER:")
      print(PFER_thr_blocks)
    }
    if (all(!is.infinite(FDP_thr_blocks))){
      print("Threshold(s) in FDP:")
      print(FDP_thr_blocks)
    }
  }
  
  # Computation of the selection proportions over Lambda
  pb=txtProgressBar(style=3)
  if (PFER_method=="MB"){
    for (k in 1:K){
      # print(k)
      setTxtProgressBar(pb, k/K)
      set.seed(k)
      s=GetSubsample(ydata=ydata, family=family, tau=tau, resampling_method=resampling_method, resampling_function=resampling_function)
      Xsub = xdata[s,]
      Ysub = ydata[s,]
      model_sub = glmnet(x=Xsub, y=Ysub, lambda=Lambda[,1], family=family, ...)
      while (is.infinite(model_sub$lambda[1])){
        print("resampling")
        s=GetSubsample(ydata=ydata, family=family, tau=tau, resampling_method=resampling_method, resampling_function=resampling_function)
        Xsub = xdata[s,]
        Ysub = ydata[s,]
        model_sub = glmnet(x=Xsub, y=Ysub, lambda=Lambda[,1], family=family, ...)
      }
      beta_sub=t(as.matrix(coef(model_sub)))
      beta_sub=beta_sub[,colnames(xdata)] # removing the intercept if included
      Beta[rownames(beta_sub),colnames(beta_sub),k]=beta_sub
    }
  }
  if (PFER_method=="SS"){
    for (k in 1:ceiling(K/2)){
      setTxtProgressBar(pb, k/K)
      set.seed(k)
      s=GetSubsample(ydata=ydata, family=family, tau=tau, resampling_method=resampling_method, resampling_function=resampling_function)
      
      # Fist subset
      Xsub = xdata[s,]
      Ysub = ydata[s,]
      model_sub = glmnet(x=Xsub, y=Ysub, lambda=Lambda[,1], family=family, penalty.factor=penalty.factor)
      beta_sub=t(as.matrix(coef(model_sub)))
      beta_sub=beta_sub[,colnames(xdata)] # removing the intercept if included
      Beta[rownames(beta_sub),colnames(beta_sub),k]=beta_sub
      
      # Complementary subset
      Xsub = xdata[seq(1,nrow(xdata))[!seq(1,nrow(xdata))%in%s],]
      Ysub = ydata[seq(1,nrow(xdata))[!seq(1,nrow(xdata))%in%s],]
      model_sub = glmnet(x=Xsub, y=Ysub, lambda=Lambda[,1], family=family, penalty.factor=penalty.factor)
      beta_sub=t(as.matrix(coef(model_sub)))
      beta_sub=beta_sub[,colnames(xdata)] # removing the intercept if included
      Beta[rownames(beta_sub),colnames(beta_sub),ceiling(K/2)+k]=beta_sub
    }
  }
  
  cat("\n")
  if (K>1){
    bigstab=apply(Beta,c(1,2),FUN=function(x){sum(x!=0)})/K # selection proportions
  }
  
  # Computation of the stability score over Lambda and pi_list
  if (K>1){
    for (k in 1:nrow(Lambda)){
      q_block=mean(apply(Beta[k,,],2,FUN=function(x){sum(x!=0)}))
      Q[k,1]=round(q_block)
      for (j in 1:length(pi_list)){
        pi=pi_list[j]
        PFER[k,j]=ComputePFER(q=q_block, pi=pi, N=N_block, K=K, method=PFER_method)
        FDP[k,j]=ComputeFDP(PFER=PFER[k,j], pi=pi, selection_proportions=bigstab[k,])
        if ((PFER[k,j]<=PFER_thr)&(FDP[k,j]<=FDP_thr)){
          loglik[k,j]=GetBinomialScore(stab_iter=bigstab[k,], pi=pi, K=K)
        }
      }
    }
    rownames(loglik)=rownames(PFER)=rownames(FDP)=rownames(bigstab)
    colnames(loglik)=colnames(PFER)=colnames(FDP)=pi_list
    
    # Add constraint
    if ((!is.infinite(PFER_thr))|(!is.infinite(FDP_thr))){
      if (!is.infinite(PFER_thr)){
        loglik=ifelse(PFER>PFER_thr, yes=NA, no=loglik)
      } else {
        loglik=ifelse(FDP>FDP_thr, yes=NA, no=loglik)
      }
    }
    
    # Computation of the best score by lambda
    block_id=1
    for (k in 1:nrow(Lambda)){
      tmp_loglik=loglik[k,]
      if (any(!is.na(tmp_loglik))){
        tmp_loglik[is.na(tmp_loglik)]=0
        myid=which.min(tmp_loglik)
        tmp_loglik[which(tmp_loglik==0)]=NA
        best_loglik[k,block_id]=tmp_loglik[myid]
        P[k,block_id]=pi_list[myid]
        Q_s[k,block_id]=sum(bigstab[k,]>=pi_list[myid])
        best_PFER[k,block_id]=PFER[k,myid]
        best_FDP[k,block_id]=FDP[k,myid]
      }
    }
  } else {
    for (k in 1:nrow(Lambda)){
      q_block=sum(myscreen$Beta[k,,1]!=0)
      Q[k,1]=round(q_block)
      for (j in 1:length(pi_list)){
        pi=pi_list[j]
        PFER[k,j]=ComputePFER(q=q_block,pi=pi,N=N_block)
      }
    }
  }
  
  # Prepare outputs
  if (nblocks==1){ # Only possible with 1 block so far
    if (K>1){
      return(list(S=-best_loglik, Lambda=Lambda, 
                  Q=Q, Q_s=Q_s, P=P,
                  PFER=best_PFER, FDP=best_FDP,
                  S_2d=-loglik, PFER_2d=PFER, FDP_2d=FDP, 
                  selprop=bigstab, Beta=Beta,
                  methods=list(family=family, resampling_method=resampling_method, PFER_method=PFER_method),
                  params=list(K=K, pi_list=pi_list, tau=tau, pk=pk, PFER_thr=PFER_thr, FDP_thr=FDP_thr, seed=seed, xdata=xdata, ydata=ydata)))
    } else {
      return(list(Q=Q, Beta=Beta, PFER_2d=PFER))
    }
  }
}


LambdaGridRegression=function(xdata, ydata, tau=0.5, seed=1, 
                              family="gaussian", penalty.factor=NULL,
                              resampling_method="subsampling", resampling_function=NULL, PFER_method="MB", PFER_thr=+Inf, FDP_thr=+Inf, 
                              lambda_path_factor=0.0001, max_density=0.3, 
                              lambda_path_refined_cardinal=100, verbose=TRUE){
  # Prepare dimensions
  pk=ncol(xdata)
  N=pk
  
  # Prepare data format
  if (is.vector(ydata)|is.factor(ydata)){
    ydata=matrix(ydata, ncol=1)
  }
  
  # Check consistency between X and Y
  if (nrow(xdata)!=nrow(ydata)){
    stop("Different numbers of observations in xdata and ydata.")
  }
  
  # Prepare penalty.factor
  if (is.null(penalty.factor)){
    penalty.factor=rep(1, sum(pk))
  }
  
  # Name columns of xdata
  if (is.null(colnames(xdata))){
    colnames(xdata)=paste0("var", 1:ncol(xdata))
  }
  
  # Get upperbound of Lambda
  set.seed(1)
  s=GetSubsample(ydata=ydata, family=family, tau=tau, resampling_method=resampling_method, resampling_function=resampling_function)
  set.seed(1)
  mycv=cv.glmnet(x=xdata[s,], y=ydata[s,], family=family, penalty.factor=penalty.factor)
  if (verbose){
    plot(mycv, las=1)
    print(paste("Minimum number of selected variables:", min(mycv$nzero)))
    print(paste("Maximum number of selected variables:", max(mycv$nzero)))
  }
  Lambda=cbind(GetLambdaPath(lmax=max(mycv$lambda), lmin=min(mycv$lambda), cardinal=lambda_path_refined_cardinal))
  
  return(Lambda)
}


CalibrateRegression=function(xdata, ydata, Lambda=NULL, pi_list=seq(0.6,0.9,by=0.01), K=100, tau=0.5, seed=1, 
                             family="gaussian", penalty.factor=NULL,
                             resampling_method="subsampling", resampling_function=NULL, PFER_method="MB", PFER_thr=+Inf, FDP_thr=+Inf, 
                             lambda_path_refined_cardinal=100,
                             verbose=TRUE, debug=FALSE, ...){
  if (is.null(Lambda)){
    # Define grid of lambda values (using glmnet implementation)
    Lambda=LambdaGridRegression(xdata=xdata, ydata=ydata, tau=tau, seed=seed, 
                                family=family, penalty.factor=penalty.factor,
                                resampling_method=resampling_method, resampling_function=resampling_function, PFER_method=PFER_method, PFER_thr=PFER_thr, FDP_thr=FDP_thr, 
                                lambda_path_factor=0.0001, max_density=0.3, 
                                lambda_path_refined_cardinal=lambda_path_refined_cardinal, verbose=verbose)
  }
  
  # Perform stability-based calibration
  out=StabilityCalibRegression(xdata=xdata, ydata=ydata, pk=NULL, Lambda=Lambda, pi_list=pi_list, K=K, tau=tau, seed=seed, 
                               family=family, penalty.factor=penalty.factor,
                               resampling_method=resampling_method, resampling_function=resampling_function, PFER_method=PFER_method, PFER_thr=PFER_thr, FDP_thr=FDP_thr, 
                               verbose=verbose, debug=debug)
  
  return(out)
}

