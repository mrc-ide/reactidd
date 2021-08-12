### Univariate function


univariateModels <- function(data, adjustments = NULL, modelname=NULL){
  
  results.df<- data.frame(var=react1_sympnames[c(1:26)], coef=NA,lower = NA,upper=NA,pval=NA,modname = modelname)
  results.adj.df<- data.frame(var=react1_sympnames[c(1:26)], coef=NA,lower = NA,upper=NA,pval=NA, modname = modelname)
  results.adj.interact.df<- data.frame(var=react1_sympnames[c(1:26)], coef=NA,lower = NA,pval=NA,upper=NA, modname = modelname)
  
  for (i in 1:26){
    print(i)
    var <- covid_yesnos[[i]]
    f <- as.formula(paste("estbinres ~ ", var))
    mod.unadj <- glm(formula = f, family = "binomial", data = data)
    
    mod.unadj.res <- jtools::summ(mod.unadj, exp = T)
    results.df[i,2:5] <- mod.unadj.res$coeftable[2,c(1,2,3,5)]
    
    if(!is.null(adjustments)){
      f2 <- as.formula(paste("estbinres ~ ", var, "+", paste(adjustments, collapse = "+")))
      mod.adj <- glm(formula = f2, family = "binomial", data = data)
      f3 <- as.formula(paste("estbinres ~ (", var, "+", paste(adjustments, collapse = "+"),")^2"))
      mod.adj.inter <- glm(formula = f3, family = "binomial", data = data)
      mod.adj.res <- jtools::summ(mod.adj, exp = T)
      mod.adj.inter.res <- jtools::summ(mod.adj.inter, exp = T)
      results.adj.df[i,2:5] <- mod.adj.res$coeftable[2,c(1,2,3,5)]
      results.adj.interact.df[i,2:5] <- mod.adj.inter.res$coeftable[2,c(1,2,3,5)]
    }
    
    
  }
  if(!is.null(adjustments)){
    return(list(crude=results.df,
                adj=results.adj.df,
                adj_inter=results.adj.interact.df
    ))
  }else{
    return(results.df)
  }
  
}


### Univariate function with interactions

# data = mod.dat
# adjustments = c("age", "gender") 
# modelname = "all"
# interactvar = "round_split"
univariateModelsInteract <- function(data, modelname=NULL, interactvar){
  
  results.df<- data.frame(var=react1_sympnames[c(1:26)], coef=NA,lower = NA,upper=NA,modname = modelname)
  results.interact.df<- data.frame(var=react1_sympnames[c(1:26)], coef=NA,lower = NA,upper=NA, pval=NA,
                                       roundsplit_coef=NA,roundsplit_lower = NA,roundsplit_upper=NA, roundsplit_pval=NA,
                                       interact_coef=NA,interact_lower = NA,interact_upper=NA, interact_pval=NA,
                                       modname = modelname)

  for (i in 1:26){
    print(i)
    var <- covid_yesnos[[i]]
    f <- as.formula(paste("estbinres ~ ", var, "*", interactvar))
    mod.unadj <- glm(formula = f, family = "binomial", data = data)
    mod.unadj.res <- jtools::summ(mod.unadj, exp = T)
    results.interact.df[i,2:13] <- c(mod.unadj.res$coeftable[2,c(1,2,3,5)],
                           mod.unadj.res$coeftable[3,c(1,2,3,5)],
                           mod.unadj.res$coeftable[4,c(1,2,3,5)])
                           
  }
    return(results.interact.df)
  
}





### Pairwise univariate
univariatePairwiseModels <- function(data){
  results.df<- data.frame(matrix(nrow = 26,ncol = 26,dimnames=list(react1_sympnames[c(1:26)],react1_sympnames[c(1:26)])))
  results.adj.df<- data.frame(matrix(nrow = 26,ncol = 26,dimnames=list(react1_sympnames[c(1:26)],react1_sympnames[c(1:26)])))

  for (i in 1:26){
    t0 <- Sys.time()
    for(j in 1:26){
    var1 <- covid_yesnos[[i]]
    var2 <- covid_yesnos[[j]]
    data[,"newvar"] <- data[,var1] * data[,var2]
    f <- as.formula(paste("estbinres ~ newvar"))
    mod.unadj <- glm(formula = f, family = "binomial", data = data)
    mod.unadj.res <- jtools::summ(mod.unadj, exp = T)
    results.df[i,j] <- mod.unadj.res$coeftable[2,1]
    }
    print(paste(round(100*i/26,0),"% complete"))
    t1 <- Sys.time()
    looptime <- t1-t0
    print(paste("Estimated time remaining:", round((26-i) * looptime, 2)))
    
  }
    return(results.df)
}




### Pairwise prevalence
univariatePairwisePrevs <- function(data){
  results.df<- data.frame(matrix(nrow = 26,ncol = 26,dimnames=list(react1_sympnames[c(1:26)],react1_sympnames[c(1:26)])))
  for (i in 1:26){
    t0 <- Sys.time()
    for(j in 1:26){
      var1 <- covid_yesnos[[i]]
      var2 <- covid_yesnos[[j]]
      data[,"newvar"] <- data[,var1] * data[,var2]
      results.df[i,j] <- sum(data[,"newvar"]) / nrow(data)
    }
    print(paste(round(100*i/26,0),"% complete"))
    t1 <- Sys.time()
    looptime <- t1-t0
    print(paste("Estimated time remaining:", round((26-i) * looptime, 2)))
    
  }
  return(results.df)
}



# Jaccard function --------------------------------------------------------

jaccardCalc <- function(var1, var2){
  numer <- sum(var1+var2 == 2)
  denom <- length(var1)
  return(numer/denom)
}




### Pairwise correlation
univariatePairwiseSimilarity <- function(data, similarity = "jaccard"){
  results.df<- data.frame(matrix(nrow = 26,ncol = 26,dimnames=list(react1_sympnames[c(1:26)],react1_sympnames[c(1:26)])))
  for (i in 1:26){
    t0 <- Sys.time()
    for(j in 1:26){
      var1 <- as.numeric(unlist(data[,covid_yesnos[[i]]]))
      var2 <- as.numeric(unlist(data[,covid_yesnos[[j]]]))
      if(similarity == "jaccard"){
        results.df[i,j] <- jaccardCalc(var1,var2)
      }
      else if(similarity == "hamming"){
        results.df[i,j] <- sum(var1 == var2) / length(var1)
      }
      
    }
    print(paste(round(100*i/26,0),"% complete"))
    t1 <- Sys.time()
    looptime <- t1-t0
    print(paste("Estimated time remaining:", round((26-i) * looptime, 2)))
    
  }
  return(results.df)
}

