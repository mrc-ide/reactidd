rm(list=ls())

library(pROC)
library(colorspace)


### Functions

ExtrapolateByProportionToTest=function(proportions, tpr, prop_pop_to_test){
  delta=(proportions-prop_pop_to_test)
  delta_pos=delta_neg=delta
  delta_pos[delta_pos<0]=NA
  delta_neg[delta_neg>0]=NA
  ids=c(which.min(abs(delta_neg)), which.min(abs(delta_pos)))
  x=proportions[ids]
  y=tpr[ids]
  a=(y[2]-y[1])/(x[2]-x[1])
  b=y[1]-a*x[1]
  tmp=a*prop_pop_to_test+b
  if (is.na(tmp)){
    tmp=unique(y)
  }
  return(tmp)
}


### Loading predicted probabilities and test results

mydata=data.frame(readRDS("Data/Symptoms/predictions.rds"))


### Computing model performance (AUC)

# AUC in R2-R7
x=mydata[mydata$Data=="Test",]
myroc=roc(response=x$test_result_1, predictor=x$StabilityModelPreds)

# AUC in R8
x=mydata[mydata$Data=="Round_8",]
myroc=roc(response=x$test_result_1, predictor=x$StabilityModelPreds)


### Setting parameters

# Definition of age groups
age_groups=c("under_18", "from_18_to_54", "over_55")
names(age_groups)=c("Age 5-17", "Age 18-54", "Age 55+")

# Optimal proportions of tests 
prop_opt_stab_list=c(0.3474319316,0.3668497731)
names(prop_opt_stab_list)=c("Test", "Round_8")

# ID of where to stop (no more individuals with deleterious symptoms)
id_no_deleterious_sympt_1_of_4_list=c(2,2)
names(id_no_deleterious_sympt_1_of_4_list)=c("Test", "Round_8")
id_no_deleterious_sympt_logistic_4_list=c(5,5)
names(id_no_deleterious_sympt_logistic_4_list)=c("Test", "Round_8")
id_no_deleterious_sympt_stab_list=NULL
for (mysubset in c("Test", "Round_8")){
  mydata=data.frame(readRDS("Data/Symptoms/predictions.rds"))
  mydata=mydata[mydata$Data=="Test",]
  id_no_deleterious_sympt_stab_list=c(id_no_deleterious_sympt_stab_list, 
                                      max(which(cumsum(table(mydata$StabilityModelPreds,mydata$any_of_our_vars)[,"1"])==0))+1)
}
names(id_no_deleterious_sympt_stab_list)=c("Test", "Round_8")

# Figure format
pdf_output=FALSE


k=0
for (mysubset in c("Test", "Round_8")){
  print(mysubset)
  
  # Extracting prepared numbers
  prop_opt_stab=prop_opt_stab_list[mysubset]
  id_no_deleterious_sympt_1_of_4=id_no_deleterious_sympt_1_of_4_list[mysubset]
  id_no_deleterious_sympt_logistic_4=id_no_deleterious_sympt_logistic_4_list[mysubset]
  id_no_deleterious_sympt_stab=id_no_deleterious_sympt_stab_list[mysubset]
  
  
  ### Getting the age-specific thresholds in predicted probability
  
  mydata=data.frame(readRDS("Data/Symptoms/predictions.rds"))
  mydata=mydata[mydata$Data==mysubset,]
  
  # Definition of the figure name
  if (mysubset=="Test"){
    plotname="Figure3_A"
  } else {
    plotname="Figure3_B"
  }
  
  # Definition of the figure format
  if (pdf_output){
    pdf(paste0("Figures/COVID-19_symptoms/",plotname,".pdf"), width=11, height=12)
  } else {
    png(paste0("Figures/COVID-19_symptoms/",plotname,".png"), width=11, height=12, unit="in", res=200)
  }
  
  
  # For loop over full population and the three age groups
  par(mar=c(7,7,10,1))
  
  # Getting the corresponding subset of the data
  mydata=data.frame(readRDS("Data/Symptoms/predictions.rds"))
  mydata=mydata[mydata$Data==mysubset,]
  
  
  ### Computing the True Positive Rate (i.e. proportion of symptomatic positives that would be invited for a test)
  
  # Creating empty objects to be filled for all three models
  mytpr_1_of_4=mytpr_logistic_4=mytpr=NULL
  myproportion_1_of_4=myproportion_logistic_4=myproportion=NULL
  
  # Creating empty objects to be filled (specific to stability selection)
  myp_none_of_four=NULL 
  myp_none_of_four_n=NULL
  
  # Stability model
  proba_list=sort(unique(mydata$StabilityModelPreds))
  for (p in proba_list){
    tmppred=ifelse(mydata$StabilityModelPreds>=p, yes=1, no=0)
    tmppred=factor(tmppred, levels=c(0,1))
    cont=table(tmppred, mydata$test_result_1)
    myproportion=c(myproportion, (cont["1","Negative"]+cont["1","Positive"])/sum(cont))
    mytpr=c(mytpr, cont["1","Positive"]/(cont["1","Positive"]+cont["0","Positive"]))
    TP=cont["1","Positive"]
    
    # Proportion of TP with none of the four
    cont_4sympt=table(tmppred, mydata$one_of_four, mydata$test_result_1)
    cont_4sympt=cont_4sympt[,,"Positive"]
    myp_none_of_four=c(myp_none_of_four, cont_4sympt["1","0"]/(cont_4sympt["1","0"]+cont_4sympt["1","1"]))
    
    # Proportion of (TP+FP) with none of the four
    cont_4sympt_bis=table(tmppred, mydata$one_of_four)
    myp_none_of_four_n=c(myp_none_of_four_n, cont_4sympt_bis["1","0"]/(cont_4sympt_bis["1","0"]+cont_4sympt_bis["1","1"]))
  }
  myproportion=c(myproportion,0)
  mytpr=c(mytpr,0)
  
  # Remove last point (TPR=1) as the group with smallest predicted probabilities should not be sent for a test
  myproportion=myproportion
  mytpr=mytpr
  myp_none_of_four=c(myp_none_of_four,NA)
  myp_none_of_four=myp_none_of_four
  myp_none_of_four_n=myp_none_of_four_n
  
  # Storing subset and age-specific values
  assign(paste0("mytpr_",mysubset,"_",k), mytpr)
  assign(paste0("myproportion_",mysubset,"_",k), myproportion)
  assign(paste0("myp_none_of_four_",mysubset,"_",k), myp_none_of_four)
  assign(paste0("myp_none_of_four_n_",mysubset,"_",k), myp_none_of_four_n)
  
  # Logistic with four classic symptoms as predictors
  proba_list=sort(unique(mydata$fourSymptomModelPreds))
  for (p in proba_list){
    tmppred=ifelse(mydata$fourSymptomModelPreds>=p, yes=1, no=0)
    tmppred=factor(tmppred, levels=c(0,1))
    cont=table(tmppred, mydata$test_result_1)
    myproportion_logistic_4=c(myproportion_logistic_4, (cont["1","Negative"]+cont["1","Positive"])/sum(cont))
    mytpr_logistic_4=c(mytpr_logistic_4, cont["1","Positive"]/(cont["1","Positive"]+cont["0","Positive"]))
  }
  myproportion_logistic_4=c(myproportion_logistic_4,0)
  mytpr_logistic_4=c(mytpr_logistic_4,0)
  myproportion_logistic_4=myproportion_logistic_4
  mytpr_logistic_4=mytpr_logistic_4
  
  # At-least-one-of-four classic symptoms
  proba_list=sort(unique(mydata$one_of_four))
  for (p in proba_list){
    tmppred=ifelse(mydata$one_of_four>=p, yes=1, no=0)
    tmppred=factor(tmppred, levels=c(0,1))
    cont=table(tmppred, mydata$test_result_1)
    myproportion_1_of_4=c(myproportion_1_of_4, (cont["1","Negative"]+cont["1","Positive"])/sum(cont))
    mytpr_1_of_4=c(mytpr_1_of_4, cont["1","Positive"]/(cont["1","Positive"]+cont["0","Positive"]))
  }
  myproportion_1_of_4=c(myproportion_1_of_4,0)
  myproportion_1_of_4=myproportion_1_of_4
  mytpr_1_of_4=c(mytpr_1_of_4,0)
  mytpr_1_of_4=mytpr_1_of_4
  
  # Storing subset and age-specific values
  assign(paste0("mytpr_1_of_4_",mysubset,"_",k), mytpr_1_of_4)
  assign(paste0("myproportion_1_of_4_",mysubset,"_",k), myproportion_1_of_4)
  
  # Starting the figure: TPR as a function of the proportion of symptomatics tested
  plot(myproportion, mytpr, col="tomato", ylim=c(0,1), las=1, cex.lab=2,
       xlab="", ylab="", xlim=c(0,1),
       pch=19, cex=0.5, xaxt="n", yaxt="n")
  mtext("Proportion of symptomatic PCR positives detected", side=2, line=4, cex=1.5)
  mtext("Proportion of symptomatic individuals tested in England", side=1, line=4, cex=1.5)
  
  # Including X-axis
  for (tick in 1:11){
    axis(side=1, at=formatC(seq(0,1,by=0.1)[tick], format="f", digits=1), las=1, cex.axis=ifelse(k==0,yes=1.5,no=1))
  }
  
  # Including Y-axis
  axis(side=2, at=seq(0,1,by=0.1), las=1, cex.axis=ifelse(k==0,yes=1.5,no=1))
  
  # Adding the x-labels
  par(xpd=TRUE)
  
  # Adding the proportions of none-of-the-four classic symptoms among the TP for stability selection model
  hstart=1.07
  myscale=3
  polygon(x=c(max(myproportion),myproportion,rev(myproportion)),
          y=c(hstart,hstart+myp_none_of_four/myscale,
              rev(rep(hstart,length(myp_none_of_four)))), 
          col=adjustcolor(lighten("lightgrey", amount=0.5), alpha.f = 0.5), border=NA)
  
  for (b in 1:length(myproportion)){
    segments(x0=myproportion[b], x1=myproportion[b],
             y0=hstart, y1=hstart+myp_none_of_four[b]/myscale, 
             lwd=0.5, col="grey")
  }
  
  segments(x0=-0.04,x1=1.04,y0=hstart,y1=hstart,lty=1)
  segments(x0=-0.04,x1=-0.04,y0=hstart,y1=hstart+0.5/myscale,lty=1)
  segments(x0=-0.04,x1=-0.05,y0=hstart,y1=hstart,lty=1)
  text(-0.05, hstart, labels="0.0", pos=2, cex=1.5)
  segments(x0=-0.04,x1=-0.05,y0=hstart+0.5/myscale,y1=hstart+0.5/myscale,lty=1)
  text(-0.05, hstart+0.5/myscale, labels="0.5", pos=2, cex=1.5)
  
  for (tick in seq(0.1,0.4,by=0.1)){
    segments(x0=-0.04,x1=-0.05,y0=hstart+tick/myscale,y1=hstart+tick/myscale,lty=1)
    text(-0.05, hstart+tick/myscale, labels=tick, pos=2, cex=1.5)
  }
  
  text(-0.152,y=hstart+0.25/myscale, labels="None-of-four", srt=90, cex=1.5)
  par(xpd=FALSE)
  
  # Adding segments behind points from stability selection model (link to top)
  for (b in 1:length(myproportion)){
    abline(v=myproportion[b], col=lighten("lightgrey", amount=0.3))
  }
  
  # Adding the points (TPR by proportion of symptomatics to test)
  points(myproportion, mytpr, col="tomato", pch=19, cex=0.5)
  points(myproportion_1_of_4, mytpr_1_of_4, col="navy", pch=19, cex=0.5)
  points(myproportion_logistic_4, mytpr_logistic_4, col="forestgreen", pch=19, cex=0.5)
  
  # Extrapolating between the points for all three models (up to all individuals with any of the deleterious symptoms)
  lines(myproportion_1_of_4[id_no_deleterious_sympt_1_of_4:length(myproportion_1_of_4)], 
        mytpr_1_of_4[id_no_deleterious_sympt_1_of_4:length(myproportion_1_of_4)], col="navy") # at-least-one-of-four strategy
  lines(myproportion_logistic_4[id_no_deleterious_sympt_logistic_4:length(myproportion_logistic_4)],
        mytpr_logistic_4[id_no_deleterious_sympt_logistic_4:length(myproportion_logistic_4)], col="forestgreen") # logistic with four classic symptoms as predictors
  lines(myproportion[id_no_deleterious_sympt_stab:length(myproportion)],
        mytpr[id_no_deleterious_sympt_stab:length(myproportion)], col="tomato") # stability selection
  
  # Extrapolating between the points for all three models (up to all symptomatic individuals being tested)
  lines(myproportion_1_of_4[1:id_no_deleterious_sympt_1_of_4], 
        mytpr_1_of_4[1:id_no_deleterious_sympt_1_of_4], col="navy", lty=2) # at-least-one-of-four strategy
  lines(myproportion_logistic_4[1:id_no_deleterious_sympt_logistic_4], 
        mytpr_logistic_4[1:id_no_deleterious_sympt_logistic_4], col="forestgreen", lty=2) # logistic with four classic symptoms as predictors
  lines(myproportion[1:id_no_deleterious_sympt_stab], 
        mytpr[1:id_no_deleterious_sympt_stab], lty=2, col="tomato") # stability selection
  
  # Adding dotted lines for optimal proportion of symptomatics tested with stability selection 
  myproportion=eval(parse(text=paste0("myproportion_",mysubset,"_",k)))
  mytpr=eval(parse(text=paste0("mytpr_",mysubset,"_",k)))
  abline(v=myproportion[which.min(abs(myproportion-prop_opt_stab))], lty=3, col="tomato")
  abline(h=mytpr[which.min(abs(myproportion-prop_opt_stab))], lty=3, col="tomato")
  
  # Adding dotted lines for optimal proportion of symptomatics with any of the four tested
  abline(v=myproportion_1_of_4[2], lty=3, col="navy")
  abline(h=mytpr_1_of_4[2], lty=3, col="navy")
  
  if (mysubset=="Test"){
    # Adding legend
    legend("bottomright", lty=c(rep(1,3),3,2), pch=c(rep(19,3),NA,NA), pt.cex=0.5, cex=1, 
           col=c("navy","forestgreen","tomato","black","black"),
           legend=c("At-least-one-of-four classic symptoms", 
                    "Logistic model with four classic symptoms", 
                    "Stability selection", 
                    "Optimal performance",
                    "No deleterious symptoms"), bty="n")
  }
  box()
  dev.off()
}


### Extracting the True Positive Rate (TPR) for optimal proportions of symptomatic individuals to test

mysummary=NULL
for (mysubset in c("Test", "Round_8")){
  print(mysubset)
  
  # At-least-one-of-four strategy
  prop_opt_1_of_4=eval(parse(text=paste0("myproportion_1_of_4_",mysubset,"_0")))[2]
  tpr_1_of_4=eval(parse(text=paste0("mytpr_1_of_4_",mysubset,"_0")))[2]
  mysummary=rbind(mysummary, c(prop_opt_1_of_4, tpr_1_of_4))
  
  # Optimal proportion of symptomatic individuals to test for stability selection
  prop_opt_stab=prop_opt_stab_list[mysubset]
  
  # Extracting corresponding True Positive Rate with stability selection
  tpr_opt=ExtrapolateByProportionToTest(proportions=eval(parse(text=paste0("myproportion_",mysubset,"_0"))), 
                                        tpr=eval(parse(text=paste0("mytpr_",mysubset,"_0"))), 
                                        prop_pop_to_test=prop_opt_stab)
  mysummary=rbind(mysummary, c(prop_opt_stab, tpr_opt))
}

myp_none_of_four_Test_0[which.min(abs(myproportion_Test_0-prop_opt_stab_list["Test"]))]
myp_none_of_four_Round_8_0[which.min(abs(myproportion_Round_8_0-prop_opt_stab_list["Round_8"]))]

