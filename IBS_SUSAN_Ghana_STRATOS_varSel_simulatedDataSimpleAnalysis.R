library(MASS) # needed for mvrnorm() and stepAIC()
library(glmnet)
library(tidyverse)

rm(list=ls())


#################
# simulate data #
#################

set.seed(123)

## a few Gaussian random variables
simData<-function(k=12,n=500){
  df<-mvrnorm(n=n,mu=rnorm(k,sd=10),Sigma=rWishart(1,k,diag(k))[,,1])
  colnames(df)<-paste(sep="","n",1:k)
  df<-as.data.frame(df)
  
  ## a Poisson, a gamma and a categorical random variable
  df<-df %>%
    dplyr::mutate(
      p1=rpois(n=n,lambda=2),
      e1=rgamma(n=n,rate=1,shape=2),
      c1=sample(size=n,paste(sep="","lvl",1:5),prob=c(0.5,0.25,0.125,0.08,0.045),replace=TRUE)
    )
  
  ## generate the response variable
  df<-df %>%
    dplyr::mutate(
      y=n1+0.5*n2+0.25*n3+0.1*n4+0.05*p1+0.2*e1+ifelse(c1=="lvl2",5,0)+rnorm(n=n,sd=1)
      # n1, n2, n3, n4, p1, e1 are predictors as is level "lvl12" from c1
      # change the sd of the noise term to see effect on variable selection
    )
  
  return(df)
}



######################
# variable selection #
######################

# true set
corVars<-c("n1","n2","n3","n4","p1","e1","c1lvl2")

# backward and forward selection (based on AIC; p=0.157)
modFull<-lm(y~.,data=df,x=TRUE)
modEmpty<-lm(y~1,data=df)

modBwd<-stepAIC(modFull,direction="backward",trace=FALSE,scope=list(upper=modFull,lower=modEmpty)) # specify k=qchisq(0.05,1,lower.tail=F) for p=0.05 threshold
summary(modBwd)
modFwd<-stepAIC(modEmpty,direction="forward",trace=FALSE,scope=list(upper=modFull,lower=modEmpty)) # specify k=qchisq(0.05,1,lower.tail=F) for p=0.05 threshold
summary(modFwd)

# lasso
cv<-cv.glmnet(x=modFull$x[,-1],y=df$y)
modLasso<-glmnet(x=modFull$x[,-1],y=df$y,lambda=cv$lambda.min)
modLassoCoefs<-predict(modLasso, type="coefficients")
modLassoCoefs<-setdiff(names(modLassoCoefs[,1])[modLassoCoefs[,1]!=0],"(Intercept)")

# univariate selection (p=0.05)
modUniCoefs<-character(0)
for(k in 2:ncol(modFull$x)){
  modTmp<-lm(df$y~modFull$x[,k])
  if(summary(modTmp)$coeff[2,"Pr(>|t|)"]<0.05){modUniCoefs<-c(modUniCoefs,colnames(modFull$x)[k])}
}



#################################
## simulating many times       ##
## as an alternative: boostrap ##
#################################

B<-1000
allVars<-names(modFull$coefficients)[-1]
resMatBwd<-matrix(0,nrow=B,ncol=length(allVars))
colnames(resMatBwd)<-allVars
resMatFwd<-resMatBwd
resMatLasso<-resMatBwd
resMatUni<-resMatBwd

for(b in 1:B){
  dfBt<-simData(k=12,n=500)
  #dfBt<-df[sample(1:nrow(df),size=nrow(df),replace=TRUE),] # bootstrap alternative
  
  # forward selection and backward elimination
  modFullBt<-lm(y~.,data=dfBt,x=TRUE)
  modEmptyBt<-lm(y~1,data=dfBt)
  modBwd<-stepAIC(modFullBt,direction="backward",trace=FALSE,scope=list(upper=modFullBt,lower=modEmptyBt)) # specify k=qchisq(0.05,1,lower.tail=F) for p=0.05 threshold
  modFwd<-stepAIC(modEmptyBt,direction="forward",trace=FALSE,scope=list(upper=modFullBt,lower=modEmptyBt)) # specify k=qchisq(0.05,1,lower.tail=F) for p=0.05 threshold
  
  # lasso
  cv<-cv.glmnet(x=modFullBt$x[,-1],y=df$y)
  modLassoBt<-glmnet(x=modFullBt$x[,-1],y=df$y,lambda=cv$lambda.min)
  modLassoCoefsBt<-predict(modLassoBt, type="coefficients")
  modLassoCoefsBt<-setdiff(names(modLassoCoefsBt[,1])[modLassoCoefsBt[,1]!=0],"(Intercept)")
  
  # univariate with p=0.05
  modUniCoefsBt<-character(0)
  for(k in 2:ncol(modFullBt$x)){
    modTmp<-lm(df$y~modFullBt$x[,k])
    if(summary(modTmp)$coeff[2,"Pr(>|t|)"]<0.05){modUniCoefsBt<-c(modUniCoefsBt,colnames(modFullBt$x)[k])}
  }
  
  resMatBwd[b,allVars %in% names(modBwd$coefficients)]<-1
  resMatFwd[b,allVars %in% names(modFwd$coefficients)]<-1
  resMatLasso[b,allVars %in% modLassoCoefsBt]<-1
  resMatUni[b,allVars %in% modUniCoefsBt]<-1
}


##################################
## compute performance measures ##
##################################

corVarsBin<-rep(0,length(names(modFull$coefficients))-1)
corVarsBin[allVars %in% corVars]<-1

## proportions of correct variable sets being identified
propCorrectBwd<-mean(apply(MARGIN=1,X=resMatBwd,function(x){sum(x==corVarsBin)==length(corVarsBin)}))
propCorrectFwd<-mean(apply(MARGIN=1,X=resMatFwd,function(x){sum(x==corVarsBin)==length(corVarsBin)}))
propCorrectLas<-mean(apply(MARGIN=1,X=resMatLasso,function(x){sum(x==corVarsBin)==length(corVarsBin)}))
propCorrectUni<-mean(apply(MARGIN=1,X=resMatUni,function(x){sum(x==corVarsBin)==length(corVarsBin)}))

## proportion of times the correct variables are all included (even if some irrelevant ones are included)
propNomissBwd<-mean(apply(MARGIN=1,X=resMatBwd,function(x){sum(x==1 & x==corVarsBin)==sum(corVarsBin==1)}))
propNomissFwd<-mean(apply(MARGIN=1,X=resMatFwd,function(x){sum(x==1 & x==corVarsBin)==sum(corVarsBin==1)}))
propNomissLas<-mean(apply(MARGIN=1,X=resMatLasso,function(x){sum(x==corVarsBin)==length(corVarsBin)}))
propNomissUni<-mean(apply(MARGIN=1,X=resMatUni,function(x){sum(x==corVarsBin)==length(corVarsBin)}))

## proportion of times no irrelevant variables are included (even if some relevant ones are missed)
propNoirrBwd<-mean(apply(MARGIN=1,X=resMatBwd,function(x){sum(x==0 & x==corVarsBin)==sum(corVarsBin==0)}))
propNoirrFwd<-mean(apply(MARGIN=1,X=resMatFwd,function(x){sum(x==0 & x==corVarsBin)==sum(corVarsBin==0)}))
propNoirrLas<-mean(apply(MARGIN=1,X=resMatLasso,function(x){sum(x==0 & x==corVarsBin)==sum(corVarsBin==0)}))
propNoirrUni<-mean(apply(MARGIN=1,X=resMatUni,function(x){sum(x==0 & x==corVarsBin)==sum(corVarsBin==0)}))

## summarising this in a table
resDf<-data.frame(
  method=c("Univariable","Forward selection","Backward elimination","Lasso"),
  propCorrectSet=c(propCorrectUni,propCorrectFwd,propCorrectBwd,propCorrectLas),
  propNoMissSet=c(propNomissUni,propNomissFwd,propNomissBwd,propNomissLas),
  propNoIrrSet=c(propNoirrUni,propNoirrFwd,propNoirrBwd,propNoirrLas)
)

for(j in 2:4){
  resDf[,j]<-paste(sep="",format(nsmall=1,100*resDf[,j]),"%")
}

resDf %>%
  knitr::kable(align=c("l",rep("c",3)),col.names=c("","% correct set","% containing true set","% no irrelevant variables")) %>%
  kableExtra::kable_styling(full_width=FALSE)
