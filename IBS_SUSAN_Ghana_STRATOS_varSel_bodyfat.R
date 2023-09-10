library(MASS) # needed for mvrnorm() and stepAIC()
library(glmnet)
library(tidyverse)

rm(list=ls())

#############
# load data #
#############

load("../datasets/bodyfat.rda") # object bodyfat
bodyfat<-bodyfat %>%
  dplyr::select(!density) %>% # virtually linear relationship with 'fat'
  dplyr::mutate(bmi=weight/(height/100)^2)


######################
# variable selection #
######################

# backward and forward selection (based on AIC; p=0.157)
modFull<-lm(fat~.,data=bodyfat,x=TRUE)
modEmpty<-lm(fat~height+abdomen,data=bodyfat) # we fix height and abdomen in the model

modBwd<-stepAIC(modFull,direction="backward",trace=FALSE,scope=list(upper=modFull,lower=modEmpty)) # specify k=qchisq(0.05,1,lower.tail=F) for p=0.05 threshold
summary(modBwd)
modFwd<-stepAIC(modEmpty,direction="forward",trace=FALSE,scope=list(upper=modFull,lower=modEmpty)) # specify k=qchisq(0.05,1,lower.tail=F) for p=0.05 threshold
summary(modFwd)

# lasso
cv<-cv.glmnet(x=modFull$x[,-1],y=bodyfat$fat)
modLasso<-glmnet(x=modFull$x[,-1],y=bodyfat$fat,lambda=cv$lambda.min,penalty.factor=c(1,1,0,1,1,0,1,1,1,1,1,1,1)) # fix height and abdomen to be in the model
modLassoCoefs<-predict(modLasso, type="coefficients")
modLassoCoefs<-setdiff(names(modLassoCoefs[,1])[modLassoCoefs[,1]!=0],"(Intercept)")

# univariate selection (p=0.05)
modUniCoefs<-character(0)
for(k in 2:ncol(modFull$x)){
  modTmp<-lm(bodyfat$fat~modFull$x[,k])
  if(summary(modTmp)$coeff[2,"Pr(>|t|)"]<0.05 | colnames(modFull$x)[k] %in% c("height","abdomen")){modUniCoefs<-c(modUniCoefs,colnames(modFull$x)[k])} # fix height and abdomen to be in the model
}


#############################
## bootstrapping the above ##
#############################

B<-1000
allVars<-names(modFull$coefficients)[-1]
resMatBwd<-matrix(0,nrow=B,ncol=length(allVars))
colnames(resMatBwd)<-allVars
resMatFwd<-resMatBwd
resMatLasso<-resMatBwd
resMatUni<-resMatBwd

for(b in 1:B){
  bodyfatBt<-bodyfat[sample(1:nrow(bodyfat),size=nrow(bodyfat),replace=TRUE),]
  
  modFullBt<-lm(fat~.,data=bodyfatBt,x=TRUE)
  modEmptyBt<-lm(fat~height+abdomen,data=bodyfatBt) # we fix height and abdomen in the model
  
  modBwd<-stepAIC(modFullBt,direction="backward",trace=FALSE,scope=list(upper=modFullBt,lower=modEmptyBt)) # specify k=qchisq(0.05,1,lower.tail=F) for p=0.05 threshold
  modFwd<-stepAIC(modEmptyBt,direction="forward",trace=FALSE,scope=list(upper=modFullBt,lower=modEmptyBt)) # specify k=qchisq(0.05,1,lower.tail=F) for p=0.05 threshold
  
  cv<-cv.glmnet(x=modFullBt$x[,-1],y=bodyfat$fat,)
  modLassoBt<-glmnet(x=modFullBt$x[,-1],y=bodyfat$fat,lambda=cv$lambda.min,penalty.factor=c(1,1,0,1,1,0,1,1,1,1,1,1,1)) # fix height and abdomen to be in the model
  modLassoCoefsBt<-predict(modLassoBt, type="coefficients")
  modLassoCoefsBt<-setdiff(names(modLassoCoefsBt[,1])[modLassoCoefsBt[,1]!=0],"(Intercept)")
  
  modUniCoefsBt<-character(0)
  for(k in 2:ncol(modFullBt$x)){
    modTmp<-lm(bodyfat$fat~modFullBt$x[,k])
    if(summary(modTmp)$coeff[2,"Pr(>|t|)"]<0.05 | colnames(modFull$x)[k] %in% c("height","abdomen")){modUniCoefsBt<-c(modUniCoefsBt,colnames(modFullBt$x)[k])} # fix height and abdomen to be in the model
  }
  
  resMatBwd[b,allVars %in% names(modBwd$coefficients)]<-1
  resMatFwd[b,allVars %in% names(modFwd$coefficients)]<-1
  resMatLasso[b,allVars %in% modLassoCoefsBt]<-1
  resMatUni[b,allVars %in% modUniCoefsBt]<-1
}


####################################################
## summarising selection probabilities in a table ##
####################################################

resDf<-data.frame(
  variable=rownames(summary(modFull)$coefficients),
  coefFullModel=summary(modFull)$coefficients[,"Estimate"],
  seFullModel=summary(modFull)$coefficients[,"Std. Error"],
  coefSelectedModel=summary(modBwd)$coefficients[match(rownames(summary(modFull)$coefficients),rownames(summary(modBwd)$coefficients)),"Estimate"],
  seSelectedModel=summary(modBwd)$coefficients[match(rownames(summary(modFull)$coefficients),rownames(summary(modBwd)$coefficients)),"Std. Error"],
  probSelectionUni=c("100.0% (fixed)",paste(sep="",format(nsmall=1,100*colMeans(resMatUni)),"%")),
  probSelectionFwd=c("100.0% (fixed)",paste(sep="",format(nsmall=1,100*colMeans(resMatFwd)),"%")),
  probSelectionBwd=c("100.0% (fixed)",paste(sep="",format(nsmall=1,100*colMeans(resMatBwd)),"%")),
  probSelectionLasso=c("100.0 (fixed)%",paste(sep="",format(nsmall=1,100*colMeans(resMatLasso)),"%"))
)
resDf<-resDf[order(c(1,colMeans(resMatBwd)),decreasing=TRUE),]

resDf$probSelectionUni[resDf$variable %in% c("height","abdomen")]<-gsub(resDf$probSelectionUni[resDf$variable %in% c("height","abdomen")],pattern="%",replacement="% (fixed)")
resDf$probSelectionFwd[resDf$variable %in% c("height","abdomen")]<-gsub(resDf$probSelectionFwd[resDf$variable %in% c("height","abdomen")],pattern="%",replacement="% (fixed)")
resDf$probSelectionBwd[resDf$variable %in% c("height","abdomen")]<-gsub(resDf$probSelectionBwd[resDf$variable %in% c("height","abdomen")],pattern="%",replacement="% (fixed)")
resDf$probSelectionLasso[resDf$variable %in% c("height","abdomen")]<-gsub(resDf$probSelectionLasso[resDf$variable %in% c("height","abdomen")],pattern="%",replacement="% (fixed)")

resDf %>%
  knitr::kable(digits=4,row.names=FALSE,col.names=c("",rep(c("Estimate","SE"),2),c("Univariable","Forward selection","Backward elimination","Lasso"))) %>%
  kableExtra::kable_styling(full_width=FALSE) %>%
  kableExtra::add_header_above(c(" ","Full model" = 2,"Selected model (backward elimination)" = 2, "Selection probabilities" = 4))




