library(tidyverse)
library(corrplot)
library(splines) # function ns()
library(rms) # function lsp()
library(mfp) # function mfp(), fp()

rm(list=ls())


############################################################################
## data prep                                                              ##
## for workshop only; in other situations you may want to do more here... ##
############################################################################

datX<-read.csv("../datasets/NHANESI_subset_X.csv")
datY<-read.csv("../datasets/NHANESI_subset_y.csv")

colnames(datX)[1]<-"pid"
colnames(datY)[1]<-"pid"

dat<-full_join(datX,datY,by="pid")
datComplete<-dat[rowSums(is.na(dat))==0,]

sumDat<-dat %>%
  dplyr::group_by(Age) %>%
  dplyr::summarise(
    meanBMI=mean(BMI),
    seBMI=sd(BMI)/sqrt(n()),
    .groups="drop"
  ) %>%
  dplyr::mutate(
    meanBMI_low=meanBMI-qnorm(0.975)*seBMI,
    meanBMI_upp=meanBMI+qnorm(0.975)*seBMI
  )


###############################################
## sample mean BMI + 95% CI per year of age. ##
###############################################

gFull<-ggplot() +
  geom_point(data=datComplete,mapping=aes(x=Age,y=BMI),alpha=0.5)+
  geom_ribbon(data=sumDat,mapping=aes(x=Age,ymin=meanBMI_low,ymax=meanBMI_upp),fill="darkgrey",alpha=0.5) +
  geom_line(data=sumDat,mapping=aes(x=Age,y=meanBMI),col="steelblue",lwd=1.25) +
  theme_light()

gZoom<-ggplot() +
  geom_point(data=datComplete,mapping=aes(x=Age,y=BMI),alpha=0.5)+
  geom_ribbon(data=sumDat,mapping=aes(x=Age,ymin=meanBMI_low,ymax=meanBMI_upp),fill="darkgrey",alpha=0.5) +
  geom_line(data=sumDat,mapping=aes(x=Age,y=meanBMI),col="steelblue",lwd=1.5) +
  ylim(c(23,28.5)) +
  theme_light()

gridExtra::grid.arrange(gFull,gZoom,nrow=1)


############################
## piecewise constant fit ##
############################

datComplete<-datComplete %>%
  dplyr::mutate(
    ageGrp2lvls=case_when(
      Age<48~"<48",
      Age>=48~"48+"
    ),
    ageGrp4lvls=case_when(
      Age<35~"<35",
      Age>=35 & Age<48~"35-47",
      Age>=48 & Age<66~"48-65",
      Age>=66~"66+"
    )
  )

## 25-47 and 48-74
datComplete<-datComplete %>%
  dplyr::mutate(
    BMImeanPiece2steps=case_when(
      ageGrp2lvls=="<48"~mean(datComplete$BMI[datComplete$ageGrp2lvls=="<48"]),
      ageGrp2lvls=="48+"~mean(datComplete$BMI[datComplete$ageGrp2lvls=="48+"]),
    ),
    BMIsePiece2steps=case_when(
      ageGrp2lvls=="<48"~sd(datComplete$BMI[datComplete$ageGrp2lvls=="<48"])/sqrt(sum(datComplete$ageGrp2lvls=="<48")),
      ageGrp2lvls=="48+"~sd(datComplete$BMI[datComplete$ageGrp2lvls=="48+"])/sqrt(sum(datComplete$ageGrp2lvls=="48+")),
    )
  ) %>%
  dplyr::mutate(
    BMIciLowPiece2steps=BMImeanPiece2steps-qnorm(0.975)*BMIsePiece2steps,
    BMIciUppPiece2steps=BMImeanPiece2steps+qnorm(0.975)*BMIsePiece2steps
  )

## 25-34, 35-47, 48-65, 66-74
datComplete<-datComplete %>%
  dplyr::mutate(
    BMImeanPiece4steps=case_when(
      ageGrp4lvls=="<35"~mean(datComplete$BMI[datComplete$ageGrp4lvls=="<35"]),
      ageGrp4lvls=="35-47"~mean(datComplete$BMI[datComplete$ageGrp4lvls=="35-47"]),
      ageGrp4lvls=="48-65"~mean(datComplete$BMI[datComplete$ageGrp4lvls=="48-65"]),
      ageGrp4lvls=="66+"~mean(datComplete$BMI[datComplete$ageGrp4lvls=="66+"])
    ),
    BMIsePiece4steps=case_when(
      ageGrp4lvls=="<35"~sd(datComplete$BMI[datComplete$ageGrp4lvls=="<35"])/sqrt(sum(datComplete$ageGrp4lvls=="<35")),
      ageGrp4lvls=="35-47"~sd(datComplete$BMI[datComplete$ageGrp4lvls=="35-47"])/sqrt(sum(datComplete$ageGrp4lvls=="35-47")),
      ageGrp4lvls=="48-65"~sd(datComplete$BMI[datComplete$ageGrp4lvls=="48-65"])/sqrt(sum(datComplete$ageGrp4lvls=="48-65")),
      ageGrp4lvls=="66+"~sd(datComplete$BMI[datComplete$ageGrp4lvls=="66+"])/sqrt(sum(datComplete$ageGrp4lvls=="66+")),
    )
  ) %>%
  dplyr::mutate(
    BMIciLowPiece4steps=BMImeanPiece4steps-qnorm(0.975)*BMIsePiece4steps,
    BMIciUppPiece4steps=BMImeanPiece4steps+qnorm(0.975)*BMIsePiece4steps
  )

datPred<-data.frame(
  Age=seq(25,75,length=100)
) %>%
  dplyr::mutate(
    Age2=Age^2,
    Age3=Age^3,
    Age4=Age^4,
    Age5=Age^5
  ) %>%
  dplyr::mutate(
    ageGrp2lvls=case_when(
      Age<48~"<48",
      Age>=48~"48+"
    ),
    ageGrp4lvls=case_when(
      Age<35~"<35",
      Age>=35 & Age<48~"35-47",
      Age>=48 & Age<66~"48-65",
      Age>=66~"66+"
    )
  )

gPiece2lvls<-datComplete %>%
  ggplot() +
  geom_point(mapping=aes(x=Age,y=BMI),alpha=0.5)+
  geom_ribbon(mapping=aes(x=Age,ymin=BMIciLowPiece2steps,ymax=BMIciUppPiece2steps),fill="darkgrey",alpha=0.5) +
  geom_line(data=sumDat,mapping=aes(x=Age,y=meanBMI),col="steelblue",lwd=1.5) +
  geom_line(mapping=aes(x=Age,y=BMImeanPiece2steps),col="orange",lwd=1.5) +
  ylim(c(23,28.5)) +
  theme_light()

gPiece4lvls<-datComplete %>%
  ggplot() +
  geom_point(mapping=aes(x=Age,y=BMI),alpha=0.5)+
  geom_ribbon(mapping=aes(x=Age,ymin=BMIciLowPiece4steps,ymax=BMIciUppPiece4steps),fill="darkgrey",alpha=0.5) +
  geom_line(data=sumDat,mapping=aes(x=Age,y=meanBMI),col="steelblue",lwd=1.5) +
  geom_line(mapping=aes(x=Age,y=BMImeanPiece4steps),col="orange",lwd=1.5) +
  ylim(c(23,28.5)) +
  theme_light()

gridExtra::grid.arrange(gPiece2lvls,gPiece4lvls,nrow=1)


#####################
## polynomial fit  ##
#####################

datComplete<-datComplete %>%
  dplyr::mutate(
    Age2=Age^2,
    Age3=Age^3,
    Age4=Age^4,
    Age5=Age^5
  )
modPoly3<-glm(BMI~Age+Age2+Age3,data=datComplete)
modPoly5<-glm(BMI~Age+Age2+Age3+Age4+Age5,data=datComplete)

datPred<-datPred %>%
  dplyr::mutate(
    BMIpred3=predict(modPoly3,newdata=datPred,se.fit=T)$fit,
    BMIpred3se=predict(modPoly3,newdata=datPred,se.fit=T)$se.fit,
    BMIpred5=predict(modPoly5,newdata=datPred,se.fit=T)$fit,
    BMIpred5se=predict(modPoly5,newdata=datPred,se.fit=T)$se.fit
  ) %>%
  dplyr::mutate(
    BMIpred3low=BMIpred3-qnorm(0.975)*BMIpred3se,
    BMIpred3upp=BMIpred3+qnorm(0.975)*BMIpred3se,
    BMIpred5low=BMIpred5-qnorm(0.975)*BMIpred5se,
    BMIpred5upp=BMIpred5+qnorm(0.975)*BMIpred5se,
  )

gPoly3<-datComplete %>%
  ggplot() +
  geom_point(mapping=aes(x=Age,y=BMI),alpha=0.5)+
  geom_ribbon(data=datPred,mapping=aes(x=Age,ymin=BMIpred3low,ymax=BMIpred3upp),fill="darkgrey",alpha=0.5) +
  geom_line(data=sumDat,mapping=aes(x=Age,y=meanBMI),col="steelblue",lwd=1.25) +
  geom_line(data=datPred,mapping=aes(x=Age,y=BMIpred3),col="orange",lwd=1.5) +
  ylim(c(23,28.5)) +
  theme_light()

gPoly5<-datComplete %>%
  ggplot() +
  geom_point(mapping=aes(x=Age,y=BMI),alpha=0.5)+
  geom_ribbon(data=datPred,mapping=aes(x=Age,ymin=BMIpred5low,ymax=BMIpred5upp),fill="darkgrey",alpha=0.5) +
  geom_line(data=sumDat,mapping=aes(x=Age,y=meanBMI),col="steelblue",lwd=1.25) +
  geom_line(data=datPred,mapping=aes(x=Age,y=BMIpred5),col="orange",lwd=1.5) +
  ylim(c(23,28.5)) +
  theme_light()

gridExtra::grid.arrange(gPoly3,gPoly5,nrow=1)


####################
## linear splines ##
####################

modLsp2Knts<-glm(BMI~lsp(Age,parms=c(40,65)),data=datComplete)
modLsp3Knts<-glm(BMI~lsp(Age,parms=c(35,50,65)),data=datComplete)

datPred<-datPred %>%
  dplyr::mutate(
    BMIpredls2=predict(modLsp2Knts,newdata=datPred,se.fit=T)$fit,
    BMIpredls2se=predict(modLsp2Knts,newdata=datPred,se.fit=T)$se.fit,
    BMIpredls3=predict(modLsp3Knts,newdata=datPred,se.fit=T)$fit,
    BMIpredls3se=predict(modLsp3Knts,newdata=datPred,se.fit=T)$se.fit
  ) %>%
  dplyr::mutate(
    BMIpredls2low=BMIpredls2-qnorm(0.975)*BMIpredls2se,
    BMIpredls2upp=BMIpredls2+qnorm(0.975)*BMIpredls2se,
    BMIpredls3low=BMIpredls3-qnorm(0.975)*BMIpredls3se,
    BMIpredls3upp=BMIpredls3+qnorm(0.975)*BMIpredls3se
  )

gPredLs2<-datComplete %>%
  ggplot() +
  geom_point(mapping=aes(x=Age,y=BMI),alpha=0.5)+
  geom_ribbon(data=datPred,mapping=aes(x=Age,ymin=BMIpredls2low,ymax=BMIpredls2upp),fill="darkgrey",alpha=0.5) +
  geom_line(data=sumDat,mapping=aes(x=Age,y=meanBMI),col="steelblue",lwd=1.25) +
  geom_line(data=datPred,mapping=aes(x=Age,y=BMIpredls2),col="orange",lwd=1.5) +
  ylim(c(23,28.5)) +
  theme_light()

gPredLs3<-datComplete %>%
  ggplot() +
  geom_point(mapping=aes(x=Age,y=BMI),alpha=0.5)+
  geom_ribbon(data=datPred,mapping=aes(x=Age,ymin=BMIpredls2low,ymax=BMIpredls3upp),fill="darkgrey",alpha=0.5) +
  geom_line(data=sumDat,mapping=aes(x=Age,y=meanBMI),col="steelblue",lwd=1.25) +
  geom_line(data=datPred,mapping=aes(x=Age,y=BMIpredls3),col="orange",lwd=1.5) +
  ylim(c(23,28.5)) +
  theme_light()

gridExtra::grid.arrange(gPredLs2,gPredLs3,nrow=1)


##############################
## restricted cubic splines ##
##############################

modNs5Knts1<-glm(BMI~ns(Age,knots=c(35,50,65),Boundary.knots=c(25,75)),data=datComplete)
modNs5Knts2<-glm(BMI~ns(Age,knots=c(40,52,65),Boundary.knots=c(35,75)),data=datComplete)

datPred<-datPred %>%
  dplyr::mutate(
    BMIpredns1=predict(modNs5Knts1,newdata=datPred,se.fit=T)$fit,
    BMIpredns1se=predict(modNs5Knts1,newdata=datPred,se.fit=T)$se.fit,
    BMIpredns2=predict(modNs5Knts2,newdata=datPred,se.fit=T)$fit,
    BMIpredns2se=predict(modNs5Knts2,newdata=datPred,se.fit=T)$se.fit
  ) %>%
  dplyr::mutate(
    BMIpredns1low=BMIpredns1-qnorm(0.975)*BMIpredls2se,
    BMIpredns1upp=BMIpredns1+qnorm(0.975)*BMIpredls2se,
    BMIpredns2low=BMIpredns2-qnorm(0.975)*BMIpredls3se,
    BMIpredns2upp=BMIpredns2+qnorm(0.975)*BMIpredls3se
  )

gPredNs1<-datComplete %>%
  ggplot() +
  geom_point(mapping=aes(x=Age,y=BMI),alpha=0.5)+
  geom_ribbon(data=datPred,mapping=aes(x=Age,ymin=BMIpredns1low,ymax=BMIpredns1upp),fill="darkgrey",alpha=0.5) +
  geom_line(data=sumDat,mapping=aes(x=Age,y=meanBMI),col="steelblue",lwd=1.25) +
  geom_line(data=datPred,mapping=aes(x=Age,y=BMIpredns1),col="orange",lwd=1.5) +
  ylim(c(23,28.5)) +
  theme_light()

gPredNs2<-datComplete %>%
  ggplot() +
  geom_point(mapping=aes(x=Age,y=BMI),alpha=0.5)+
  geom_ribbon(data=datPred,mapping=aes(x=Age,ymin=BMIpredns2low,ymax=BMIpredns2upp),fill="darkgrey",alpha=0.5) +
  geom_line(data=sumDat,mapping=aes(x=Age,y=meanBMI),col="steelblue",lwd=1.25) +
  geom_line(data=datPred,mapping=aes(x=Age,y=BMIpredns2),col="orange",lwd=1.5) +
  ylim(c(23,28.5)) +
  theme_light()

gridExtra::grid.arrange(gPredNs1,gPredNs2,nrow=1)


############################
## fractional polynomials ##
############################

modFP<-mfp(BMI~fp(Age),data=datComplete)

datPred<-datPred %>%
  dplyr::mutate(
    BMIpredfp=predict(modFP,newdata=datPred,se.fit=T)$fit,
    BMIpredfpse=predict(modFP,newdata=datPred,se.fit=T)$se.fit,
  ) %>%
  dplyr::mutate(
    BMIpredfplow=BMIpredfp-qnorm(0.975)*BMIpredfpse,
    BMIpredfpupp=BMIpredfp+qnorm(0.975)*BMIpredfpse
  )

gPredFp<-datComplete %>%
  ggplot() +
  geom_point(mapping=aes(x=Age,y=BMI),alpha=0.5)+
  geom_ribbon(data=datPred,mapping=aes(x=Age,ymin=BMIpredfplow,ymax=BMIpredfpupp),fill="darkgrey",alpha=0.5) +
  geom_line(data=sumDat,mapping=aes(x=Age,y=meanBMI),col="steelblue",lwd=1.25) +
  geom_line(data=datPred,mapping=aes(x=Age,y=BMIpredfp),col="orange",lwd=1.5) +
  ylim(c(23,28.5)) +
  theme_light()
