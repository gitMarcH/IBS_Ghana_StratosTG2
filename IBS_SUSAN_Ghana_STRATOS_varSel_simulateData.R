library(tidyverse)

#################
# simulate data #
#################

## a few Gaussian random variables
k<-12
n<-500
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
