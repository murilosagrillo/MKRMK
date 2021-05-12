source("functions.R") 

#------------------------------------------------------------------
# MK functions MK and RMK
#------------------------------------------------------------------

#################
#MK DISTRIBUTION#
#################

#cumulative distribution function of MK distribution with parameters alpha=0.2 and beta=0.5
dmk(y=0.3,alpha=0.2,beta=0.5)

#probability density function of MK distribution with parameters alpha=0.2 and beta=0.5
pmk(q=0.5,alpha=0.2,beta=0.5)

#Quantile function of MK distribution with parameters alpha=0.2 and beta=0.5, for u=1/2 
qmk(u=1/2,alpha=0.2,beta=0.5)

#Generation of n=200 pseudo-random numbers of MK distribution with parameters alpha=0.2 and beta=0.5
pseudo_MK<-rmk(n=200,alpha=0.2,beta=0.5)

#MK fit in pseudo_MK data
mkfit(pseudo_MK)

#################
#RMK DISTRIBUTION#
#################

#cumulative distribution function of RMK distribution with parameters alpha=0.2 and beta=0.5
prmk(z=0.5,alpha=0.2,beta=0.5)

#probability density function of RMK distribution with parameters alpha=0.2 and beta=0.5
drmk(q=0.7,alpha=0.2,beta=0.5)

#Quantile function of RMK distribution with parameters alpha=0.2 and beta=0.5, for u=1/2 
qrmk(u=1/2,alpha=0.2,beta=0.5)

#Generation of n=200 pseudo-random numbers of RMK distribution with parameters alpha=0.2 and beta=0.5
pseudo_RMK<-rrmk(n=200,alpha=0.2,beta=0.5)

#RMK fit in pseudo_RMK data
rmkfit(pseudo_RMK)

