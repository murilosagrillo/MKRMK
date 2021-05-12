source("functions.R") 

#------------------------------------------------------------------
# MK functions MK and RMK
#------------------------------------------------------------------

#################
#MK DISTRIBUTION#
#################

#cumulative distribution function of MK distribution with parameters alpha=0.2 and beta=0.3
dmk(0.3,0.2,0.3)

#probability density function of MK distribution with parameters alpha=0.2 and beta=0.3
pmk(q=0.5,alpha=0.2,beta=0.3)

#Quantile function of MK distribution with parameters alpha=0.2 and beta=0.3, for u=1/2 
qmk(u=1/2,alpha=0.2,beta=0.3)

#Generation of n=200 pseudo-random numbers of MK distribution with parameters alpha=0.2 and beta=0.3
pseudo_MK<-rmk(n=200,alpha=0.2,beta=0.3)

#MK fit in pseudo_MK data
mkfit(pseudo_MK)

#################
#RMK DISTRIBUTION#
#################

#cumulative distribution function of RMK distribution with parameters alpha=0.2 and beta=0.3
prmk(0.5,0.2,0.3)

#probability density function of RMK distribution with parameters alpha=0.2 and beta=0.3
drmk(0.7,0.2,0.3)

#Quantile function of RMK distribution with parameters alpha=0.2 and beta=0.3, for u=1/2 
qrmk(u=1/2,alpha=0.2,beta=0.3)

#Generation of n=200 pseudo-random numbers of RMK distribution with parameters alpha=0.2 and beta=0.3
pseudo_RMK<-rrmk(200,0.2,0.3)

#RMK fit in pseudo_RMK data
rmkfit(pseudo_RMK)

