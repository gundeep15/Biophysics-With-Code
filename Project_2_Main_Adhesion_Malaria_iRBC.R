##########################
#### SINGLE BOND #########
##########################

############ SLIP BOND ########

### Reading single slip data file ###
single.slip.data = read.table("single_slip_bena.txt", header = T)

### beta = 1/ KbT | T = 298K | Kb = Boltzmann Constant ###
beta = 1/4.114; ## in Units : 1/(nm. pN) ##


### NLL Function with Gaussian Noise ###
single.slip.func.norm = function( x0, k.s, standDev, Force, ts ){
  Foo = -1*Force*(10^12)  ## Converting force from N to pN ## 
  tau = (1/k.s)*(exp (-beta*Foo*x0)) ## expected value of lifetime tau ##
  -sum( dnorm (ts, mean = tau, sd = standDev, log = TRUE ))
}

### NLL Function with Exponential Noise ###
single.slip.func.expo = function( x0, k.s, Force, ts ){
  Foo = -1*Force*(10^12)   
  tau = (1/k.s)*(exp (-beta*Foo*x0))
  -sum( dexp (ts, rate = 1/tau, log = TRUE ))
}

library(bbmle)

## Running the two NLL functions ##
slip.mle.norm = mle2(single.slip.func.norm, start=list (x0 =0.05, k.s = 0.3, standDev = 1), data = list(Force = single.slip.data$Force , ts = single.slip.data$Lifetime), method = "L-BFGS-B", lower =0.000001  )
slip.mle.expo = mle2(single.slip.func.expo, start=list (x0 =0.05, k.s = 0.3), data = list(Force = single.slip.data$Force , ts = single.slip.data$Lifetime), method = "L-BFGS-B" , lower = 0.000001)

## Comparing the two models ##
AICtab(slip.mle.expo,slip.mle.norm)

## > AICtab(slip.mle.expo,slip.mle.norm)
##                 dAIC  df
## slip.mle.expo   0.0   2 
## slip.mle.norm   625.1 3 
## Therefore, for the rest of the analysis we will use exponential noise
## This result is also aligned with the experimental setup
## as we expect the probability of detatchment to be the same at all time intervals, i.e constant rate.


### Parameter Estimates: x0, and k.s ###
coef(slip.mle.expo)

# > coef(slip.mle.expo)
#         x0        k.s 
#  0.05636251 0.30783295 


### Confidence interval of param 
confint(profile(slip.mle.expo))

## > confint(profile(slip.mle.expo))
##        2.5 %     97.5 %
## x0  0.01457108 0.09743211
## k.s 0.26655931 0.35478147


## Profile plots ##
## Note : For an unknown reason, we were not getting the proper x-axes in the profile plots
## Hence, we made two separate for each of the two parameters.
par(mfrow=c(2,2))
plot(profile(slip.mle.expo), xlim = c(0.01,0.1))
plot(profile(slip.mle.expo))

## Plotting ##

par(mfrow=c(1,2))

#ALL Data#
plot(-1*single.slip.data$Force*10^12, single.slip.data$Lifetime, xlab = 'Force (pN)', ylab = 'Lifetime (s) ')

#Model Estimate#
x0.mle = coef(slip.mle.expo)[1]
k.s.mle = coef(slip.mle.expo)[2]
F.vec = xvec = seq(0, 50, by = 0.1)   
lines( F.vec, (1/k.s.mle)*(exp (-beta*F.vec*x0.mle)), type = "l", lty = 1, col = 'Red' )

#Binning : To regenerate Paper Fig, with average lifetimes of bins + standard error of means#
#Need package for binning : install.packages("fields")#
library(fields)
out<-stats.bin(single.slip.data$Force, single.slip.data$Lifetime, N = 20)

#eliminating non-significant bins#
fvec = -1*out$centers[-(1:10)]*10^12
avg = out$stats["mean",][-(1:10)]
sdev = out$stats["Std.Dev.",][-(1:10)]
num = out$stats["N",][-(1:10)]

#standard error of means#
sem = sdev/ sqrt(num)

#plotting bin data with error bars#
plot(fvec, avg, xlab = 'Force (pN)', ylab = 'Lifetime (s)')
sem = sdev/ sqrt(num)
arrows(fvec, avg-sem, fvec, avg+sem, length=0.05, angle=90, code=3)

#Model estimate#
x0.mle = coef(slip.mle.expo)[1]
k.s.mle = coef(slip.mle.expo)[2]

F.vec = xvec = seq(0, 50, by = 0.1)   
lines( F.vec, (1/k.s.mle)*(exp (-beta*F.vec*x0.mle)), type = "l", lty = 1, col = 'Red' )


############ CATCH BOND ###############

#Reading Data#
single.catch.data = read.table("single_catch_bena.txt", header = T)
beta = 1/4.114;

#creating NLL function with exponential error#
single.catch.func.expo = function( ksi, k.c, Force, F0,ts ){
  Foo = -1*Force*(10^12)   
  tau = (1/k.c)*(exp (-beta*((Foo-F0)^2)*ksi))
  -sum( dexp (ts, rate = 1/tau, log = TRUE ))
}

#running mle analysis#
catch.mle.expo = mle2(single.catch.func.expo, start=list (ksi =0.02, k.c = 0.32, F0 = 14), data = list(Force = single.catch.data$Force, ts = single.catch.data$Lifetime), method = "L-BFGS-B", lower = 0.000001 )

### Parameter Estimates : ksi, k.c, F0 ###
coef(catch.mle.expo)

# ksi         k.c          F0 
# 0.02757103  0.31686357 13.79788543 

# Confidence intervals of param #

confint(profile(catch.mle.expo))
#       2.5 %       97.5 %
#  ksi  0.01702403  0.03792804
#  k.c  0.28975044  0.34584268
#  F0   12.87507739 15.01138668


plot(profile(catch.mle.expo))

## Plotting ##

par(mfrow=c(1,2))

#ALL DATA#
plot(-1*single.catch.data$Force*10^12, single.catch.data$Lifetime,xlab = 'Force (pN)', ylab = 'Lifetime (s)')

ksi.mle =coef(catch.mle.expo)[1]
k.c.mle = coef(catch.mle.expo)[2]
F0.mle = coef(catch.mle.expo)[3]

F.vec = xvec = seq(0, 50, by = 0.1) 
lines( F.vec,  (1/k.c.mle)*(exp (-beta*((F.vec - F0.mle)^2)*ksi.mle)), type = "l", lty = 1, col = 'Red' )

# Binning : To regenerate Paper Fig, with average lifetimes of bins + standard error of means #

out2<-stats.bin(single.catch.data$Force, single.catch.data$Lifetime, N = 20)

#eliminating non-significant bins#
fvec2 = -1*out2$centers*10^12
avg2 = out2$stats["mean",]
sdev2 = out2$stats["Std.Dev.",]
num2 = out2$stats["N",]

#standard error of means#
sem2 = sdev2/ sqrt(num2)

#plotting bin data with error bars#
plot(fvec2, avg2, xlab = 'Force (pN)', ylab = 'Lifetime (s)')
arrows(fvec2, avg2-sem2, fvec2, avg2+sem2, length=0.05, angle=90, code=3)

#Model estimate#

ksi.mle =coef(catch.mle.expo)[1]
k.c.mle = coef(catch.mle.expo)[2]
F0.mle = coef(catch.mle.expo)[3]

F.vec = xvec = seq(0, 50, by = 0.1) 
lines( F.vec,  (1/k.c.mle)*(exp (-beta*((F.vec - F0.mle)^2)*ksi.mle)), type = "l", lty = 1, col = 'Red' )


### FINAL TEST ON SINGLE MOLECULE DATA : 
#Fitting slip data on slip and catch model, running AIC# 
#Fitting catch data on slip and catch model, running AIC# 

## Is SLIP AND CATCH Justified ##

slip.mle.catch.data = mle2(single.slip.func.expo, start=list (x0 =0.05, k.s = 0.3), data = list(Force = single.catch.data$Force , ts = single.catch.data$Lifetime ), method = "L-BFGS-B" )
catch.mle.slip.data = mle2(single.catch.func.expo, start=list (ksi =0.02, k.c = 0.32, F0 = 14 ), data = list(Force = single.slip.data$Force, ts = single.slip.data$Lifetime), method = "L-BFGS-B" )


AICtab(slip.mle.expo , catch.mle.slip.data)
#                     dAIC df
# slip.mle.expo       0.0  2 
# catch.mle.slip.data 6.9  3 
#Hence slip model is justified for slip data#


AICtab(catch.mle.expo , slip.mle.catch.data )
#                     dAIC df
# catch.mle.expo       0.0 3 
# slip.mle.catch.data 23.4 2 
#Hence catch model is justified for catch data#





#################################
####### MULTIPLE BONDS ##########
#################################

################# SLIP BOND #############

#Reading Data
multiple.slip.data = read.table("multiple_slip_bena.txt", header = T)
beta = 1/4.114;
library(bbmle)

# Known : Domain of N = 2-10
# Creating a for loop from N 2 to 10 - Run mle for each N
# Compare nll for optimum N

# Initializing nll vector, and mle object output#

nll.n.slip = rep(0, 9)
mle.n.slip= c(new("mle2"),new("mle2"),new("mle2"),new("mle2"),new("mle2"),new("mle2"),new("mle2"),new("mle2"),new("mle2"))

for(i in 2:10){
  multiple.slip.func.expo4 = function(sheer, kon.s,ts){
    sigabc = 2.39;
    radius = 3.5;
    N=i;
    sheer.foo = sheer
    k.s = 0.3095; #Estimated from single bond study
    x0= 0.0613;  #Estimated from single bond study
    force.foo = ((28.4*pi*radius^3)/(sigabc))*sheer.foo; # in pN if sheer in Pa
    tau = (1/k.s) * (factorial(N-1)/N) * (kon.s/k.s)^(N-1) * exp( -beta * force.foo* x0* sum(1/1:N))
    -sum( dexp (ts, rate = 1/tau, log = TRUE ))
  }
  slip.mle.expo = mle2(multiple.slip.func.expo4, start=list (kon.s = 10),
                       data = list(sheer = multiple.slip.data$SheerStress, 
                                   ts = multiple.slip.data$Lifetime), method = "L-BFGS-B", lower = c( kon.s = 0.000001) )
  kon = coef(slip.mle.expo)  
  print (i)
  print (kon)
  nll.n.slip[i-1] = -logLik(slip.mle.expo)
  mle.n.slip[i-1] = slip.mle.expo
}

#Running AIC analysis for finding optimum N
AICtab(mle.n.slip, weights = TRUE)

#N dAIC df weight
#2 model1  0.0 1  0.8205 <<<<----- OPTIMUM N = 2 for slip
#3 model2  3.2 1  0.1628
#4 model3  7.9 1  0.0155
#5 model4 13.1 1  0.0012
#6 model5 18.3 1  <0.001
#7 model6 23.5 1  <0.001
#8 model7 28.4 1  <0.001
#9 model8 33.3 1  <0.001
#10 model9 37.9 1  <0.001


# With N=2, Finding estimates for other parameters : kon.s, k.s, x0 

multiple.slip.func.expo5 = function(sheer, kon.s, k.s, x0, ts){
  sigabc = 2.39;
  radius = 3.5;
  N=2;
  sheer.foo = sheer ;  
  force.foo = ((28.4*pi*radius^3)/(sigabc))*sheer.foo; #in PN if sheer in Pa
  tau = (1/k.s) * (factorial(N-1)/N) * (kon.s/k.s)^(N-1) * exp( -beta * force.foo* x0* sum(1/1:N))
  -sum( dexp (ts, rate = 1/tau, log = TRUE ))
}

slip.mle.multi = mle2(multiple.slip.func.expo5, start=list (kon.s = 1.25, k.s = 0.31, x0 = 0.061),
                      data = list(sheer = multiple.slip.data$SheerStress, 
                                  ts = multiple.slip.data$Lifetime), 
                      method = "L-BFGS-B", lower = c( kon.s = 1, k.s = 0.2, x0 = 0.03), 
                      upper =  c(kon.s =2, k.s = 0.6, x0 = 0.08))

#
# Coefficients:
#   kon.s        k.s         x0 
# 1.21654640 0.31943957 0.05548954 
#
# Log-likelihood: -510.83 
#

# Confidence Interval
# With 3 free parameters kon.s, k.s, x0, the function 
# we do not get the output with confint(profile(fit))

confint(profile(slip.mle.multi))

#             2.5 %     97.5 %
# kon.s         NA         NA
# k.s   0.24281539         NA
# x0    0.03535238 0.07570132
# with an error saying of "flat profiles"


########
#### REASON FOR ERROR IN confint : #####

# When trying to manually estimate confidence interal
# by fixing N=2 
# giving k.s as data vector ks.vector
# varying kon.s and x0 are parameters
# we get ks.profile with SAME Neg Log Likelihood values

multiple.slip.func.expo25 = function(sheer, kon.s, k.s , x0, ts){
  sigabc = 2.39;
  radius = 3.5;
  N=2;
  sheer.foo = sheer ;  
  force.foo = ((28.4*pi*radius^3)/(sigabc))*sheer.foo; #in PN if sheer in Pa
  tau = (1/k.s) * (factorial(N-1)/N) * (kon.s/k.s)^(N-1) * exp( -beta * force.foo* x0* sum(1/1:N))
  -sum( dexp (ts, rate = 1/tau, log = TRUE ))
}

ksvec = seq(0.1, 0.6, by = 0.05)
ks.profile = matrix(ncol = 3, nrow = length(ksvec))
colnames(ks.profile) = c("kon", "x0", "NLL")
for (i in 1:length(ksvec)) {
slip.mle.multi = mle2(multiple.slip.func.expo25, start=list (kon.s = 1.25, x0 = 0.061),
                      data = list(sheer = multiple.slip.data$SheerStress, 
                                  ts = multiple.slip.data$Lifetime,  k.s = ksvec[i]  ), 
                      method = "L-BFGS-B", lower = c( kon.s = 0.00001, x0 = 0.000003))

ks.profile[i, ] = c(coef(slip.mle.multi), -logLik(slip.mle.multi))
}

# > ks.profile
#         kon         x0      NLL
#[1,] 0.1192701 0.05551224 510.8304
#[2,] 0.2682997 0.05550060 510.8304
#[3,] 0.4769064 0.05549259 510.8304
#[4,] 0.7451542 0.05549147 510.8304
# AND SO ON


#### Continuing Slip bond analysis #### 

# BUT we already know the limits of k.s and x0 from single bond model,
# AND we know how they they both monotonically affect the tau function
# SO we can find the limit of k.on using the two extremes of k.s and x0.

par(mfrow=c(1,2))

multiple.slip.func.expo6 = function(sheer, kon.s, ts){
  sigabc = 2.39;
  radius = 3.5;
  sheer.foo = sheer ;
  k.s = 0.09743211
  x0 = 0.01457108 ;
  force.foo = ((28.4*pi*radius^3)/(sigabc))*sheer.foo; #in PN if sheer in Pa
  tau = (1/k.s^2) * (1/2) * (kon.s) * exp( -beta * force.foo* x0* (3/2) )
  -sum( dexp (ts, rate = 1/tau, log = TRUE ))
}

slip.mle.multi = mle2(multiple.slip.func.expo6, start=list (kon.s = 1.25),
                      data = list(sheer = multiple.slip.data$SheerStress, 
                                  ts = multiple.slip.data$Lifetime), 
                      method = "L-BFGS-B", lower = c( kon.s = 0.0000001), 
                      upper =  c(kon.s =1000))

#confint(profile(slip.mle.multi),  std.err = c(0.4,0.2,0.04) ,try_harder=TRUE)

confint(profile(slip.mle.multi))
plot(profile(slip.mle.multi))

#       kon.s
#       2.5 %     97.5 % 
#  0.05365173 0.06811167 <<<----- MINIMUM KON.S


multiple.slip.func.expo6 = function(sheer, kon.s, ts){
  sigabc = 2.39;
  radius = 3.5;
  sheer.foo = sheer ;
  k.s = 0.09743211; # lowest extremum
  x0 = 0.35478147;
  force.foo = ((28.4*pi*radius^3)/(sigabc))*sheer.foo; #in PN if sheer in Pa
  tau = (1/k.s^2) * (1/2) * (kon.s) * exp( -beta * force.foo* x0* (3/2) )
  -sum( dexp (ts, rate = 1/tau, log = TRUE ))
}

slip.mle.multi = mle2(multiple.slip.func.expo6, start=list (kon.s = 1.25),
                      data = list(sheer = multiple.slip.data$SheerStress, 
                                  ts = multiple.slip.data$Lifetime), 
                      method = "L-BFGS-B", lower = c( kon.s = 0.0000001), 
                      upper =  c(kon.s =1000))


confint(profile(slip.mle.multi))
plot(profile(slip.mle.multi))

#     2.5 %   97.5 % 
#  42.04367 53.37498 <<<------ MAXIMUM KON.S


### Hence range of k.on = 0.05365173 - 53.37498 
# Note this might not correspond to 95% condidence interval. 


## PLOTTING ##
par(mfrow=c(1,2))
sheer = multiple.slip.data$SheerStress 
outt = tapply (multiple.slip.data$Lifetime,multiple.slip.data$SheerStress, mean)
xvec = unique(multiple.slip.data$SheerStress)*1000
plot(xvec,outt, xlim = c(0,50), ylim = c(0,7), xlab = "Sheer Stress (mPa)", ylab = "Lifetime (sec)")


sdev = tapply (multiple.slip.data$Lifetime,multiple.slip.data$SheerStress, sd)
sem = sdev/ sqrt(table(multiple.slip.data$SheerStress))
arrows(xvec, outt-sem, xvec, outt+sem, length=0.05, angle=90, code=3)
N.s.mle = 2
kon.s.mle = coef(mle.n.slip[[1]])  

sigabc = 2.39;
radius = 3.5;
k.s = 0.3095;
x0= 0.0613 ;  

sheer.vec = seq(0, 50, by = 0.1) 
force.vec = ((28.4*pi*radius^3)/(sigabc))*(sheer.vec/1000)

tau.s = (1/k.s)*(factorial(N.s.mle-1)/N.s.mle) * (kon.s.mle /k.s)^(N.s.mle-1) * exp( -beta * force.vec* x0* sum(1/1:N.s.mle))

lines(sheer.vec, tau.s, type = "l", lty = 1, col = 'Red')

plot(multiple.slip.data$SheerStress*1000, multiple.slip.data$Lifetime, xlim = c(0,50),xlab = "Sheer Stress (mPa)", ylab = "Lifetime (sec)")


N.s.mle = 2
kon.s.mle = coef(mle.n.slip[[1]])  

sigabc = 2.39;
radius = 3.5;
k.s = 0.3095;
x0= 0.0613 ;  

sheer.vec = seq(0, 50, by = 0.1) 
force.vec = ((28.4*pi*radius^3)/(sigabc))*(sheer.vec/1000)

tau.s = (1/k.s)*(factorial(N.s.mle-1)/N.s.mle) * (kon.s.mle /k.s)^(N.s.mle-1) * exp( -beta * force.vec* x0* sum(1/1:N.s.mle))

lines(sheer.vec, tau.s, type = "l", lty = 1, col = 'Red' )



#################### CATCH BOND #####

multiple.catch.data = read.table("multiple_catch_bena.txt", header = T)

#Find optimum N#
nll.n.catch = rep(0, 9)
mle.n.catch = c(new("mle2"),new("mle2"),new("mle2"),new("mle2"),new("mle2"),new("mle2"),new("mle2"),new("mle2"),new("mle2"))
for(i in 2:10){
  multiple.catch.func.expo = function(sheer, kon.c, ts ){
    sigabc = 15.97 ;
    radius = 3.5;
    k.c = 0.3192;
    psi= 0.0263 ;  
    F0 = 13.87
    N=i;
    sheer.foo = sheer
    force.foo = ((28.4*pi*radius^3)/(sigabc))*sheer.foo; #in PN if sheer in Pa
    tau = (1/k.c) * (factorial(N-1)/N) * (kon.c/k.c)^(N-1) * exp( -beta * ((force.foo* sum(1/1:N) - F0)^2)  *psi)
    -sum( dexp (ts, rate = 1/tau, log = TRUE ))
  }
  catch.mle.expo = mle2(multiple.catch.func.expo, start=list (kon.c = 0.3),
                        data = list(sheer = multiple.catch.data$SheerStress, 
                                    ts = multiple.catch.data$Lifetime) , method = "L-BFGS-B", lower = 0.000001)
  kon = coef(catch.mle.expo)  
  mleval = logLik(catch.mle.expo)
  print (i)
  print (kon)
  print (mleval)
  nll.n.catch[i-1] = -logLik(catch.mle.expo)
  mle.n.catch[i-1] = catch.mle.expo
}


AICtab(mle.n.catch, weights = TRUE)

#N dAIC df weight
#6 model5  0.0 1  0.5288 <<<<-------- 'OPTIMUM' N
#7 model6  0.9 1  0.3377 <<<<-------- dAIC < 1 => Also a valid value of N that nicely fits the data
#5 model4  3.3 1  0.1006
#8 model7  5.7 1  0.0304
#4 model3 11.0 1  0.0022
#9 model8 14.1 1  <0.001
#3 model2 22.7 1  <0.001
#10 model9 25.7 1  <0.001
#2 model1 36.7 1  <0.001


multiple.catch.func.expo = function(sheer, kon.c, k.c , psi, F0, ts ){
  sigabc = 15.97 ;
  radius = 3.5;
  N=6;
  sheer.foo = sheer
  force.foo = ((28.4*pi*radius^3)/(sigabc))*sheer.foo; #in PN if sheer in Pa
  tau = (1/k.c) * (factorial(N-1)/N) * (kon.c/k.c)^(N-1) * exp( -beta * ((force.foo* sum(1/1:N) - F0)^2)  *psi)
  -sum( dexp (ts, rate = 1/tau, log = TRUE ))
}


catch.mle.expo = mle2(multiple.catch.func.expo, start=list (kon.c = 0.3, k.c = 0.3192,  psi= 0.0263,   F0 = 13.87),
                      data = list(sheer = multiple.catch.data$SheerStress, 
                                  ts = multiple.catch.data$Lifetime) , method = "L-BFGS-B", lower = 0.000001)

#Coefficients:
#  kon.c         k.c         psi          F0 
#0.24155461  0.46355107  0.03005871 13.87074834 
#
#Log-likelihood: -320.3 


# PLOTTING THINGS  #

par(mfrow=c(1,2))
outt = tapply (multiple.catch.data$Lifetime,multiple.catch.data$SheerStress, mean)
xvec = unique(multiple.catch.data$SheerStress)*1000
plot(xvec,outt, xlim = c(0,50), ylim = c(0,3),  xlab = "Sheer Stress (mPa)", ylab = "Lifetime (sec)")


sdev = tapply (multiple.catch.data$Lifetime,multiple.catch.data$SheerStress, sd)
sem = sdev/ sqrt(table(multiple.catch.data$SheerStress))
arrows(xvec, outt-sem, xvec, outt+sem, length=0.05, angle=90, code=3)

N.mle = 7
kon.c.mle = coef(mle.n.catch[[7-1]])
sigabc = 15.97 ;
radius = 3.5;
k.c = 0.3192;
psi= 0.0263 ;  
F0 = 13.87

sheer.vec = seq(0, 50, by = 0.1) 
force.vec = ((28.4*pi*radius^3)/(sigabc))*(sheer.vec/1000)
tau = (1/k.c) * (factorial(N.mle-1)/N.mle) * (kon.c.mle/k.c)^(N.mle-1) * exp( -beta * ((force.vec* sum(1/1:N.mle) - F0)^2)  *psi)
lines( sheer.vec ,tau, type = "l", lty = 1, col = 'Blue' )


plot(multiple.catch.data$SheerStress*1000, multiple.catch.data$Lifetime, xlim = c(0,40),xlab = "Sheer Stress (mPa)", ylab = "Lifetime (sec)")
lines( sheer.vec ,tau, type = "l", lty = 1, col = 'Blue' )



#######THE###END#######

