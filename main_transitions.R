# main function for code on distinguishing transitions

#default values for May's model for stochastic transitions
r=1
k=10
v0=1
h=2
sigma_h=3 #value of noise
sigma_rate=0.1

#run model to obtain one time series of transition
ts_may_gradrest = simulate_maystochastic_gradualrestoration(r,k,v0,h,sigma_h,sigma_rate)
ts_may_stochrest = simulate_maystochastic_stochasticrestoration(r,k,v0,h,sigma_h,sigma_rate)

#obtain ews
##Detrend data
#residual_gradrest=preprocess(ts_may_gradrest,bw=25)
#residual_stochrest=preprocess(ts_may_stochrest,bw=25)

##Calcualte ews

#plot results
par(mfrow=c(3,2))

##plot biomass data
plot(ts_may_gradrest,ylim=c(0,2*k),ylab="Biomass",xlab="time")
plot(ts_may_stochrest,ylim=c(0,2*k),col='red',ylab="Biomass",xlab="time")

##plot 



