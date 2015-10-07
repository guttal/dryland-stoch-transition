#Written by Vishu Guttal on 6th October 2015 (Tuesday)
#different models are fun

#Libraries
require('moments')
require('Kendall')
require('KernSmooth')

#One variable vegn model under grazing: May's model.
may <- function(v, r, k, v0, h)
{
  return(r*v*(1-v/k) - h*grazingt2(v,v0))
}

grazingt2 <- function(v,v0)
{
  return(v^2/(v0^2+v^2))
}

simulate_maystochastic <- function(r=1,k=10,v0=1,h=2,sigma_h=0.5,sigma_rate=0.1)
{

  tmax=1000
  dt=0.01
  
  v=numeric(tmax);
  v[1] = k;
  
  for (t in 2:tmax)
  {
    v[t] = v[t-1] + dt*may(v[t-1], r, k, v0, h) + sigma_h*sqrt(dt)*grazingt2(v[t-1],v0)*rnorm(1,0,1)   
  }
  
  return(v)
  
}

simulate_maystochastic_gradualrestoration <- function(r=1,k=10,v0=1,h0=3,sigma_h=0.5,sigma_rate=0.1)
{
  
  tmax=1000
  dt=0.01
  
  v=numeric(tmax);
  v[1] = 1;
  h[1] = h0;
  hrate=1;
  
  for (t in 2:tmax)
  {
    v[t] = v[t-1] + dt*may(v[t-1], r, k, v0, h[t-1]) + sigma_h*sqrt(dt)*grazingt2(v[t-1],v0)*rnorm(1,0,1)   
    h[t] = h[t-1] - dt*hrate
  }
  
  return(v)
  
}

simulate_maystochastic_stochasticrestoration <- function(r=1,k=10,v0=1,h0=3,sigma_h0=0.5,sigma_rate=0.1)
{
  
  tmax=1000
  dt=0.01
  
  v=numeric(tmax);
  v[1] = 1;
  sigma_h[1] = sigma_h0;
  sigma_hrate=1;
  
  for (t in 2:tmax)
  {
    v[t] = v[t-1] + dt*may(v[t-1], r, k, v0, h) + sigma_h[t-1]*sqrt(dt)*grazingt2(v[t-1],v0)*rnorm(1,0,1)   
    sigma_h[t] = sigma_h[t-1] - dt*sigma_hrate
  }
  
  return(v)
  
}

preprocess <- function(precrash_tsdata,bw)
{
  
  #This function takes precrash time series data (precrash_tsdata)
  #Then smoothens it using a Gaussian kernel. 
  #It returns a list of smoothened and detrended/residual data.
  
  smoothdata=ksmooth(1:N,precrash_tsdata,kernel=c("box","normal"),bandwidth=bw)
  
  smooth=smoothdata$y
  #Find residual data
  residuals=precrash-smooth
  return(list(ts_smooth=smooth,ts_residuals=residuals))
}

ews_kendall <- function(ts_residuals, l_rw, l_ktau) 
  {
  #A function to calculate early warning indicators and their respective Kendall-tau. 
  
  #Declare arrays   
  N=length(residuals)
  var_residuals=numeric(N)
  skew_residuals=numeric(N)
  acf_residuals=numeric(N)
  avgspec_residuals=numeric(N)
  
  #Apply rolling window to residuals.
  #Residual data analysis.
  for (i in 1:(N-l_rw+1))
  {
    rolldata = residuals[i:(i+l_rw-1)];
    
    var_residuals[i+l_rw-1] = var(rolldata);
    skew_residuals[i+l_rw-1] = skewness(rolldata);
    acf_residuals[i+l_rw-1] = acf(rolldata,plot=FALSE)$acf[2];
    
    spec_residuals = spectrum(rolldata,plot=FALSE)$spec;
    avgspec_residuals[i+l_rw-1] = mean(spec_residuals[2:floor(l_rw/8)]);
  }
  
  init_index=N-l_ktau 
  
  kendaltau_var = Kendall(1:l_ktau, var_residuals[init_index:(N-1)])$tau[1]
  kendaltau_acf = Kendall(1:l_ktau, acf_residuals[init_index:(N-1)])$tau[1]
  kendaltau_avgspec = Kendall(1:l_ktau, avgspec_residuals[init_index:(N-1)])$tau[1]
  kendaltau_skew = Kendall(1:l_ktau, skew_residuals[init_index:(N-1)])$tau[1]
  #Output results
  ews_n_kendall = list(acf_residuals=acf_residuals, var_residuals=var_residuals, spec_residuals=avgspec_residuals,
                       ktau_acf=kendaltau_acf,  ktau_var=kendaltau_var, ktau_spec=kendaltau_avgspec)
  return(ews_n_kendall)
}
