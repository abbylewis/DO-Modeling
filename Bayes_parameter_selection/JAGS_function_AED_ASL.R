# plug and play scripts 
# JAZ 2019-02-15
# WB Updates
# MEL updates for seasonal for-loop 30JUL19
# ASL adapted 3 Feb 2020

jags_plug_ins <- function(model_name){

#JAGS Plug-ins: Add each separate model here 
#variable.names are variables you would like to plot for model convergence (e.g., excludes mu)
#variable.names.out are all variables you would like to monitor in the jags run 
#init are a range of initial conditions for parameters in each of 3 chains 

#Seasonal_DO_model
  data.Seasonal_DO_model_AED <- list(y=y, year_no = year_no, beta.m_R20=log(.7),beta.m_theta=1.08,beta.m_ko2=log(.2), beta.m_sss_scalar=2.2, beta.v_R20=100, beta.v_theta=5000,beta.v_ko2=50,beta.v_sss_scalar=4, SSS = SSS, temp = temp, season_days=season_days,tau_ic = 100,m_proc = .1,sd_proc = .02, m_obs = .05, sd_obs = .02)
  variable.names.Seasonal_DO_model_AED <- c("tau_proc", "tau_obs","R20", "theta", "ko2","sss_scalar")
  variable.namesout.Seasonal_DO_model_AED <- c("tau_proc", "tau_obs","mu","R20", "theta", "ko2","sss_scalar")
  init.Seasonal_DO_model_AED <- list(list(tau_proc=0.001, tau_obs = 0.01, R20=.7, theta=1.08, ko2=.2,sss_scalar =2.2), list(tau_proc=0.01,  tau_obs = .1, R20=.7, theta=1.08, ko2=.2,sss_scalar =2.2), list(tau_proc=1, tau_obs = .1, R20=.7, theta=1.08, ko2=.2,sss_scalar =2.2))
  params.Seasonal_DO_model_AED <- c("tau_proc","R20", "theta","ko2","sss_scalar","tau_obs")
  

  data = eval(parse(text = paste0('data.', model_name)))
  variable.names = eval(parse(text = paste0('variable.names.', model_name)))
  variable.namesout = eval(parse(text = paste0('variable.namesout.', model_name)))
  init = eval(parse(text = paste0('init.', model_name)))
  params = eval(parse(text = paste0('params.', model_name)))

  return(list(data.model = data, variable.names.model = variable.names, variable.namesout.model = variable.namesout, init.model = init, params.model = params)) 
}


# converting mmol/m3 to mg/l o2
# 31.25 mmol/m3 (K from AED) * 10^-3 m3 / L * 31.998 g/mol = 1 mg/L