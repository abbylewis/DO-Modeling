# plug and play scripts 
# JAZ 2019-02-15
# WB Updates
# MEL updates for seasonal for-loop 30JUL19
# ASL adapted 3 Feb 2020
#2020-06-11: Trying to not have o2 go above 100 mg/l.....

jags_plug_ins <- function(model_name){

#JAGS Plug-ins: Add each separate model here 
#variable.names are variables you would like to plot for model convergence (e.g., excludes mu)
#variable.names.out are all variables you would like to monitor in the jags run 
#init are a range of initial conditions for parameters in each of 3 chains 

#Seasonal_DO_model
  data.Seasonal_DO_model_noO2 <- list(y=y, year_no = year_no, beta.m_R20=.7,beta.m_theta=1.08,beta.m_ko2=log(.2), beta.m_sss_scalar=2.2, beta.v_R20=10, beta.v_theta=1000,beta.v_ko2=10,beta.v_sss_scalar=4, SSS = SSS, temp = temp, o2 = o2, season_days=season_days,tau_ic = 100,k_proc = .2,t_proc = 1, k_obs = .1, t_obs = 1)
  variable.names.Seasonal_DO_model_noO2 <- c("tau_proc", "tau_obs","R20", "theta", "ko2","sss_scalar")
  variable.namesout.Seasonal_DO_model_noO2 <- c("tau_proc", "tau_obs","mu","R20", "theta", "ko2","sss_scalar")
  init.Seasonal_DO_model_noO2 <- list(list(R20=.7, sss_scalar =2.2), list(R20=.7,sss_scalar =2.2), list(R20=.7, sss_scalar =2.2))
  params.Seasonal_DO_model_noO2 <- c("tau_proc","R20", "theta","ko2","sss_scalar","tau_obs")
  

  data = eval(parse(text = paste0('data.', model_name)))
  variable.names = eval(parse(text = paste0('variable.names.', model_name)))
  variable.namesout = eval(parse(text = paste0('variable.namesout.', model_name)))
  init = eval(parse(text = paste0('init.', model_name)))
  params = eval(parse(text = paste0('params.', model_name)))

  return(list(data.model = data, variable.names.model = variable.names, variable.namesout.model = variable.namesout, init.model = init, params.model = params)) 
}


# converting mmol/m3 to mg/l o2
# 31.25 mmol/m3 (K from AED) * 10^-3 m3 / L * 31.998 g/mol = 1 mg/L