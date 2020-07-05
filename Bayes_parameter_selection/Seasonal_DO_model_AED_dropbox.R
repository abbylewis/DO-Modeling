model{
  
  for(k in 1:max(year_no)){
    
    for(j in 2:max(season_days)){
      #this fits the latent DO to your observed DO
      y[k,j] ~ dnorm(mu[k,j],tau_obs)
      #process model for DO
      mu[k,j]~dnorm(lambda[k,j],tau_proc) 
      lambda[k,j] <- mu[k,j-1]+SSS[k,j-1]-R20*theta^(temp[k,j-1]-20)*mu[k,j-1]/(mu[k,j-1]+ko2)*ifelse(SSS[k,j]>0, sss_scalar,1)
    }
    
    
    #Loops through items in seasonal for-loop and defines initial conditions
    mu[k,1] ~ dnorm(y[k,1],tau_ic)
    #mu[k,1] ~ dnorm(x_ic,tau_ic) #keep in mind you'll need to index like a matrix 
    #mu_S[k,1]~dnorm(x_S_ic,tau_S_ic) 
    
  }
  #### Priors
  tau_proc ~ dgamma(sh_proc,ra_proc)
  R20 ~ dlnorm(beta.m_R20,beta.v_R20) 
  theta ~ dnorm(beta.m_theta,beta.v_theta) 
  ko2 ~ dlnorm(beta.m_ko2,beta.v_ko2) 
  sss_scalar ~ dnorm(beta.m_sss_scalar,beta.v_sss_scalar)
  tau_obs ~ dgamma(sh_obs,ra_obs)
  #tau_S_obs ~ dgamma(0.01, 0.01) 
  #tau_S_proc ~ dgamma(0.24, 300)
  #tau_S_proc ~ dnorm(0,.01)
  
  
  #calculations for gamma parameterized by mean (m) and standard deviation (sd)
  sh_proc <- pow(m_proc,2) / pow(sd_proc,2)
  ra_proc <- m_proc / pow(sd_proc,2)
  
  sh_obs <- pow(m_obs,2) / pow(sd_obs,2)
  ra_obs <- m_obs / pow(sd_obs,2)
  
}