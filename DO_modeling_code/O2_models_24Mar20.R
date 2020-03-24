
O2_model <- function(times, states, parms, inputs){
  
  #STATES 
  O2_mgL <- states[1]
  
  #PARMS 
  R20 <- as.numeric(parms[1])
  theta <- as.numeric(parms[2])
  ko2 <- as.numeric(parms[3])
  sss_scalar <- as.numeric(parms[4])
  
  #INPUT
  SSS <- inputs[1]
  temp <- inputs[2]
  
  #CALCULATIONS
  R = R20*theta^(temp-20)*O2_mgL/(O2_mgL+ko2)
  scalar = 1
  if(SSS>0){
    scalar = sss_scalar
  }
  OD = R*scalar
  
  #DERIVATIVE
  dO2_mgLdt <- SSS-OD
  
  #RETURN A LIST OF THE DERIVATIVES
  return(list(c(dO2_mgLdt), c(parms)))
}

no_O2_model <- function(times, states, parms, inputs){
  
  #STATES 
  state_O2 <- states[1]
  
  #PARMS 
  R20 <- as.numeric(parms[1])
  theta <- as.numeric(parms[2])
  ko2 <- as.numeric(parms[3])
  sss_scalar <- as.numeric(parms[4])
  
  #INPUT
  SSS <- inputs[1]
  temp <- inputs[2]
  O2_mgL <- inputs[3]
  
  #CALCULATIONS
  R = R20*theta^(temp-20)*O2_mgL/(O2_mgL+ko2)
  scalar = 1
  if(SSS>0){
    scalar = sss_scalar
  }
  OD = R*scalar
  
  #DERIVATIVE
  dO2_mgLdt <- SSS-OD
  
  #RETURN A LIST OF THE DERIVATIVES
  return(list(c(dO2_mgLdt), c(parms)))
}