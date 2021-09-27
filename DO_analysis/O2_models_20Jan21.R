#Created by ASL 20 Jan 21 to make it so temp and o2 models don't use the mean 
#they just don't include the other variable at all

O2_model <- function(times, states, parms, inputs, inc_metals){
  #Used for full model
  #STATES 
  O2_mgL <- states[1]
  
  #PARMS 
  R20 <- as.numeric(parms[1])
  theta <- as.numeric(parms[2])
  ko2 <- as.numeric(parms[3])
  sss_scalar <- as.numeric(parms[4])
  #fe_use <- as.numeric(parms[5])
  
  #INPUT
  SSS <- inputs[1]
  temp <- inputs[2]
  #Fe <- inputs[4]
  
  #CALCULATIONS
  BOD = R20*theta^(temp-20)*O2_mgL/(O2_mgL+ko2)
  if(inc_metals){
    COD = fe_use*Fe*O2_mgL/(O2_mgL+ko2)
  }else{
    COD = 0
  }
  R = BOD+COD
  scalar = 1
  if(SSS>0){
    scalar = sss_scalar
  }
  OD = R*scalar
  
  #DERIVATIVE
  dO2_mgLdt <- SSS-OD
  
  if(O2_mgL+dO2_mgLdt < 0){dO2_mgLdt = -O2_mgL}
  
  #RETURN A LIST OF THE DERIVATIVES
  return(list(c(dO2_mgLdt), c(parms)))
}

no_O2_model <- function(times, states, parms, inputs, inc_metals){
  #Used for realDrivers = T
  
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
  O2_mgL <- inputs[3] #average o2
  
  #CALCULATIONS
  R = R20*theta^(temp-20)*O2_mgL/(O2_mgL+ko2)
  scalar = 1
  if(SSS>0){
    scalar = sss_scalar
  }
  OD = R*scalar
  
  #DERIVATIVE
  dO2_mgLdt <- SSS-OD
  
  if(O2_mgL+dO2_mgLdt < 0){dO2_mgLdt = -O2_mgL}
  
  #RETURN A LIST OF THE DERIVATIVES
  return(list(c(dO2_mgLdt), c(parms)))
}

temp_model <- function(times, states, parms, inputs, inc_metals){
  #Used for temp-only model
  
  #STATES 
  O2_mgL <- states[1]
  
  #PARMS 
  R20 <- as.numeric(parms[1])
  theta <- as.numeric(parms[2])
  #ko2 <- as.numeric(parms[3])
  sss_scalar <- as.numeric(parms[3])
  
  #INPUT
  SSS <- inputs[1]
  temp <- inputs[2]
  
  #CALCULATIONS
  R = R20*theta^(temp-20)
  scalar = 1
  if(SSS>0){
    scalar = sss_scalar
  }
  OD = R*scalar
  
  #DERIVATIVE
  dO2_mgLdt <- SSS-OD
  
  if(O2_mgL+dO2_mgLdt < 0){dO2_mgLdt = -O2_mgL}
  
  #RETURN A LIST OF THE DERIVATIVES
  return(list(c(dO2_mgLdt), c(parms)))
}

no_temp_model <- function(times, states, parms, inputs, inc_metals){
  #Used for O2-only model
  
  #STATES 
  O2_mgL <- states[1]
  
  #PARMS 
  R20 <- as.numeric(parms[1])
  #theta <- as.numeric(parms[2])
  ko2 <- as.numeric(parms[2])
  sss_scalar <- as.numeric(parms[3])
  
  #INPUT
  SSS <- inputs[1]
  #temp <- inputs[2]
  
  #CALCULATIONS
  R = R20*O2_mgL/(O2_mgL+ko2)
  scalar = 1
  if(SSS>0){
    scalar = sss_scalar
  }
  OD = R*scalar
  
  #DERIVATIVE
  dO2_mgLdt <- SSS-OD
  
  if(O2_mgL+dO2_mgLdt < 0){dO2_mgLdt = -O2_mgL}
  
  #RETURN A LIST OF THE DERIVATIVES
  return(list(c(dO2_mgLdt), c(parms)))
}


sss_model <- function(times, states, parms, inputs, inc_metals){
  #Used for sss-only (null) model
  
  #STATES 
  O2_mgL <- states[1]
  
  #PARMS 
  R20 <- as.numeric(parms[1])
  sss_scalar <- as.numeric(parms[2])
  
  #INPUT
  SSS <- inputs[1]
  
  #CALCULATIONS
  R = R20
  scalar = 1
  if(SSS>0){
    scalar = sss_scalar
  }
  OD = R*scalar
  
  #DERIVATIVE
  dO2_mgLdt <- SSS-OD
  
  if(O2_mgL+dO2_mgLdt < 0){dO2_mgLdt = -O2_mgL}
  
  #RETURN A LIST OF THE DERIVATIVES
  return(list(c(dO2_mgLdt), c(parms)))
}
