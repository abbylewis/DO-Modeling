#Created by ASL 05Apr20
#Dealing with uncertinaty in EnKF functions

#' Sumarize variance for uncertainty analysis
#' 
#' @param sd_df full dataframe of variance
#' @param name type of uncertainty
#' @param year year these data are from
#' @return data.frame with summarized variance
#' 
sd_sum <- function(sd_df,name,year){ 
  out <- sd_df%>%
    filter(!is.na(TODAY))
  out <- out[,31:44]
  out <- out%>%
    summarize_all(~mean(.,na.rm=T))%>%
    mutate(uncert = paste(name),
           year = year)
  out
}



#' Calculate the average temperature, oxygen concentration, and hypolimnetic volume as inputs to the model
#' 
#' @param start vector of start dates for each year
#' @param stop vector of stop dates for each year
#' @param CTD CTD data file
#' @param SSS SSS driver 
#' @return average temperature, oxygen concentration, and hypolimnetic volume as a list
#' 
calc_avg <- function(start,stop,CTD,SSS){
  dates_init <- tibble(start,stop)
  dates <- dates_init %>%
    mutate(year = year(start))%>%
    dplyr::select(year,start,stop)
  
  allDays<-tibble(seq(min(start),max(stop),by = "days"))
  colnames(allDays) <- "Date"
  inputs = allDays %>%
    left_join(CTD %>%
                dplyr::select(Date,Conc,Temp,hypoVolume,thermo_depth,Chla_ugL,SA_m2), by = "Date")%>%
    left_join(SSS %>% 
                dplyr::select(time, scfm), by = c("Date"="time"))%>%
    arrange(Date)
  inputs = inputs%>%
    group_by(Date)%>%
    summarise(Conc = mean(Conc),
              Temp = mean(Temp),
              hypoVolume = mean(hypoVolume),
              thermo_depth = mean(thermo_depth),
              scfm = mean(scfm),
              Chla = mean(Chla_ugL),
              SA = mean(SA_m2))
  inputs_inRange <- inputs%>%
    mutate(year = year(Date))%>%
    left_join(dates)%>%
    filter(Date>=start,
           Date<=stop)%>%
    dplyr::select(-start,-stop,-year)
  
  toCalcAvg <- inputs_inRange%>%
    fill(Conc,Temp,hypoVolume,thermo_depth,Chla,SA)
  avg_temp <- mean(toCalcAvg$Temp)
  avg_o2 <- mean(toCalcAvg$Conc)
  avg_hypoVolume <- mean(toCalcAvg$hypoVolume)
  return(list(avg_temp,avg_o2,avg_hypoVolume))
}



#' Creates an input dataframe for this year
#' 
#' @param start start date for this year
#' @param stop stop date for this year
#' @param CTD CTD data file
#' @param SSS SSS driver file
#' @return dataframe of inputs for this year, including oxygen concentration, temperature, thermocline depth, oxygenation (scfm)
#' chlorophyll-a, and surface area
#' 
inputs_year<- function(start, stop, CTD, SSS){
  dates = data.frame(seq(start,stop,by = "days"))
  colnames(dates) <- "Date"
  inputs = dates %>%
    left_join(CTD %>%
                dplyr::select(Date,Conc,Temp,hypoVolume,thermo_depth,Chla_ugL,SA_m2), by = "Date")%>%
    left_join(SSS %>% 
                dplyr::select(time, scfm), by = c("Date"="time"))%>%
    arrange(Date)
  if(all(is.na(inputs$Chla_ugL))){
    inputs$Chla_ugL<-1
  }
  inputs = inputs%>%
    group_by(Date)%>%
    summarise(Conc = mean(Conc),
              Temp = mean(Temp),
              hypoVolume = mean(hypoVolume),
              thermo_depth = mean(thermo_depth),
              scfm = mean(scfm),
              Chla = mean(Chla_ugL),
              SA = mean(SA_m2))
  inputs
}



#' Epic function to run the EnKF at a given time point
#' 
#' @param inputs inputs file for this year 
#' @param obs all observations
#' @param today today's date
#' @param n_days forecast horizon
#' @param model_name model scenario ("normal","temp","o2","sss")
#' @param uncert type of uncertainty ("all", "driver","param","init","proc")
#' @return enKF results
#' 
run_do_hindcast <- function(inputs, obs, today, n_days = 14, model_name = "normal", uncert = "all", parms, start, stop, realDrivers){
  obs = obs%>%
    filter(datetime<=today)
  if(realDrivers == F){
    inputs[inputs$Date > today,c(2,3,4,5,7,8)] <- NA
  }
  inputs = inputs%>% #rule = 2 makes points outside the range of observed points equal to the last observation
    mutate(Conc = na.approx(Conc, rule = 2), ### O2 is only used for temp and SSS models (set to average) and when realDrivers == T
           Temp = na.approx(Temp, rule = 2), ### Recent values will be overwritten below
           hypoVolume = na.approx(hypoVolume, rule = 2),
           thermo_depth = na.approx(thermo_depth, rule = 2),
           scfm = ifelse(Date>=as.Date("2013-05-15"),scfm,0),
           SSS_add_conc = scfm*50*1000000/hypoVolume,
           Chla = na.approx(Chla, rule = 2),
           SA = na.approx(SA, rule = 2))
  inputs = inputs%>%
    filter(Date<=stop, 
           Date>=start,
           Date<=today+n_days)
  if(realDrivers == F){
  inputs$Temp[inputs$Date > max(obs$datetime)]<- inputs$Temp[inputs$Date == max(obs$datetime)]+temp_slope*seq(1,n_days)
  }
  simulation_time <- as.numeric(difftime(min(today+n_days, stop), start, unit = "days")+1) #days
  #Assemble driver data
  model_inputs <- list(datetime = inputs$Date,
                       SSS = inputs$SSS_add_conc,
                       temp = inputs$Temp,
                       O2_mgL = inputs$Conc) #included for temp and SSS models (o2 = average o2)
  model_inputs<- data.frame(model_inputs)
  
  #Set initial conditions
  #yini <- c(
  #  O2_mgL = inputs$Conc[inputs$Date == start]
  #)
  if(model_name %in% c("temp","SSS")|realDrivers == T){
    est_out = EnKF(n_en = n_en, 
                   start = start,
                   stop = min(stop, today+n_days),
                   time_step = "days",
                   obs_file = obs, 
                   driver_file = model_inputs, 
                   n_states_est = 1,
                   n_params_est = 4,
                   n_params_obs = 0,
                   n_drivers = 3,
                   parm_init = parms, 
                   obs_cv = obs_cv,
                   param_cv = param_cv,
                   driver_cv = driver_cv, 
                   proc_sd = proc_sd,
                   init_cond_cv = init_cond_cv,
                   model = no_O2_model,
                   today = today,
                   uncert = uncert)
  }else{
    #Run EnKF
    est_out = EnKF(n_en = n_en, 
                   start = start,
                   stop = min(stop, today+n_days),
                   time_step = "days",
                   obs_file = obs, 
                   driver_file = model_inputs, 
                   n_states_est = 1,
                   n_params_est = 4,
                   n_params_obs = 0,
                   n_drivers = 2,
                   parm_init = parms, 
                   obs_cv = obs_cv,
                   param_cv = param_cv,
                   driver_cv = driver_cv, 
                   proc_sd = proc_sd,
                   init_cond_cv = init_cond_cv,
                   model = O2_model,
                   today = today,
                   uncert = uncert)
  }
  return(est_out)
}



#' Function to plot the parameter values over time
#' 
#' @param est_out output of the EnKF
#' 
plot_hindcast <- function(est_out){
  par(mfrow = c(2,2))
  param_names = c("R20","theta","ko2","sss_scalar")
  for(i in seq(1,length(param_names))){
    plot_param(est_out,i)
  }
}



#' Function to plot oxygen concentrations over time
#' 
#' @param est_out output of the EnKF
#' @param today today's date
#' @param start start date for this year
#' @param stop stop date for this year
#' @param n_days forecast horizon
#' 
plot_o2 = function(est_out, today, start, stop, n_days = 14){
  mean_o2_est = apply(est_out$Y[1,,], 1, FUN = mean)
  plot(mean_o2_est ~ est_out$dates, type ='l', 
       ylim = c(-1,25),#max(c(est_out$Y[1,,], est_out$obs[1,,]),na.rm = T)),
       col = 'grey', ylab = 'Dissolved oxygen (mg L-1)', xlab = '',main = year(est_out$dates[1]), xlim = c(start, stop))
  for(i in 2:n_en){
    lines(est_out$Y[1,,i] ~ est_out$dates, 
          col = 'grey')
  }
  lines(mean_o2_est ~ est_out$dates, col = 'black', lwd =2 )
  future_obs <- obs[obs$datetime>today&obs$datetime<=(today+n_days),]
  obs <- obs[obs$datetime<=today,]
  points(obs$O2_mgL~obs$datetime, pch = 16, col = rgb(1,0,0))
  points(future_obs$O2_mgL~future_obs$datetime, pch = 16, col = rgb(0,0,1))
  abline(v = today)
  legend('topright',c('Model','Observations'),lty=c(1,NA),pch=c(NA,16))
  #arrows(est_out$dates, est_out$obs[1,,] - 
  #        est_out$state_sd, 
  #     est_out$dates, est_out$obs[1,,] +
  #      est_out$state_sd, 
  #   code = 3, length = 0.1, angle = 90, col = 'red')
}



#' Function to plot a given parameter over time
#' 
#' @param est_out results of kalman filter
#' @param num parameter number
#' @param name of parameter
#' 
plot_param = function(est_out,num,name = "Parameter value"){
  num = num+1
  mean_param_est = apply(est_out$Y[num,,], 1, FUN = mean)
  plot(mean_param_est ~ est_out$dates, type ='l', 
       ylim = range(est_out$Y[num,,]),
       col = 'grey', ylab = name, xlab ='', main = year(est_out$dates[1]))
  for(i in 2:n_en){
    lines(est_out$Y[num,,i] ~ est_out$dates, col = 'grey')
  }
  lines(mean_param_est ~ est_out$dates, col = 'black', lwd = 2)
}



#' Function to create an empty dataframe to hold predctions (forecasted), future observations, and 
#' variance of ensemples at each forecast horizon
#' 
#' @param start start date for this year
#' @param stop stop date for this year
#' @return empty results dataframe to hold predctions (forecasted), future observations, and 
#' variance of ensemples at each forecast horizon
#' 
createResultsDF <- function(start, stop){
  days <- as.numeric(difftime(stop,start, units = "days"))
  rows <- floor(days/run_space)+1
  cols <- n_days*3+1
  results <- matrix(NA, nrow = rows, ncol = cols)
  n <- rep(seq(1,n_days),3)
  type <- rep(c("_pred","_obs","_var"),each = n_days)
  col_names <- c("TODAY",paste("PLUS",n,type,sep = ""))
  colnames(results)<-col_names
  return(results)
}





#' Function to add all dates to the observation file (unsampled dates should be NA)
#' 
#' @param start start date for the summer
#' @param stop stop date for the summer
#' @param obs observation dataframe (ctd)
#' @return file with all dates this summer and oxygen observations where applicable
#' 
extendObsDF <- function(start,stop,obs){
  dates = data.frame(seq(start,stop,by = "days"))
  colnames(dates) <- "datetime"
  obs_allDates = dates %>%
    left_join(obs)
  return(obs_allDates)
}



#' Function to create a dataframe of RMSE, leaving out the first 14 days
#' 
#' @param n_days forecast horizon
#' @param results results from enkf
#' @param cali calibration period
#' @return RMSE dataframe
#' 
createRmseDF <- function(n_days,results, cali = 14){
  n <- seq(1,n_days)
  val <- rep(NA,n_days)
  class(val) <- "numeric"
  rmse_thisYear <- data.frame(n,val)
  results <- tail(results, n = -cali)
  results <- as.data.frame(results)
  for(i in seq(1:n_days)){
    predicted <- results[i+1]
    observed <- results[i+1+n_days]
    rmse_thisYear$val[rmse_thisYear$n == i] <- rmse(observed[!is.na(observed)], predicted[!is.na(observed)])
  }
  rmse_thisYear
}

#' Function to create a dataframe of NMAE, leaving out the first 14 days
#' 
#' @param n_days forecast horizon
#' @param results results from enkf
#' @param cali calibration period
#' @return RMSE dataframe
#' 
createNmaeDF <- function(n_days,results, cali = 14){
  n <- seq(1,n_days)
  val <- rep(NA,n_days)
  class(val) <- "numeric"
  nmae_thisYear <- data.frame(n,val)
  results <- tail(results, n = -cali)
  results <- as.data.frame(results)
  for(i in seq(1:n_days)){
    predicted <- results[i+1]
    observed <- results[i+1+n_days]
    nmae_thisYear$val[nmae_thisYear$n == i] <- compute.nmae(observed[!is.na(observed)], predicted[!is.na(observed)])
  }
  nmae_thisYear
}

#' Function to calculate one rmse value, leaving out the first 14 days
#' 
#' @param n_days forecast horizon
#' @param results results from enkf
#' @param cali calibration period
#' @return RMSE dataframe
#' 
calcRmseTotal <- function(n_days,results, cali = 14){
  #n <- seq(1,n_days)
  #val <- rep(NA,n_days)
  #class(val) <- "numeric"
  #rmse_thisYear <- data.frame(n,val)
  results <- tail(results, n = -cali)
  results <- as.data.frame(results)
  predicted = c()
  observed = c()
  for(i in seq(1:n_days)){
    predicted <- c(predicted, results[i+1])
    observed <- c(observed, results[i+1+n_days])
  }
  predicted = unlist(predicted)
  observed = unlist(observed)
  rmse(observed[!is.na(observed)], predicted[!is.na(observed)])
}


#' Take a results file and calculate RMSE for each forecast horizon starting from an observation
#' 
#' @param n_days number of days forecasted
#' @param results results file from EnKF
#' @param obs
#' @reutrn data.frame with RMSE for each forecast horizon this year
#' 
createRmseDF_filt <- function(n_days, results, obs, cali = 14){
  results <- as.data.frame(results)
  obs_dates <- as.numeric(obs$datetime)
  results_filtered <- results[results$TODAY %in% obs_dates,]
  rmse_filt <- createRmseDF(n_days, results_filtered, cali)
  return(rmse_filt)
}


#' Epic function to run everything!
#' 
#' @param start start date for the summer
#' @param stop stop date for the summer 
#' @param n_days forecast horizon
#' @param run_space days between runs (1 = run every day)
#' @param obs observation file
#' @param gif should images be stored and create a gif?
#' @param archiveForecasts should forecasts be archived?
#' @param remove should images be removed after making a gif?
#' @param delay how much time between images in the gif
#' @param model_name which model? ("normal","temp","o2","sss")
#' @param avg_o2 calculated average amount of oxygen
#' @param avg_temp calculated average temperature
#' @param uncert what type of uncertainty? ("all", "driver","param","init")
#' @return results dataframe with predctions (forecasted), future observations, and 
#' variance of ensemples at each forecast horizon
#' 
runForecasts <- function(start, stop, n_days, run_space, obs, gif = TRUE, archiveForecasts = FALSE, remove = FALSE, delay = 30, model_name = "full",avg_o2,avg_temp,uncert,parms, realDrivers){
  #Model types = full (""), temp ("_temp"), o2 ("_o2")
  if(model_name == "full"){model = ""}
  if(model_name == "temp"){model = "_temp"}
  if(model_name == "o2"){model = "_o2"}
  if(model_name %in% c("SSS","sss")){model = "_sss"}
  
  #Create results dataframe
  results<- createResultsDF(start, stop)
  #Create obs file with all dates (filling dates for NAs)
  obs_allDates <- extendObsDF(start,stop,obs)
  #Create inputs (drivers) for this year
  inputs_thisYear <- inputs_year(start,stop, CTD, SSS)
  if(model_name == "temp"|model_name == "SSS"){inputs_thisYear$Conc <- avg_o2}  ### Need to adjust this for param est
  if(model_name == "o2"|model_name == "SSS"){inputs_thisYear$Temp <- avg_temp}
  today <- start
  if(gif == TRUE){
    dateRun <- format(Sys.Date(),"%d%b%y")
    year <- year(start)
    dir1 <- paste("../DO_modeling_figures",dateRun,sep = "/")
    dir2<-paste(dir1,"/",year,model,sep = "")
    dir.create(dir1)
    dir.create(dir2)
    #setwd(dir2)
    png(file=paste(dir2,"/forecast%03d.png",sep=""), width=480, height=480)
  }
  while(today < stop){
    row<- as.numeric(difftime(today,start))/run_space+1
    cols <- n_days*2+1
    est_thisYear <- run_do_hindcast(inputs_thisYear, obs, today, n_days, model_name = model_name,uncert = uncert, parms = parms, start = start, stop = stop, realDrivers = realDrivers)
    mean_o2_est <- apply(est_thisYear$Y[1,,], 1, FUN = mean)
    mean_o2_est_short <- mean_o2_est[(run_space*(row-1)+2):min(length(mean_o2_est),(run_space*(row-1)+1+n_days))]
    var_o2_est <- apply(est_thisYear$Y[1,,], 1, FUN = var)
    var_o2_est_short = var_o2_est[(run_space*(row-1)+2):min(length(var_o2_est),(run_space*(row-1)+1+n_days))]
    obs_toAdd <- obs_allDates$O2_mgL[obs_allDates$datetime>today & obs_allDates$datetime<=today+n_days]
    if(length(obs_toAdd)<n_days){
      obs_toAdd<-c(obs_toAdd, rep(NA,n_days-length(obs_toAdd)))
    }
    results[row,(2+n_days):(1+2*n_days)]<- obs_toAdd
    results[row,2:(1+n_days)]<-c(mean_o2_est_short,rep(NA,(n_days-length(mean_o2_est_short))))
    results[row,(2+2*n_days):(1+3*n_days)]<-c(var_o2_est_short,rep(NA,(n_days-length(var_o2_est_short))))
    results[row,1]<-today
    if(gif == TRUE){
      plot_o2(est_thisYear, today, start, stop)
    }
    if(archiveForecasts == TRUE){
      dateRun <- format(Sys.Date(),"%d%b%y")
      year <- year(start)
      dir1.1 <- paste("../Archived_forecasts",dateRun,sep = "/")
      dir2.1<-paste(dir1.1,"/",year,model,sep = "")
      dir.create(dir1.1)
      dir.create(dir2.1)
      write.csv(est_thisYear$Y[1,,], paste(dir2.1,"/",format(today,"%d%b%y"),".csv",sep = ""))
      dir1.1 <- paste("../Archived_temps",dateRun,sep = "/")
      dir2.1<-paste(dir1.1,"/",year,model,sep = "")
      dir.create(dir1.1)
      dir.create(dir2.1)
      write.csv(est_thisYear$drivers[,2,], paste(dir2.1,"/",format(today,"%d%b%y"),".csv",sep = ""))
    }
    today <- today+run_space
  }
  if (gif == TRUE){
    convert <- paste("convert -delay ",delay," ",dir2,"/forecast*.png ",dir2,"/animated_forecast.gif", sep = "")
    system(convert)
    dev.off()
    jpeg(paste(dir2,"/params.jpeg",sep = ""), width = 6, height = 5, units = "in", res= 300)
    plot_hindcast(est_thisYear)
    dev.off()
    # Remove png files
    if(remove == TRUE){
      file.remove(list.files(pattern=".png"))
    }
  }
  saved_R20 <-   as.data.frame(est_thisYear$Y[2,,], col.names = seq(1,n_en))
  saved_theta <- as.data.frame(est_thisYear$Y[3,,], col.names = seq(1,n_en))
  saved_ko2 <-   as.data.frame(est_thisYear$Y[4,,], col.names = seq(1,n_en))
  saved_sss <-   as.data.frame(est_thisYear$Y[5,,], col.names = seq(1,n_en))
  if(uncert == "all"){
    collumn_names = colnames(saved_R20 %>% mutate(Var = "R20",Date = row.names(saved_R20)))
    saved_params <- saved_R20 %>% mutate(Var = "R20",Date = row.names(saved_R20)) %>%
      full_join(saved_theta   %>% mutate(Var = "theta",Date = row.names(saved_R20)), by = collumn_names)%>%
      full_join(saved_ko2     %>% mutate(Var = "ko2",Date = row.names(saved_R20)), by = collumn_names)%>%
      full_join(saved_sss     %>% mutate(Var = "sss",Date = row.names(saved_R20)), by = collumn_names)%>%
      pivot_longer(1:100, "Sim")
    date <- format(Sys.Date(),"%d%b%y")
    dir.create(paste("../DO_modeling_results/",date,sep = ""))
    year <- year(start)
    write.csv(saved_params,paste("../DO_modeling_results/",date,"/params_",year,"_",model_name,".csv", sep = ""))
  }
  return(results)
}



#' Function to run a persistence forecast 
#' 
#' @param start start date
#' @param stop stop date
#' @param n_days forecast horizon
#' @param run_space how often to run forecast (days; 1 = run every day)
#' @param obs observations data file
#' @return results dataframe with persistence forecasts
#' 
persistenceForecast <- function(start,stop,n_days,run_space,obs){
  #Create results dataframe
  results<- createResultsDF(start, stop)
  #Create obs file with all dates (filling dates for NAs)
  obs_allDates <- extendObsDF(start,stop,obs)
  today <- start
  while(today < stop){
    row<- as.numeric(difftime(today,start))/run_space+1
    cols <- n_days*2+1
    current_obs <- obs[obs$datetime<=today,]
    most_recent_obs <- current_obs$O2_mgL[which.max(current_obs$datetime)]
    forecast <- rep(most_recent_obs,n_days)
    obs_toAdd <- obs_allDates$O2_mgL[obs_allDates$datetime>today & obs_allDates$datetime<=today+n_days]
    if(length(obs_toAdd)<n_days){
      obs_toAdd<-c(obs_toAdd, rep(NA,n_days-length(obs_toAdd)))
    }
    results[row,(2+n_days):cols]<- obs_toAdd
    results[row,2:(1+n_days)]<-forecast
    results[row,1]<-today
    today <- today+run_space
  }
  results 
}
