calc_avg <- function(start,stop,CTD,SSS){
  dates_init <- data_frame(start,stop)
  dates <- dates_init %>%
    mutate(year = year(start))%>%
    select(year,start,stop)
  
  allDays<-data.frame(seq(min(start),max(stop),by = "days"))
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
    select(-start,-stop,-year)
  
  toCalcAvg <- inputs_inRange%>%
    fill(Conc,Temp,hypoVolume,thermo_depth,Chla,SA)
  avg_temp <- mean(toCalcAvg$Temp)
  avg_o2 <- mean(toCalcAvg$Conc)
  avg_hypoVolume <- mean(toCalcAvg$hypoVolume)
  return(list(avg_temp,avg_o2,avg_hypoVolume))
}

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

run_do_hindcast <- function(inputs, obs, today, n_days = 14, model_name = "normal"){
  obs = obs%>%
    filter(datetime<=today)
  inputs = inputs%>%
    filter(Date<=stop, 
           Date>=start,
           Date<=today+n_days)
  inputs[inputs$Date > today,c(2,3,4,5,7,8)] <- NA
  inputs = inputs%>%
    mutate(Conc = na.approx(Conc, rule = 2),
           Temp = na.approx(Temp, rule = 2),
           hypoVolume = na.approx(hypoVolume, rule = 2),
           thermo_depth = na.approx(thermo_depth, rule = 2),
           scfm = ifelse(Date>=as.Date("2013-05-15"),scfm,0),
           SSS_add_conc = scfm*50*1000000/hypoVolume,
           Chla = na.approx(Chla, rule = 2),
           SA = na.approx(SA, rule = 2))
  simulation_time <- as.numeric(difftime(min(today+n_days, stop), start, unit = "days")+1) #days
  #Assemble driver data
  model_inputs <- list(datetime = inputs$Date,
                       SSS = inputs$SSS_add_conc,
                       temp = inputs$Temp,
                       O2_mgL = inputs$Conc)
  model_inputs<- data.frame(model_inputs)
  
  #Set initial conditions
  yini <- c(
    O2_mgL = inputs$Conc[inputs$Date == start]
  )
  if(model_name %in% c("temp","SSS")){
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
                   init_cond_cv = init_cond_cv,
                   model = no_O2_model)
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
                   init_cond_cv = init_cond_cv,
                   model = O2_model)
  }
  return(est_out)
}

plot_hindcast <- function(est_out){
  par(mfrow = c(2,3))
  param_names = c("R20","theta","ko2","sss_scalar","chla_vmax")
  plot_param = function(est_out,num,name = "Parameter value"){
    num = num+1
    mean_param_est = apply(est_out$Y[num,,], 1, FUN = mean)
    plot(mean_param_est ~ est_out$dates, type ='l', 
         ylim = range(est_out$Y[num,,]),
         col = 'grey', ylab = name, xlab ='', main = param_names[num-1])
    for(i in 2:n_en){
      lines(est_out$Y[num,,i] ~ est_out$dates, col = 'grey')
    }
    lines(mean_param_est ~ est_out$dates, col = 'black', lwd = 2)
  }
  for(i in seq(1,length(param_names))){
    plot_param(est_out,i)
  }
}


plot_o2 = function(est_out, today, start, stop, n_days = 14){
  mean_o2_est = apply(est_out$Y[1,,], 1, FUN = mean)
  plot(mean_o2_est ~ est_out$dates, type ='l', 
       ylim = c(-1,25),#max(c(est_out$Y[1,,], est_out$obs[1,,]),na.rm = T)),
       col = 'grey', ylab = 'O2 (mg L-1)', xlab = '',main = year(est_out$dates[1]), xlim = c(start, stop))
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
  #arrows(est_out$dates, est_out$obs[1,,] - 
  #        est_out$state_sd, 
  #     est_out$dates, est_out$obs[1,,] +
  #      est_out$state_sd, 
  #   code = 3, length = 0.1, angle = 90, col = 'red')
}

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

createResultsDF <- function(start, stop){
  days <- as.numeric(difftime(stop,start, units = "days"))
  rows <- floor(days/run_space)+1
  cols <- n_days*2+1
  results <- matrix(NA, nrow = rows, ncol = cols)
  n <- rep(seq(1,n_days))
  type <- rep(c("_pred","_obs"),each = n_days)
  col_names <- c("TODAY",paste("PLUS",n,type,sep = ""))
  colnames(results)<-col_names
  return(results)
}

extendObsDF <- function(start,stop,obs){
  dates = data.frame(seq(start,stop,by = "days"))
  colnames(dates) <- "datetime"
  obs_allDates = dates %>%
    left_join(obs)
  return(obs_allDates)
}

createRmseDF <- function(n_days,results){
  n <- seq(1,n_days)
  val <- rep(NA,n_days)
  class(val) <- "numeric"
  rmse_thisYear <- data.frame(n,val)
  results <- as.data.frame(results)
  for(i in seq(1:n_days)){
    predicted <- results[i+1]
    observed <- results[i+1+n_days]
    rmse_thisYear$val[rmse_thisYear$n == i] <- rmse(observed[!is.na(observed)], predicted[!is.na(observed)])
  }
  rmse_thisYear
}

runForecasts <- function(start, stop, n_days, run_space, obs, gif = TRUE, archiveForecasts = FALSE, remove = FALSE, delay = 30, model_name = "full",avg_o2,avg_temp){
  #Model types = full (""), temp ("_temp"), o2 ("_o2")
  if(model_name == "full"){model = ""}
  if(model_name == "temp"){model = "_temp"}
  if(model_name == "o2"){model = "_o2"}
  if(model_name == "SSS"){model = "_sss"}
  
  #Create results dataframe
  results<- createResultsDF(start, stop)
  #Create obs file with all dates (filling dates for NAs)
  obs_allDates <- extendObsDF(start,stop,obs)
  #Create inputs (drivers) for this year
  inputs_thisYear <- inputs_year(start,stop, CTD, SSS)
  if(model_name == "temp"|model_name == "SSS"){inputs_thisYear$Conc <- avg_o2}
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
    est_thisYear <- run_do_hindcast(inputs_thisYear, obs, today, n_days, model_name = model_name)
    mean_o2_est <- apply(est_thisYear$Y[1,,], 1, FUN = mean)
    obs_toAdd <- obs_allDates$O2_mgL[obs_allDates$datetime>today & obs_allDates$datetime<=today+n_days]
    if(length(obs_toAdd)<n_days){
      obs_toAdd<-c(obs_toAdd, rep(NA,n_days-length(obs_toAdd)))
    }
    results[row,(2+n_days):cols]<- obs_toAdd
    results[row,2:(1+n_days)]<-tail(mean_o2_est, n = n_days)
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
      write.csv(est_thisYear$Y[1,,],paste(dir2,"/",format(today,"%d%b%y"),".csv",sep = ""))
    }
    today <- today+run_space
  }
  if (gif == TRUE){
    convert <- paste("convert -delay ",delay," ",dir2,"/forecast*.png ",dir2,"/animated_forecast.gif", sep = "")
    system(convert)
    dev.off()
    # Remove png files
    if(remove == TRUE){
      file.remove(list.files(pattern=".png"))
    }
  }
  return(results)
}


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