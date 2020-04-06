# DO-Modeling

Models were created by ASL September 2019 - March 2020
Forecasting system developed by ASL March 2020 - April 2020

##Abstract
Warming temperatures and increased nutrient loads have led to decreased oxygen concentrations in lakes around the world. This presents a key management challenge for drinking water reservoirs, as decreased oxygen concentrations near the sediments can lead to persistent water quality problems. A thorough understanding of the factors governing bottom-water oxygen consumption would allow managers to anticipate and preempt water quality concerns, helping provide safe drinking water to people around the world. In this project, I am working to develop an operational system that can forecast oxygen concentrations in the hypolimnion of a drinking water reservoir using biweekly (two times per week) temperature and oxygen measurements. I will then use the ensemble Kalman filter to assimilate data and I will determine how the dominant sources of uncertainty change over time.

##Personnel
Abigail Lewis, Virginia Tech, aslewis@vt.edu

##Principal Investigator
Cayelan Carey, Virginia Tech, cayelan@vt.edu

##Github file structure

*DO_modeling_code* folder contains all files required to run the forecast. The most recent version is DO_forecasting_05Apr20.Rmd. 

*DO_modeling_data* folder contains the data used to run the model (automatically referenced in model code)

*DO_modeling_figures* folder contains figures and gifs from forecasts

*DO_modeling_results* contains forecasting results (RMSE, as well as mean and variation for each forecast run)

Forecasts are archived in the *Archived_forecasts* folder