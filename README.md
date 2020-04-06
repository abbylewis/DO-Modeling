# DO-Modeling

Models were created by ASL September 2019 - March 2020
Forecasting system developed by ASL March 2020 - April 2020

## Abstract

Warming temperatures and increased nutrient loads have led to decreased oxygen concentrations in lakes around the world. This presents a key management challenge for drinking water reservoirs, as decreased oxygen concentrations near the sediments can lead to persistent water quality problems. A thorough understanding of the factors governing bottom-water oxygen consumption would allow managers to anticipate and preempt water quality concerns, helping provide safe drinking water to people around the world. While a variety of factors are known to influence the rate of oxygen demand (e.g., organic matter quantity and quality), relatively little work has studied the influence of oxygen concentrations or the interaction between oxygen concentrations and temperature. Furthermore, few studies have investigated these processes on the whole-ecosystem scale, which is most relevant to ecosystem management. In this study, we used seven years (2013-2019) of whole-ecosystem oxygenation experiments to determine the separate impacts of temperature and oxygen concentration on the rate of oxygen consumption in a small, eutrophic, drinking water reservoir. Whole-ecosystem oxygenation was performed using an engineered system which can add oxygen to the bottom-waters of the reservoir. We found that 98% of the variation in oxygen demand was explained by temperature and the rate of hypolimnetic oxygenation, with increasing temperature and oxygenation both increasing oxygen demand. Based on the data from our oxygenation experiments, I am working to develop an operational system that can forecast oxygen concentrations in the hypolimnion of the reservoir using biweekly (two times per week) temperature and oxygen measurements. I will then use the ensemble Kalman filter for data assimilation and determine how the dominant sources of uncertainty change over time.

## Personnel

Abigail Lewis, Virginia Tech, aslewis@vt.edu

## Principal Investigator

Cayelan Carey, Virginia Tech, cayelan@vt.edu

## Github file structure

*DO_modeling_code* folder contains all files required to run the forecast. The most recent version is DO_forecasting_05Apr20.Rmd. 

*DO_modeling_data* folder contains the data used to run the model (automatically referenced in model code)

*DO_modeling_figures* folder contains figures and gifs from forecasts

*DO_modeling_results* contains forecasting results (RMSE, as well as mean and variation for each forecast run)

Forecasts are archived in the *Archived_forecasts* folder
