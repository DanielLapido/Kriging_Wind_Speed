# Wind prediction based on aircraft data using Kriging

The goal of this project is to explore the applicability of a geostatistical
technique called Kriging for short term wind speed prediction. The predictions are based on spatio-temporal correlation among wind observations
collected by aircraft acting as meteorological sensors. The objective is to
employ data from the last 3 minutes and provide reliable predictions for the
next minutes.

Data collected along the aircraft trajectories:
<img src="https://github.com/DanielLapido/Kriging_Wind_Speed/blob/main/Figures/datapoints.jpeg" width=50% height=50%>

Visualization of the wind vector field:
<img src="https://github.com/DanielLapido/Kriging_Wind_Speed/blob/main/Figures/arrows.jpg" width=50% height=50%>

The Wind is decomposed into a North speed component and an East speed component by transforming the polar coordinates into cartesian coordinates.
For each one of them, it is being assumed that the Speed is a random function that can be expressed as a sum of a trend component and a stationary residual component with zero mean.

The model is a combination of a Machine Learning model (Random Forest) that estimates the trend component and a geostatistical model (Kriging) that estimate the residual component at unknown locations.

The Random Forest model is responsible for the large scale variation of the wind. If the trend estimated with the Random Forest is removed from the data, we obtain the residuals at the known locations. The residuals at nearby locations tend to be more similar that those far apart. This is called spatio-temporal correlation and Kriging takes advantage of it to get more precise predictions.

Spatio-temporal correlation of the residuals:
![](https://github.com/DanielLapido/Kriging_Wind_Speed/blob/main/Figures/north_variogram.jpeg)


The spatio-temporal model trained with just 3 minutes provides good predictions for the next minutes:

First 10 seconds:
![](https://github.com/DanielLapido/Kriging_Wind_Speed/blob/main/Figures/presentation_rfwindrose10.jpeg)

Next minute:
![](https://github.com/DanielLapido/Kriging_Wind_Speed/blob/main/Figures/rf_1mrose.jpeg)

Next 7 minutes:
![](https://github.com/DanielLapido/Kriging_Wind_Speed/blob/main/Figures/rf5mrose.jpeg)
