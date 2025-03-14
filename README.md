# cartilage-impact-analysis
MATLAB scripts for cartilage impact data analysis
This MATLAB script calculates impact parameters for mechanical loading experiments on tissue samples. It processes sensor data from loadcell and acceleration measurements to extract key mechanical properties of the impact event. Additionally, it generates plots for force, pressure, acceleration, velocity, deformation, and hysteresis curves.
# Required dependencies
1.MATLAB version R2018a or later
2.Signal Processing Toolbox 
3.Curve Fitting Toolbox
# Brief Summary
There are 2 inputs for this code (voltage from load cell and acceleoemeter). First, it applies sensor scaling factors to convert raw voltage readings into physical units. The acceleration scale (accelscale) converts voltage readings to acceleration in meters per second squared (m/s²), while the force scale (loadscale) converts voltage readings to force in Newtons (N).
The script generates several key mechanical parameters that describe the impact event. It calculates the peak force in both raw and filtered forms, along with the peak acceleration during impact. The impact duration is determined based on the time interval between the start and end of the impact event.

# How to Use
Place the sensor data files (.txt) in the working directory.
Run the script in MATLAB.
The script will automatically process all .txt files and display the results.
Plots will be generated for each dataset.
