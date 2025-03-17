# Cartilage-impact-analysis
The Cartilage Impact Analysis code is a MATLAB script designed to analyze data from mechanical impact experiments on bovine cartilage samples. It processes sensor measurements from a load cell and accelerometer to calculate key mechanical properties of the impact event, including impact energy, stress, impulse, deformation, and loading rate.

# Required dependencies
1.MATLAB version R2018a or later
2.Signal Processing Toolbox 
3.Curve Fitting Toolbox

# Brief Summary
The script reads raw sensor data, including time and voltage from the load cell and accelerometer, from .txt files and uses parameters such as impactor properties (mass and diameter), sample properties (diameter and thickness), Impactor mass and sensor scaling factors for force and acceleration conversion. First, it applies these scaling factors to convert raw voltage data into physical measurements, where the acceleration scale (accelscale) converts voltage into acceleration (m/s²) and the force scale (loadscale) converts voltage into force (N). For each dataset, the script calculates and reports Peak Force (raw and filtered), Peak Acceleration, Impact Duration, Impulse (N·s), Work (J), Impact Energy (J), Maximum Stress (MPa), Loading Rate (MPa/ms), and Maximum Deformation (mm). 
 
# How to Use
Place the sensor data files (.txt) in the working directory. Run the script in MATLAB. The script will automatically process all .txt files and display the results. Plots will be generated for each dataset.

# Getting started 
•	Run one of the examples (test 1- test 4)!

For more information, please refer to the paper:

A Reproducible Cartilage Impact Model to Generate Post-Traumatic Osteoarthritis in the Rabbit (https://app.jove.com/v/64450/a-reproducible-cartilage-impact-model-to-generate-post-traumatic)
