Contents: 
	MATLAB Function: Basic_ESN_and_PRC_Mix.m
	MATLAB data file: SampleDataJF41.Mat

System Requirements:
Has been run Successfully MATLAB Version 2024a. 
Requires MATLAB License and "Statistics and Machine Learning" toolbox which can be procured from the Mathworks website.

Installation Guide: 
Requires a standard MATLAB software and is run with both files in the same folder. 

Demo and Instructions for use: 
To run the data simply open the Basic_ESN_and_PRC_Mix.m in MATLAB and click run or execute "Basic_ESN_and_PRC_Mix" in the command window. 
The expected output will be 2 plots and 4 data structures in the caller workspace that includes target Jellyfish data motion data and a prediction made with the Hybrid RC. 
	The first plot will show the coefficient of determination for x, y, and z between the predicted and actual velocity of one Jellyfish during stimulation, over a future prediction horizon. 
	The second plot shows the predicted and actual trajectories of the Jellyfish's velocity over time for predictions made at 3 different times into the future. 
	The data structures include the target data and predictions used to make these plots in addition to the predictions made with just the ESN and just the PRC.

