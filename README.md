# Code for the paper "A Low Rank and Sparse Penalized Signal Decomposition Model with Constraints: Anomaly Detection in PV Systems".

This repository provides a MATLAB implementation of the paper. 

- `signal_decomposition.m`: Function to implement the signal decomposition for each week of data
- `signal_decomposition_online_monitoring.m`: Codes to implement the online monitoring using signal decomposition approach proposed
- folder 'input': 
	- folder 'input1': contains all scenarios of faults and normal data (used in 'signal_decomposition.m')
	- folder 'input2': contains one scenario of fault for the demonstration of online monitoring (used in 'signal_decomposition_online_monitoring.m')
	- 'normal6.mat': first 6 normal days used for implementing online monitoring (used in 'signal_decomposition_online_monitoring.m')



