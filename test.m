% test the analysis code for meropenem treated S. aureus
clear all; close all; 
p = parpool(12);
foldPath = './20240206_apt_SA_meropenem/';
numChan = 16;
channel_index = 1:16;
col = [0 0 0; 0 0 0; .9 0 0; .9 0 0; 0 0 .9; 0 0 .9; 0 0 .9; 0 .9 0; ...
	0 .9 0; 0 .9 0; .7 .7 0; .7 .7 0; .7 .7 0; 1 0 1; 1 0 1; 1 0 1];
time_index = 1:125;
main
main_DRT_minDelta
close(p)
