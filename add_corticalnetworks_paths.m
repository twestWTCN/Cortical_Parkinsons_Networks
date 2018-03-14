function [] = add_corticalnetworks_paths()
addpath(genpath('C:\Users\Tim\Documents\Work\Cortical_Networks\MRI_processing'))
addpath(genpath('C:\Users\Tim\Documents\Work\Cortical_Networks\phase_amp_analysis'))
addpath(genpath('C:\Users\Tim\Documents\Work\Cortical_Networks\statistics'))
addpath(genpath('C:\Users\Tim\Documents\Work\Cortical_Networks\plotting'))

% Grangerpath
addpath('C:\Users\Tim\Documents\MATLAB_ADDONS\mvgc_v1.0')
run startup.m
addpath('C:\Users\Tim\Documents\MATLAB_ADDONS\gramm-master\@gramm')