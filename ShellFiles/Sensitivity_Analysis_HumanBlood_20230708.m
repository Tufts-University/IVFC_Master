%% Human Blood Sensitivity 07/08/23_DataAnalysis
% Owner: Nilay Vora
% Data Type: Full Human Dataset Performance
% Run Date: 07/08/2023
%% Initialization
clear
clc
addpath 'C:/Users/nvora01/documents/matlab/ivfc_master'
%% Calling Script
data_path = 'T:\Taras\IVFC\Acquired Data\Human Studies';
cd(data_path)
%% Data Files Names
performance = zeros(2,3);
performance(:,1) = [1,2]';

header = 'NEW_peak_values';

[sens,pur] = ScatvFLRCalculator_DeepPeak_Human(data_path,header);
performance(2,2:3) = [sens,pur];
performance_final = sortrows(performance,1);
