%% 03/10/23_DataAnalysis
% Owner: Nilay Vora
% Data Type: Full Dataset Threshold Performance
% Run Date: 06/4/2023
%% Initialization
clear
clc
addpath 'C:/Users/nvora01/documents/matlab/ivfc_master'
%% Calling Script
data_path = 'T:\Nilay\IVFC\Acquired Data\Blood Cell Data\DeepPeak2023';
cd(data_path)
%% Data Files Names
performance = zeros(2,3);
performance(:,1) = [1,2]';

header = 'Old_peak_values';

% [sens,pur] = ScatvFLRCalculator_DeepPeak(data_path,header);
% performance(1,2:3) = [sens,pur];
% cd(data_path)

header = 'NEW_peak_values';

[sens,pur] = ScatvFLRCalculator_DeepPeak(data_path,header);
performance(2,2:3) = [sens,pur];
performance_final = sortrows(performance,1);
