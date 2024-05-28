%% 231MBACellsSensitivity 05/24/224_DataAnalysis
% Owner: Taras Hanulia
% Data Type: CElls
% Run Date: 05/16/2024
%% Initialization
clear
clc
addpath 'C:\Users\thanul01\Documents\MATLAB\ivfc_master'
%% Calling Script
data_path = 'T:\Taras\IVFC\Acquired Data\Cell\Test';
cd(data_path)
%% Data Files Names
performance = zeros(2,3);
performance(:,1) = [1,2]';

header = 'NEW_peak_values';

[sens,pur] = ScatvFLRCalculator_DeepPeak_Human(data_path,header);
performance(2,2:3) = [sens,pur];
performance_final = sortrows(performance,1);
