%% 03/06/25 DataAnalysis and 03/20/25
% Owner: Taras Hanulia
% Data Type: Cells   Measurements
% Flow Date: 03/05/2025
%% Experiment Note
% BiRFP  Cells in Blood +  beads

% Cells in Blood
%% Initialization
clear
clc
addpath 'C:\Users\thanul01\Documents\MATLAB\IVFC_Master'
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'R:\Taras\IVFC\Blood Cell Data\2025\BiRFP\TH_030525_Blood_Cell';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
% output=Labview_convert_rawdata_batch_6PMT(filepath,Fs,Window_Low,Window_High);
% disp(output)
close all
%%
%Peak Detection
outputfile= 'NEW_peak_values_03_05_25';
file_range=(1:8);
analysisvals=[1,2,4,5,6,7];
sample_type= 'Blood';
exp_num=[150];
std_threshold=3;
Spectralon_tail= '';
FWMH_threshold=0;
intensity_threshold= 0.2;
bead_flag=1;
output=SamplePeakDetection6PMTRF(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
disp(output)
