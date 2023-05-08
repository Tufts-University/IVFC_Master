%% 05/05/2023 DataAnalysis
% Owner: Nilay Vora
% Data Type: Calibration Measurements
% Flow Date: 05/05/2023
%% Experiment Note
% Added New beads
%% Initialization
clear
clc
addpath 'C:\Users\nvora01\Documents\MATLAB\ivfc_master'
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'T:\Taras\IVFC\Acquired Data\Bead Calibration Data\2023\TH_050523_Calibration';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
output=Labview_convert_rawdata_batch(filepath,Fs,Window_Low,Window_High);
disp(output)
close all

% Peak Detection
outputfile= 'NEW_peak_values_05_05_23';
file_range= (1:3);
analysisvals=(1:3);
sample_type= 'Beads';
exp_num=[];
std_threshold=4;
Spectralon_tail= '';
FWMH_threshold=0;
intensity_threshold= 0.15;
bead_flag=0;
output=SamplePeakDetection(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
disp(output)

% Calibration
files={'NEW_peak_values_02_28_23';...
    'NEW_peak_values_03_01_23';...
    'NEW_peak_values_04_11_23';...
    'NEW_peak_values_04_24_23';...
    'NEW_peak_values_05_04_23';...
    'NEW_peak_values_05_05_23'};
output= DailyCalibrationScript(files);
disp(output)