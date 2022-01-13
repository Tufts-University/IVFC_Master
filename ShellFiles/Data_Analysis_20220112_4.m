%% 011222_DataAnalysis
% Owner: Nilay Vora
% Data Type: Calibration Measuremnts
% Flow Date: 01/12/2022
%% Experiment Note
% Limited number of beads were detected, the sample volume is also getting
% low, I will add more water and beads tomorrow to correct concentration
%% Initialization
clear
clc
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'U:\Nilay\IVFC\Acquired Data\Bead Calibration Data\2022\NV_011222_Calibration';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
output=Labview_convert_rawdata_batch(filepath,Fs,Window_Low,Window_High);
disp(output)
close all

Peak Detection
outputfile= 'NEW_peak_values_01_12_22';
file_range= (1:2);
analysisvals=(1:3);
sample_type= 'Beads';
exp_num=[];
std_threshold=4;
Spectralon_tail= '';
FWMH_threshold=0;
intensity_threshold= 0.75;
output=SamplePeakDetection(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold);
disp(output)

% Calibration
files={'NEW_peak_values_11_22_21';...
    'NEW_peak_values_11_29_21';...
    'NEW_peak_values_11_30_21';...
    'NEW_peak_values_12_14_21';...
    'NEW_peak_values_12_20_21';...
    'NEW_peak_values_01_12_22'};
output= DailyCalibrationScript(files);
disp(output)