%% 03/03/23 DataAnalysis
% Owner: Nilay Vora
% Data Type: Calibration Measurements
% Flow Date: 03/03/2023
%% Experiment Note
% Channel Clogged
%% Initialization
clear
clc
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'T:\Nilay\IVFC\Acquired Data\Bead Calibration Data\2023\NV_030323_Calibration';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
output=Labview_convert_rawdata_batch(filepath,Fs,Window_Low,Window_High);
disp(output)
close all

%Peak Detection
outputfile= 'NEW_peak_values_03_03_23';
file_range= (1:2);
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
files={'NEW_peak_values_01_11_23';...
    'NEW_peak_values_02_28_23';...
    'NEW_peak_values_03_01_23';...
    'NEW_peak_values_03_02_23';...
    'NEW_peak_values_03_03_23'};
output= DailyCalibrationScript(files);
disp(output)