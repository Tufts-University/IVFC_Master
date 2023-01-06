%% 09/24/22_DataAnalysis
% Owner: Nilay Vora
% Data Type: Calibration Measuremnts
% Flow Date: 09/24/2022
%% Experiment Note
%
%% Initialization
clear
clc
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'T:\Nilay\IVFC\Acquired Data\Bead Calibration Data\2022\NV_092422_Calibration';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
output=Labview_convert_rawdata_batch(filepath,Fs,Window_Low,Window_High);
disp(output)
close all
%%
%Peak Detection
outputfile= 'NEW_peak_values_09_24_22';
file_range= (1:2);
analysisvals=(1:3);
sample_type= 'Beads';
exp_num=[];
std_threshold=4;
Spectralon_tail= '';
FWMH_threshold=0;
intensity_threshold= 0.75;
bead_flag=0;
output=SamplePeakDetection(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
disp(output)

% Calibration
files={'NEW_peak_values_06_07_22';...
    'NEW_peak_values_06_08_22';...
    'NEW_peak_values_07_06_22';...
    'NEW_peak_values_07_07_22';...
    'NEW_peak_values_09_24_22'};
output= DailyCalibrationScript(files);
disp(output)