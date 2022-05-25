%% 05/24/22_DataAnalysis
% Owner: Nilay Vora
% Data Type: Calibration Measuremnts
% Flow Date: 05/24/2022
%% Experiment Note
% NEW Detection parameters
%   PMT 1,2,5 Gain increased to 2.4 with ND~0.6 added to each detector
%   PMT 3,4 Gain reduced to 3.2, 2.2 respectively with ND=0.3 removed
%% Initialization
clear
clc
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'T:\Nilay\IVFC\Acquired Data\Bead Calibration Data\2022\NV_052422_Calibration';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
output=Labview_convert_rawdata_batch(filepath,Fs,Window_Low,Window_High);
disp(output)
close all
%%
%Peak Detection
outputfile= 'NEW_peak_values_05_24_22';
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

%% Calibration
files={'NEW_peak_values_04_13_22';...
    'NEW_peak_values_04_14_22';...
    'NEW_peak_values_05_10_22';...
    'NEW_peak_values_05_12_22';...
    'NEW_peak_values_05_23_22';...
    'NEW_peak_values_05_24_22'};
output= DailyCalibrationScript(files);
disp(output)