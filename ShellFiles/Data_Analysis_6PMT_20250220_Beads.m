%% 02/20/25 DataAnalysis
% Owner: Taras Hanulia
% Data Type: Beads Calibration 
% Flow Date: 02/20/2025 
%% Experiment Note
% beads ND 0.9 on 405 and 1.3 on 633, 1.5 on 488 slit of 405 shorter 
% new ND filter for spactralon 2.3 befor 
%% Initialization
clear
clc
addpath 'C:\Users\thanul01\Documents\MATLAB\IVFC_Master'
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'R:\Taras\IVFC\Bead Calibration Data\2025\TH_022025_Calibration';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
output=Labview_convert_rawdata_batch_6PMT(filepath,Fs,Window_Low,Window_High);
disp(output)
close all
%%
% Peak Detection
outputfile= 'NEW_peak_values_02_20_25';
file_range= (1:2);
analysisvals=(1:4);
sample_type= 'Beads';
exp_num=[];
std_threshold=4;
Spectralon_tail= '';
FWMH_threshold=0;
intensity_threshold= 0.15;
bead_flag=1;
output=SamplePeakDetection6PMT(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
disp(output)

%% Calibration
files={
    'NEW_peak_values_02_04_25';...
    'NEW_peak_values_02_06_25';...
    'NEW_peak_values_02_14_25';...
    'NEW_peak_values_02_20_25'};
output= DailyCalibrationScript6PMT(files);
disp(output)