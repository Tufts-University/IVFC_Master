%% 07/10/24 DataAnalysis
% Owner: Taras Hanulia
% Data Type: Beads Calibration 
% Flow Date: 07/10/2024 
%% Experiment Note
% beads ND 1.2 on 405 and 1.3 on 633, 1.5 on 488 slit of 405 shorter 
% new ND filter for spactralon 2.3 befor new ex ND second run new channel
%% Initialization
clear
clc
addpath 'C:\Users\thanul01\Documents\MATLAB\ivfc_master'
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'T:\Taras\IVFC\Acquired Data\Bead Calibration Data\2024\TH_071024_Calibration2';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
output=Labview_convert_rawdata_batch_6PMT(filepath,Fs,Window_Low,Window_High);
disp(output)
close all
%%
%Peak Detection
outputfile= 'NEW_peak_values_07_10_24_2';
file_range= (1:3);
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

% Calibration
files={
    'NEW_peak_values_07_08_24';...
    'NEW_peak_values_07_09_24';...
    'NEW_peak_values_07_09_24_2';...
    'NEW_peak_values_07_10_24';...
    'NEW_peak_values_07_10_24_2'};
output= DailyCalibrationScript6PMT(files);
disp(output)