%% 11/29/23 DataAnalysis
% Owner: Taras Hanulia
% Data Type: Beads Calibration 
% Flow Date: 11/29/2023
%% Experiment Note
% beads 
%% Initialization
clear
clc
addpath 'C:\Users\thanul01\Documents\MATLAB\ivfc_master'
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'T:\Taras\IVFC\Acquired Data\Bead Calibration Data\2023\6PMT\TH_112923_Calibration';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
output=Labview_convert_rawdata_batch_6PMT(filepath,Fs,Window_Low,Window_High);
disp(output)
close all
%%
%Peak Detection
outputfile= 'NEW_peak_values_11_29_23';
file_range= (1:2);
analysisvals=(1:4);
sample_type= 'Beads';
exp_num=[];
std_threshold=4;
Spectralon_tail= '';
FWMH_threshold=0;
intensity_threshold= 0.15;
bead_flag=0;
output=SamplePeakDetection6PMT(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
disp(output)

% Calibration
files={
    'NEW_peak_values_11_09_23';...
    'NEW_peak_values_11_16_23';...
    'NEW_peak_values_11_27_23';...
    'NEW_peak_values_11_28_23';...
    'NEW_peak_values_11_29_23'};
output= DailyCalibrationScript6PMT(files);
disp(output)