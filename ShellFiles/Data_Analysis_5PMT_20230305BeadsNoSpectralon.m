%% 09/11/24 DataAnalysis
% Owner: Taras Hanulia
% Data Type: Beads Calibration 
% Flow Date: 03/05/2023 
%% Experiment Note
% beads 
% without Spectralon normalisation
%% Initialization
clear
clc
addpath 'C:\Users\thanul01\Documents\MATLAB\ivfc_master'
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'T:\Nilay\IVFC\Acquired Data\Bead Calibration Data\2023\NV_030523_Calibration';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
% output=Labview_convert_rawdata_batch_6PMT(filepath,Fs,Window_Low,Window_High);
% disp(output)
close all
%%
%Peak Detection
outputfile= 'NEW_peak_values_03_05_23_NS';
file_range= [1];
analysisvals=(1:4);
sample_type= 'Beads';
exp_num=[];
std_threshold=4;
Spectralon_tail= '';
FWMH_threshold=0;
intensity_threshold= 0.15;
bead_flag=1;
output=SamplePeakDetectionNS(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
disp(output)

% Calibration
files={
  
    'NEW_peak_values_03_03_23_NS';...
    'NEW_peak_values_03_04_23_NS';...
    'NEW_peak_values_03_05_23_NS'};
output= DailyCalibrationScript(files);
disp(output)