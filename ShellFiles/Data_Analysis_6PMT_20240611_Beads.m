%% 06/11/24 DataAnalysis
% Owner: Taras Hanulia
% Data Type: Beads Calibration 
% Flow Date: 06/11/2024 
%% Experiment Note
% beads ND 0.6 on 405 and 633, 1.0 on 488 slit of 405 shorter 
% new ND filter for spactralon 2.3 befor was 3.0 new aligmentof 405
% telescope
%% Initialization
clear
clc
addpath 'C:\Users\thanul01\Documents\MATLAB\ivfc_master'
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'T:\Taras\IVFC\Acquired Data\Bead Calibration Data\2024\TH_061124_Calibration';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
output=Labview_convert_rawdata_batch_6PMT(filepath,Fs,Window_Low,Window_High);
disp(output)
close all
%%
%Peak Detection
outputfile= 'NEW_peak_values_06_11_24';
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

% Calibration
files={
  
    'NEW_peak_values_04_26_24';...
    'NEW_peak_values_05_08_24';...
    'NEW_peak_values_06_10_24';...
    'NEW_peak_values_06_11_24'};
output= DailyCalibrationScript6PMT(files);
disp(output)