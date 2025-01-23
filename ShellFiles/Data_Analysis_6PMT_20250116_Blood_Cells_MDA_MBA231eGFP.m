%% 01/16/25 DataAnalysis 
% Owner: Taras Hanulia
% Data Type: Cells   Measurements
% Flow Date: 01/16/2025
%% Experiment Note
% MDA_MB231eGFP Cells in Blood +  beads 
% Cells in Media
%% Initialization
clear
clc
addpath 'C:\Users\thanul01\Documents\MATLAB\IVFC_Master'
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'R:\Taras\IVFC\Blood Cell Data\2025\TH_011625_Blood_Cell';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
% output=Labview_convert_rawdata_batch_6PMT(filepath,Fs,Window_Low,Window_High);
% disp(output)
close all
%%
%Peak Detection
outputfile= 'NEW_peak_values_01_16_25_7';
file_range= [7];
analysisvals=[1,2,4,5,6,7];
sample_type= 'Blood';
exp_num=[];
std_threshold=3;
Spectralon_tail= '7';
FWMH_threshold=0;
intensity_threshold= 0.15;
bead_flag=1;
output=SamplePeakDetection6PMTRF(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
disp(output)
