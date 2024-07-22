%% 07/18/2024
% Owner: Taras Hanulia
% Data Type: Blood    Measurements
% Flow Date: 07/15/2023
%% Experiment Note
% MDA_MB231 GFP Cells 
% flow faster than expected
%% Initialization
clear
clc
addpath 'C:\Users\thanul01\Documents\MATLAB\ivfc_masterr'
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'T:\Taras\IVFC\Acquired Data\Blood Cell Data\SNR\TH_071524_Blood_Cell';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
output=Labview_convert_rawdata_batch_6PMT(filepath,Fs,Window_Low,Window_High);
disp(output)
close all
%%
%Peak Detection
outputfile= 'NEW_peak_values_07_15_23_5s';
file_range= (1:4);
analysisvals=[1,2,3,6];
sample_type= 'Blood';
exp_num=[];
std_threshold=5;
Spectralon_tail= '';
FWMH_threshold=0;
intensity_threshold= 0.15;
bead_flag=1;
output=SamplePeakDetection6PMTRF_PCA_PN(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
disp(output)
