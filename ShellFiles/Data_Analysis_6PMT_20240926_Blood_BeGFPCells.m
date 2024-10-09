%% 09/26/2024
% Owner: Taras Hanulia
% Data Type: Blood    Measurements
% Flow Date: 09/26/2023
%% Experiment Note
% B GFP Cells 
% flow faster than expected
%% Initialization
clear
clc
addpath 'C:\Users\thanul01\Documents\MATLAB\ivfc_master'
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'T:\Taras\IVFC\Acquired Data\Blood Cell Data\SNR\TH_092624_Blood_Cell';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
% output=Labview_convert_rawdata_batch_6PMT(filepath,Fs,Window_Low,Window_High);
% disp(output)
close all
%%
%Peak Detection
outputfile= 'NEW_peak_values_09_26_24_PCA';
file_range= (1:6);
analysisvals=[1,2,3,6];
sample_type= 'Blood';
exp_num=[];
std_threshold=5;
Spectralon_tail= '';
FWMH_threshold=0;
intensity_threshold= 0.15;
bead_flag=0;
output=SamplePeakDetection6PMTRF_PCA_PN(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
disp(output)
