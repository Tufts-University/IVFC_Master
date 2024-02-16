%% 02/14/23 DataAnalysis 
% Owner: Taras Hanulia
% Data Type: Cells   Measurements
% Flow Date: 02/14/2024
%% Experiment Note
% B Cells +iRFP720 change intensity treshhold 0.1 from 0.15 from green 5 to
% red 4
% Cells in Media
%% Initialization
clear
clc
addpath 'C:\Users\thanul01\Documents\MATLAB\ivfc_masterr'
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'T:\Taras\IVFC\Acquired Data\Cell\TH_021524_Cell';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
% output=Labview_convert_rawdata_batch_6PMT(filepath,Fs,Window_Low,Window_High);
% disp(output)
close all
%%
%Peak Detection
outputfile= 'NEW_peak_values_02_15_23';
file_range= (1:6);
analysisvals=(1:6);
sample_type= 'Cells';
exp_num=[];
std_threshold=3;
Spectralon_tail= '';
FWMH_threshold=0;
intensity_threshold= 0.1;
bead_flag=0;
output=SamplePeakDetection6PMTRF_PCA_PN(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
disp(output)
