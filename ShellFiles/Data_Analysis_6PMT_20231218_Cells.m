%% 12/18/23 DataAnalysis 
% Owner: Taras Hanulia
% Data Type: Cells   Measurements
% Flow Date: 12/18/2023
%% Experiment Note
% B Cells 
% Cells in Media B cells+ iRFP 1-4 , B cells 5-8
%% Initialization
clear
clc
addpath 'C:\Users\thanul01\Documents\MATLAB\ivfc_masterr'
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'T:\Taras\IVFC\Acquired Data\Cell\TH_121823_Cell';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
% output=Labview_convert_rawdata_batch_6PMT(filepath,Fs,Window_Low,Window_High);
% disp(output)
close all
%%
%Peak Detection
outputfile= 'NEW_peak_values_12_18_23';
file_range= (1:8);
analysisvals=[4,5];
sample_type= 'Cells';
exp_num=[];
std_threshold=3;
Spectralon_tail= '';
FWMH_threshold=0;
intensity_threshold= 0.3;
bead_flag=0;
output=SamplePeakDetection6PMTRF_PCA_PN(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
disp(output)
