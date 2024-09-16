%% 08/30/24 DataAnalysis 
% Owner: Taras Hanulia
% Data Type: Cells   Measurements
% Flow Date: 08/30/2024
%% Experiment Note
% BeGFP Cells
% Cells in Media
%% Initialization
clear
clc
addpath 'C:\Users\thanul01\Documents\MATLAB\IVFC_Master'
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'T:\Taras\IVFC\Acquired Data\Cell\TH_083024_Cell';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
output=Labview_convert_rawdata_batch_6PMT(filepath,Fs,Window_Low,Window_High);
disp(output)
close all
%%
%Peak Detection
outputfile= 'NEW_peak_values_08_15_24_T1_2';
file_range= [1,2];
analysisvals=[1,2,3,6];
sample_type= 'Cells';
exp_num=[];
std_threshold=4;
Spectralon_tail= '_1';
FWMH_threshold=0;
intensity_threshold= 0.15;
bead_flag=0;
output=SamplePeakDetection6PMTRF(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
disp(output)
