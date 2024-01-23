%% 01/23/23 DataAnalysis 
% Owner: Taras Hanulia
% Data Type: Cells   Measurements
% Flow Date: 01/18/2024
%% Experiment Note
% B Cells +iRFP720 Cells5 minus 1 ND filter on 633 exitation line
% Cells in Media
%% Initialization
clear
clc
addpath 'C:\Users\thanul01\Documents\MATLAB\ivfc_masterr'
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'T:\Taras\IVFC\Acquired Data\Cell\TH_011824_Cell\BCells_iRFP\BiRFP5';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
% output=Labview_convert_rawdata_batch_6PMT(filepath,Fs,Window_Low,Window_High);
% disp(output)
close all
%%
%Peak Detection
outputfile= 'NEW_peak_values_01_23_23';
file_range= [1];
analysisvals=(1:6);
sample_type= 'Cells';
exp_num=[];
std_threshold=3;
Spectralon_tail= '';
FWMH_threshold=0;
intensity_threshold= 0.15;
bead_flag=0;
output=SamplePeakDetection6PMTRF_PCA_PN(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
disp(output)
