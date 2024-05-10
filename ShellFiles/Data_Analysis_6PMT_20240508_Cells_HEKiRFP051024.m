%% 05/10/23 DataAnalysis 
% Owner: Taras Hanulia
% Data Type: Cells   Measurements
% Flow Date: 05/08/2024
%% Experiment Note
% HEKiRFP ,first two experiment T1and T2 is HEKGFP , the rest is HEK iRFP
% Cells in Media, new function 6PMTRF
%% Initialization
clear
clc
addpath 'C:\Users\thanul01\Documents\MATLAB\ivfc_masterr'
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'T:\Taras\IVFC\Acquired Data\Cell\TH_050824_Cell';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
% output=Labview_convert_rawdata_batch_6PMT(filepath,Fs,Window_Low,Window_High);
% disp(output)
close all
%%
%Peak Detection
outputfile= 'NEW1_peak_values_05_10_24';
file_range= (3:7);
analysisvals=[1,2,3,4,5,6];
sample_type= 'Cells';
exp_num=[1500];
std_threshold=3;
Spectralon_tail= '';
FWMH_threshold=0;
intensity_threshold= 0.1;
bead_flag=0;
output=SamplePeakDetection6PMTRF(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
disp(output)
