%% 01/10/25 DataAnalysis 
% Owner: Taras Hanulia
% Data Type: Cells   Measurements
% Flow Date: 01/09/2025
%% Experiment Note
% BiRFP and GFP 1:5 ratio 
% Cells in Media T1-4 BiRFP_CARTeGFP 20k 100k ,T5-10 BiRFP_CARTeGFP 20k 300k
%% Initialization
clear
clc
addpath 'C:\Users\thanul01\Documents\MATLAB\IVFC_Master'
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'R:\Taras\IVFC\CellMedia\BiRFP_CARTeGFP\TH_010925_Cell';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
% output=Labview_convert_rawdata_batch_6PMT(filepath,Fs,Window_Low,Window_High);
% disp(output)
close all
%%
%Peak Detection
outputfile= 'NEW_peak_values_01_06_25_5_10';
file_range= (5:10);
analysisvals=[1,2,4,5,6,7];
sample_type= 'Cells';
exp_num=[];
std_threshold=5;
Spectralon_tail= '';
FWMH_threshold=0;
intensity_threshold= 0.2;
bead_flag=0;
output=SamplePeakDetection6PMTRF(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
disp(output)
