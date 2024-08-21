%% 08/17/24 DataAnalysis 
% Owner: Taras Hanulia
% Data Type: Cells   Measurements
% Flow Date: 08/17/2024
%% Experiment Note
% B16F10 CART ratio 1:3 T9 after 4 hour of mixture
% Cells in Media
%% Initialization
clear
clc
addpath 'C:\Users\thanul01\Documents\MATLAB\ivfc_master'
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'T:\Taras\IVFC\Acquired Data\Cell\TH_081724_Cell';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
% output=Labview_convert_rawdata_batch_6PMT(filepath,Fs,Window_Low,Window_High);
% disp(output)
close all
%%
%Peak Detection
outputfile= 'NEW_peak_values_08_09_24_T9';
file_range= [9];
analysisvals=[1,2,3];
sample_type= 'Cells';
exp_num=[];
std_threshold=4;
Spectralon_tail= '_2';
FWMH_threshold=0;
intensity_threshold= 0.15;
bead_flag=0;
output=SamplePeakDetection6PMTRF(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
disp(output)
