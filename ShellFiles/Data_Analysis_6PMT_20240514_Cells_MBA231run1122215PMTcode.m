%% 05/14/23 DataAnalysis 051424
% Owner: Taras Hanulia
% Data Type: Cells   Measurements
% Flow Date: 11/22/2021
%% Experiment Note
% MBA231 GFP+
% Cells in Media 5 PMT code run 5PMT for old data
%% Initialization
clear
clc
addpath 'C:\Users\thanul01\Documents\MATLAB\IVFC_Master'
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'T:\Nilay\IVFC\Acquired Data\231GFP+ Cell Data\NV_112221_231GFP+';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
% output=Labview_convert_rawdata_batch_6PMT(filepath,Fs,Window_Low,Window_High);
% disp(output)
close all
%%
%Peak Detection
outputfile= 'NEW_peak_values_05_14_24_5PMTPCA';
file_range= (1);
analysisvals=[1,2,3,4];
sample_type= 'Cells';
exp_num=[];
std_threshold=4;
Spectralon_tail= '_1';
FWMH_threshold=0;
intensity_threshold= 0.1;
bead_flag=0;
output=SamplePeakDetection_PCA_PN(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
disp(output)
