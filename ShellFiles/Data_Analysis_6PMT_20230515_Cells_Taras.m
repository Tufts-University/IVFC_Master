%% 05/15/23 DataAnalysis
% Owner: Taras Hanulia
% Data Type: Cell Measurements
% Flow Date: 05/15/2023
%% Experiment Note
% T_cell and B Cell
%% Initialization
clear
clc
addpath 'C:\Users\thanul01\Documents\MATLAB\ivfc_master'
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'T:\Nilay\IVFC\Acquired Data\6PMT test\NV_051523_Cell_CART';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
output=Labview_convert_rawdata_batch_6PMT(filepath,Fs,Window_Low,Window_High);
disp(output)
close all
%%
%Peak Detection
outputfile= 'NEW_peak_values_05_15_23';
file_range= (1:2);
analysisvals=(1:3);
sample_type= 'Cells';
exp_num=[];
std_threshold=4;
Spectralon_tail= '';
FWMH_threshold=0;
intensity_threshold= 0.15;
bead_flag=0;
output=SamplePeakDetection6PMT(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
disp(output)

% Calibration
files={'NEW_peak_values_09_24_22';...
    'NEW_peak_values_01_06_23';...
    'NEW_peak_values_01_09_23';...
    'NEW_peak_values_01_10_23';...
    'NEW_peak_values_01_11_23';...
    'NEW_peak_values_02_28_23';...
    'NEW_peak_values_03_01_23';...
    'NEW_peak_values_05_04_23';...
    'NEW_peak_values_05_15_23'};
output= DailyCalibrationScript6PMT(files);
disp(output)