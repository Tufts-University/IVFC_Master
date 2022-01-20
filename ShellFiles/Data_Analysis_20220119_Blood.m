%% 011222_DataAnalysis
% Owner: Nilay Vora
% Data Type: Calibration Measuremnts
% Flow Date: 01/13/2022
%% Experiment Note

%% Initialization
clear
clc
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'U:\Nilay\IVFC\Acquired Data\Blood Data\NV_011922_Blood_Cell';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
output=Labview_convert_rawdata_batch(filepath,Fs,Window_Low,Window_High);
disp(output)
close all

%Peak Detection
outputfile= 'NEW_peak_values_01_19_22';
file_range= (2:4);
analysisvals=(1:4);
sample_type= 'Blood';
exp_num=[];
std_threshold=3;
Spectralon_tail= '_1';
FWMH_threshold=0;
intensity_threshold= 0.1;
flag=1;
output=SamplePeakDetection(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold,flag);
disp(output)