%% 040722_DataAnalysis
% Owner: Nilay Vora
% Data Type: Blood Cell Data
% Flow Date: 04/07/2022
%% Experiment Note
% Two groups of control mice data collected: Lean and Obese Diets no tumor
% cells.
% Lean mouse number is 245, Obese mouse number is 180
%% Initialization
clear
clc
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'T:\Nilay\IVFC\Acquired Data\SydneyMice\NV_040722_HealthyControl\Lean';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
output=Labview_convert_rawdata_batch(filepath,Fs,Window_Low,Window_High);
disp(output)
close all

%% Peak Detection
outputfile= 'NEW2_peak_values_04_07_22';
file_range= (1:4);
analysisvals=(1:4);
sample_type= 'Blood';
exp_num=[];
std_threshold=5;
Spectralon_tail= '_1';
FWMH_threshold=0;
intensity_threshold= 0.1;
flag=1;
output=SamplePeakDetection(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold,flag);
disp(output)

%% Bead Calibration
% Labview Conversion
filepath = 'T:\Nilay\IVFC\Acquired Data\SydneyMice\NV_040722_HealthyControl\Obese';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
output=Labview_convert_rawdata_batch(filepath,Fs,Window_Low,Window_High);
disp(output)
close all

%% Peak Detection
outputfile= 'NEW2_peak_values_04_07_22';
file_range= (1:4);
analysisvals=(1:4);
sample_type= 'Blood';
exp_num=[];
std_threshold=5;
Spectralon_tail= '_1';
FWMH_threshold=0;
intensity_threshold= 0.1;
flag=1;
output=SamplePeakDetection(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold,flag);
disp(output)