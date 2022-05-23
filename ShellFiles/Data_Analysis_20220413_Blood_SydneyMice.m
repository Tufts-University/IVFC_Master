%% 041322_DataAnalysis
% Owner: Nilay Vora
% Data Type: Blood Cell Data from Obese Mice
% Flow Date: 04/06/2022
%% Experiment Note
% Obese mouse number is 250, 172, 171 
%% Initialization
clear
clc
%% Calling Script
%%
% Labview Conversion
filepath = 'T:\Nilay\IVFC\Acquired Data\Mouse Studies\NV_051022_BALBc3F16week\NV_051022_BALBc3F16week_Vessel2';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
output=Labview_convert_rawdata_batch(filepath,Fs,Window_Low,Window_High);
disp(output)
close all
%%
%Peak Detection
outputfile= 'NEW_peak_values_04_06_22';
file_range= (1:4);
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
%% 
% Labview Conversion
filepath = 'T:\Nilay\IVFC\Acquired Data\SydneyMice\NV_040622_HealthyControl\Obese';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
output=Labview_convert_rawdata_batch(filepath,Fs,Window_Low,Window_High);
disp(output)
close all

%Peak Detection
outputfile= 'NEW_peak_values_04_06_22';
file_range= (1:4);
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