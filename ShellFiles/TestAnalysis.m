%% 013122_DataAnalysis
% Owner: Nilay Vora
% Data Type: Blood Cell Data
% Flow Date: 02/15/2022
%% Experiment Note

%% Initialization
clear
clc
%% Calling Script
%% Bead Calibration
% Labview Conversion
filepath = 'T:\Nilay\IVFC\Acquired Data\Blood Cell Data\NV_021522_Blood_Cell';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
% output=Labview_convert_rawdata_batch(filepath,Fs,Window_Low,Window_High);
% disp(output)
% close all
%%
%Peak Detection
outputfile= 'Test_NEW_peak_values_02_15_22';
file_range= (1:3);
analysisvals=(1:4);
sample_type= 'Blood';
exp_num=[];
std_threshold=3;
Spectralon_tail= '_1';
FWMH_threshold=0;
intensity_threshold= 0.1;
bead_flag=1;
output=SamplePeakDetection_PCA_PN(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
disp(output)