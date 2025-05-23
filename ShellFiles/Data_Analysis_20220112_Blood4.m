%% 011222_DataAnalysis
% Owner: Nilay Vora
% Data Type: Cells in Whole Blood with Beads
% Flow Date: 11/30/2021
%% Experiment Note
% Clotting was seen in T3, could affect results
%% Initialization
clear
clc
%% Calling Script
% Labview Conversion
filepath = 'U:\Nilay\IVFC\Acquired Data\Blood Data\NV_113021_Blood_Cell';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
output=Labview_convert_rawdata_batch(filepath,Fs,Window_Low,Window_High);
disp(output)

%Peak Detection
outputfile= 'NEW_peak_values_11_30_21';
file_range= (1:4);
analysisvals=(1:4);
sample_type= 'Blood';
exp_num=[];
std_threshold=3;
Spectralon_tail= '_1';
FWMH_threshold=0;
intensity_threshold= 0.1;
output=SamplePeakDetection(filepath,outputfile,file_range,Window_Low,...
    Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
    Spectralon_tail,FWMH_threshold,intensity_threshold);
disp(output)