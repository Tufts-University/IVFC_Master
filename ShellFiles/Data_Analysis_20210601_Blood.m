%% 06012021_DataAnalysis
% Owner: Nilay Vora
% Data Type: Blood Cell Data Reprocessing for DeepPeak
%% Initialization
clear
clc
%% Calling Script
mainpath = 'T:\Nilay\IVFC\Acquired Data\Blood Cell Data\';
cd(mainpath)
filepath = 'T:\Nilay\IVFC\Acquired Data\Blood Cell Data\NV_060121_BloodCell';
Fs=60e3;
Window_Low= 50;
Window_High= 6000;
output=Labview_convert_rawdata_batch(filepath,Fs,Window_Low,Window_High);
disp(output)
%% Peak Detection
outputfile = 'NEW_peak_values_06_01_21';
file_range = (1:4);
analysisvals =(1:4);
sample_type = 'Blood';
exp_num = 6000/9;
std_threshold = 3;
Spectralon_tail = '_1';
FWMH_threshold = 0;
intensity_threshold = 0.1;
bead_flag =1;
output = SamplePeakDetection_PCA_PN(filepath,outputfile,file_range,Window_Low,...
        Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
        Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
disp(output)