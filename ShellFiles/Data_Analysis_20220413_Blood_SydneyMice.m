%% 041322_DataAnalysis
% Owner: Nilay Vora
% Data Type: Blood Cell Data from Obese
% Flow Date: 04/13/2022
%% Experiment Note
% 8 Week Diet Mice: Obese mouse number is 250, 172, 171
%% Initialization
clear
clc
%% Calling Script
%% Obese Mouse Loop
for i =1:3
    % Labview Conversion
    if i==1
        filepath = 'T:\Nilay\IVFC\Acquired Data\SydneyMice\NV_041322_ObeseCancer';
    else
        filepath = ['T:\Nilay\IVFC\Acquired Data\SydneyMice\NV_041322_ObeseCancer',num2str(i)];
    end
    Fs=60e3;
    Window_Low= 50;
    Window_High= 6000;
    output=Labview_convert_rawdata_batch(filepath,Fs,Window_Low,Window_High);
    disp(output)
    close all
    
    cd(filepath)
    dirinfo = dir();
    dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
    dirinfo(ismember( {dirinfo.name}, {'.', '..'})) = [];  %remove . and ..
    %Peak Detection
    outputfile= 'NEW_peak_values_04_13_22';
    file_range= (1:length(dirinfo));
    clear dirinfo
    analysisvals=(1:4);
    sample_type= 'Blood';
    exp_num=[];
    std_threshold=5;
    Spectralon_tail= '_1';
    FWMH_threshold=0;
    intensity_threshold= 0.1;
    flag=0;
    output=SamplePeakDetection(filepath,outputfile,file_range,Window_Low,...
        Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
        Spectralon_tail,FWMH_threshold,intensity_threshold,flag);
    disp(output)
end