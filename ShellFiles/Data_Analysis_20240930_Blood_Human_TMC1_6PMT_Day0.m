%% 100124_DataAnalysis
% Owner: Taras Hanulia
% Data Type: Blood CART Cell Data from Human Patient
% Flow Date: 09/30/24
%% Notes
% TMC-1 Day0 Blood
%% Initialization
clear
clc
addpath 'C:\Users\thanul01\Documents\MATLAB\ivfc_master'
%% Calling Script
    filepath = 'T:\Taras\IVFC\Acquired Data\Human studies\TH_093024_Patient1TMC_0Day';
    % Labview Conversion
    Fs=60e3;
    Window_Low= 50;
    Window_High= 10000;
%     output=Labview_convert_rawdata_batch_6PMT(filepath,Fs,Window_Low,Window_High);
%     disp(output)
    close all
    
    cd(filepath)
    dirinfo = dir();
    dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
    dirinfo(ismember( {dirinfo.name}, {'.', '..'})) = [];  %remove . and ..
    exp_name = dirinfo(1).name;
    date = exp_name(4:9);
    date = strread(date,'%2s');
    %%Peak Detection
    outputfile= ['NEW_peak_values_',date{1},'_',date{2},'_',date{3},'_1_3'];
    dirinfo = dir('*TMC*');
    dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
    dirinfo(ismember( {dirinfo.name}, {'.', '..'})) = [];  %remove . and ..
    file_range= (1:3); %file_range= (1:length(dirinfo));
    analysisvals=[1,2,3,6];
    sample_type= 'Blood';
    exp_num=[];
    std_threshold=5;
    Spectralon_tail='_1';
    FWMH_threshold=0;
    intensity_threshold=0.1;
    bead_flag=0;
    output=SamplePeakDetection6PMTRF_PCA_PN(filepath,outputfile,file_range,Window_Low,...
        Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
        Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
    disp(output)