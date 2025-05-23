%% 062723_DataAnalysis
% Owner: Taras Hanulia
% Data Type: Blood Cell Data from Human Patient
% Flow Date: 06/09/23
%% Experiment Note
% TMC-03
%% Initialization
clear
clc
%% Calling Script
mainpath = 'T:\Taras\IVFC\Acquired Data\Human Studies\TMC-03\TH_060923_Patient3Tiro_TMC';
cd(mainpath)
folders = dir('*TMC');
for i = 1:length(folders)
    filepath = [mainpath,'\',folders(i).name];
    % Labview Conversion
    Fs=60e3;
    Window_Low= 50;
    Window_High= 10000;
    output=Labview_convert_rawdata_batch(filepath,Fs,Window_Low,Window_High);
    disp(output)
    close all
    
    cd(filepath)
    dirinfo = dir();
    dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
    dirinfo(ismember( {dirinfo.name}, {'.', '..'})) = [];  %remove . and ..
    exp_name = dirinfo(1).name;
    date = exp_name(4:9);
    date = strread(date,'%2s');
    %Peak Detection
    outputfile= ['NEW_peak_values_',date{1},'_',date{2},'_',date{3}];
    file_range= (1:length(dirinfo));
    clear dirinfo
    analysisvals=(1:4);
    sample_type= 'Blood';
    exp_num=[];
    std_threshold=3;
    Spectralon_tail= '_1';
    FWMH_threshold=0;
    intensity_threshold= 0.1;
    flag=0;
    output=SamplePeakDetection6PMT(filepath,outputfile,file_range,Window_Low,...
        Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
        Spectralon_tail,FWMH_threshold,intensity_threshold,flag);
    disp(output)
end