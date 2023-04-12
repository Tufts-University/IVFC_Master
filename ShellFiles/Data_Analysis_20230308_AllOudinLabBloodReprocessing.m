%% 03/08/23_DataAnalysis
% Owner: Nilay Vora
% Data Type: Full Dataset Reprocessing for Oudin Lab
% Run Date: 03/08/23
% Modified: 04/07/23-> Reprocessing RFLR Data without cluster analysis
%% Initialization
clear
clc
addpath(genpath('C:\Users\nvora01\Documents\MATLAB\IVFC'))
addpath(genpath('C:\Users\nvora01\Documents\MATLAB\ivfc_master\'))
%% Calling Script
mainpath = 'T:\Nilay\IVFC\Acquired Data\SydneyMice';
cd(mainpath)
%% Healthy Controls
T = dir('*Control');
for i=1:length(T)
    disp(['Reprocessing Day ',num2str(i),' of ' num2str(length(T))])
    exp_name=T(i).name;
    cd(exp_name)
    outpath=cd;
    Fs=60e3;
    Window_Low= 50;
    Window_High= 10000;
    date=exp_name(4:10);
    date = strread(date,'%2s');
    %% Peak Detection
    outputfile= ['Corrected_030823_peak_values_',date{1},'_',date{2},'_',date{3}];
    dirinfo = dir();
    dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
    dirinfo(ismember( {dirinfo.name}, {'.', '..','Local Plots'})) = [];  %remove . and ..
    for f = 1:2
        data_type = dirinfo(f).name;
        cd(data_type)
        filepath=cd;
        subdirinfo = dir('*Control*');
        subdirinfo(~[subdirinfo.isdir]) = [];  %remove non-directories
        subdirinfo(ismember( {subdirinfo.name}, {'.', '..','Local Plots'})) = [];  %remove . and ..
        file_range= (1:length(subdirinfo));
        analysisvals=[4];
        sample_type= 'Blood';
        exp_num=round(6000/9);
        std_threshold=10;
        Spectralon_tail= '_1';
        FWMH_threshold=0;
        intensity_threshold= 0.0;
        bead_flag=0;
        output=SamplePeakDetection_PCA_PN(filepath,outputfile,file_range,Window_Low,...
            Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
            Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
        disp(output)
        cd(outpath)
    end
    cd(mainpath)
end
%% Cancer Studies
T = dir('*Cancer*');
parfor i=1:length(T)
    disp(['Reprocessing Day ',num2str(i),' of ' num2str(length(T))])
    exp_name=T(i).name;
    exp_path = [T(i).folder,'\',T(i).name];
    cd(exp_path)
    filepath=cd;
    Fs=60e3;
    Window_Low= 50;
    Window_High= 10000;
    date=exp_name(4:10);
    date = strread(date,'%2s');
    %% Peak Detection
    outputfile= ['Corrected_030823_peak_values_',date{1},'_',date{2},'_',date{3}];
    dirinfo = dir('*Cancer*');
    dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
    dirinfo(ismember( {dirinfo.name}, {'.', '..','Local Plots'})) = [];  %remove . and ..
    file_range= (1:length(dirinfo));
    analysisvals=[4];
    sample_type= 'Blood';
    exp_num=round(6000/9);
    std_threshold=10;
    Spectralon_tail= '_1';
    FWMH_threshold=0;
    intensity_threshold= 0.0;
    bead_flag=0;
    output=SamplePeakDetection_PCA_PN(filepath,outputfile,file_range,Window_Low,...
        Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
        Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
    disp(output)
    cd(mainpath)
end