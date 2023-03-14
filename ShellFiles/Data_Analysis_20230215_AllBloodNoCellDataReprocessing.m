%% 02/15/23_DataAnalysis
% Owner: Nilay Vora
% Data Type: Full Dataset Reprocessing
% Run Date: 02/27/2023
%% Initialization
clear
clc
%% Calling Script
%% Data Files Names
mainpath = 'T:\Nilay\IVFC\Acquired Data\Blood Cell Data';
cd(mainpath)
T = dir('*_noCell');
for i=1:length(T)
    disp(['Reprocessing Day ',num2str(i),' of ' num2str(length(T))])
    exp_name=T(i).name;
    cd(exp_name)
    filepath=cd;
    Fs=60e3;
    Window_Low= 50;
    Window_High= 10000;
    date=exp_name(4:10);
    date = strread(date,'%2s');
    %% Peak Detection
    outputfile= ['Corrected_012923_peak_values_',date{1},'_',date{2},'_',date{3}];
    dirinfo = dir('*Cell*');
    dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
    dirinfo(ismember( {dirinfo.name}, {'.', '..','Local Plots'})) = [];  %remove . and ..
    file_range= (1:length(dirinfo));
    analysisvals=[1];
    sample_type= 'Blood';
    exp_num=round(6000/9);
    std_threshold=3;
    Spectralon_tail= '_1';
    FWMH_threshold=0;
    intensity_threshold= 0.1;
    bead_flag=0;
    output=SamplePeakDetection_PCA_PN(filepath,outputfile,file_range,Window_Low,...
        Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
        Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
    disp(output)
    cd(mainpath)
end