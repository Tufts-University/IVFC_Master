%% 09/30/22_DataAnalysis
% Owner: Nilay Vora
% Data Type: Full Dataset Retraining
% Run Date: 02/17/2023
%% Initialization
clear
clc
%% Calling Script
%% Data Files Names
% Labview Conversion
mainpath = 'T:\Nilay\IVFC\Acquired Data\Blood Cell Data\CNNData';
cd(mainpath)
[~,T] = xlsread('dict_key.csv');
T(1)=[];
for i=31:length(T)
    disp(['Reprocessing Day ',num2str(i),' of 30'])
    exp_name=T{i};
    cd(exp_name(2:end))
    filepath=cd;
    Fs=60e3;
    Window_Low= 50;
    Window_High= 10000;
    date=exp_name(5:10);
    date = strread(date,'%2s');
    %% Peak Detection
    outputfile= ['Corrected_012923_peak_values_',date{1},'_',date{2},'_',date{3}];
    dirinfo = dir('*Cell*');
    dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
    dirinfo(ismember( {dirinfo.name}, {'.', '..'})) = [];  %remove . and ..
    file_range= (1:length(dirinfo));
    analysisvals=[1];
    sample_type= 'Blood';
    exp_num=round(6000/9);
    std_threshold=3;
    if i==21 || i==23
         Spectralon_tail= '_2';
    elseif i==22
         Spectralon_tail= '_3';
    else
         Spectralon_tail= '_1';
    end
    FWMH_threshold=0;
    intensity_threshold= 0.1;
    if i>=6 && i<=13
        bead_flag=0;
        analysisvals=[1];
    else
        bead_flag=1;
    end
    output=SamplePeakDetection_PCA_PN(filepath,outputfile,file_range,Window_Low,...
        Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
        Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
    disp(output)
    cd(mainpath)
end