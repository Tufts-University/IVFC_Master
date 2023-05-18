%% 051023_DataAnalysis
% Owner: Nilay Vora
% Data Type: Blood Cell Data from Human Patient
% Flow Date: 03/14/22
%% Notes
% TMC-02 Patient Blood + Beads
%% Initialization
clear
clc
addpath 'C:\Users\nvora01\Documents\MATLAB\ivfc_master'
%% Calling Script
data_path = 'T:\Nilay\IVFC\Acquired Data\Human Studies\TMC-02';
cd(data_path)
Final_Fname = {'NV_051123_Patient2Tiro_TMC'}';

mainpath = 'T:\Nilay\IVFC\Acquired Data\Human Studies\TMC-02';
cd(mainpath)

for i=1:length(Final_Fname)
    disp(['Processing Day ',num2str(i),' of ',num2str(length(Final_Fname))])
    exp_name = Final_Fname{i};
    mydir = fullfile([mainpath,'\',exp_name])
    cd(mydir)
    filepath=cd;
    Fs=60e3;
    Window_Low= 50;
    Window_High= 10000;
    % Labview Conversion
%     output=Labview_convert_rawdata_batch(filepath,Fs,Window_Low,Window_High);
%     disp(output)
    close all
    date=exp_name(4:9);
    date = strread(date,'%2s');
    %% Peak Detection
    outputfile= ['NEW_peak_values_',date{1},'_',date{2},'_',date{3}];
    dirinfo = dir('*TMC*');
    dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
    dirinfo(ismember( {dirinfo.name}, {'.', '..'})) = [];  %remove . and ..
    file_range= (1:length(dirinfo));
    analysisvals=[1:4];
    sample_type= 'Blood';
    exp_num=round(6000/9);
    std_threshold=3;
    Spectralon_tail= '_1';
    FWMH_threshold=0;
    intensity_threshold= 0.1;
    bead_flag = 1;
    output=SamplePeakDetection_PCA_PN(filepath,outputfile,file_range,Window_Low,...
        Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
        Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
    disp(output)
end