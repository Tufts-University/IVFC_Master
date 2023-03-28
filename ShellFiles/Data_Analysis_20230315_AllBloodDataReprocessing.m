%% 031522_DataAnalysis
% Owner: Nilay Vora
% Data Type: Blood Cell Data Reprocessing for DeepPeak
%% Initialization
clear
clc
%% Calling Script
mainpath = 'T:\Nilay\IVFC\Acquired Data\Blood Cell Data\DeepPeak2023';
cd(mainpath)
fname_dict = readtable("dict_key.csv",'Delimiter','''');
T = table2cell(fname_dict(:,2));
parfor i=1:length(T)
    disp(['Reprocessing Day ',num2str(i),' of ',num2str(length(T))])
    exp_name=T{i};
    mydir = fullfile([mainpath,'\',exp_name])
    cd(mydir)
    filepath=cd;
    Fs=60e3;
    Window_Low= 50;
    Window_High= 10000;
    date=exp_name(4:9);
    date = strread(date,'%2s');
    %% Peak Detection
    outputfile= ['NEW_peak_values_',date{1},'_',date{2},'_',date{3}];
    dirinfo = dir('*Cell*');
    dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
    dirinfo(ismember( {dirinfo.name}, {'.', '..'})) = [];  %remove . and ..
    file_range= (1:length(dirinfo));
    analysisvals=[1,2,4];
    sample_type= 'Blood';
    exp_num=round(6000/9);
    std_threshold=3;
    Spectralon_tail= '_1';
    FWMH_threshold=0;
    intensity_threshold= 0.1;
    if i>=6 && i<=10 || i>31
        bead_flag=0;
    else
        bead_flag=1;
    end
    output=SamplePeakDetection_PCA_PN(filepath,outputfile,file_range,Window_Low,...
        Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
        Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
    disp(output)
end