%% 042325_DataAnalysis
% Owner: Taras Hanulia
% Data Type: BiRFP Data Reprocessing for normalization of raw data
% of the peak
%% Initialization
clear
clc
%% Calling Script
data_path = 'R:\Taras\IVFC\Blood Cell Data\2025\BiRFP';
cd(data_path)
T = readtable("TH_04222025_DatasetSummary_BiRFP.xlsx");
fname_dict = T.FolderName;
fname_dict = cellfun(@(x) x(2:end-1), fname_dict, 'UniformOutput', false);
MouseFlag = T.MouseFlag;
BeadFlag = T.BeadFlag;
idx = find(MouseFlag==1 & T.inc==1);
Final_Fname = {fname_dict{idx}}';
Final_BeadFlag = BeadFlag(idx);

mainpath = 'R:\Taras\IVFC\Blood Cell Data\2025\BiRFP';
cd(mainpath)

parfor i=1:length(Final_Fname)
    disp(['Reprocessing Day ',num2str(i),' of ',num2str(length(Final_Fname))])
    exp_name=Final_Fname{i};
    mydir = fullfile([mainpath,'\',exp_name]);
    cd(mydir)
    filepath=cd;
    Fs=60e3;
    Window_Low= 50;
    Window_High= 10000;
    date=exp_name(4:9);
    date = strread(date,'%2s');
    %% Peak Detection
%     outputfile= ['NTH_peak_values_',date{1},'_',date{2},'_',date{3}];
    dirinfo = dir('*Cell*');
    dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
    dirinfo(ismember( {dirinfo.name}, {'.', '..'})) = [];  %remove . and ..
    file_range= (1:length(dirinfo));
%     analysisvals=[1,2,4];
    sample_type= 'Blood';
    exp_num=round(6000/9);
%     std_threshold=3;
    Spectralon_tail= '';
%     FWMH_threshold=0;
%     intensity_threshold= 0.1;
%     bead_flag = Final_BeadFlag(i);
    output=SamplePeakDetection_PCA_PN_TH_onlynormalization6PMT(filepath,file_range,Window_Low,...
        Window_High,Fs,sample_type,exp_num,Spectralon_tail);
    disp(output)
end