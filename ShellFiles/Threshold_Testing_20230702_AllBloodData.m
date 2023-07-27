%% 07/02/23_DataAnalysis
% Owner: Nilay Vora
% Data Type: Full Dataset Threshold Performance
% Run Date: 07/02/2023
%% Initialization
clear
clc
addpath 'C:\Users\nvora01\Documents\MATLAB\ivfc_master'
%% Calling Script
%% Data Files Names
mainpath=['T:\Nilay\IVFC\Acquired Data\Blood Cell Data\DeepPeak2023\'];
cd(mainpath)
file_folders= '*Blood*';
folders = dir(file_folders);
folders(~[folders.isdir]) = [];  %remove non-directories

file_name = ['T:\Nilay\IVFC\Acquired Data\Blood Cell Data\',...
                                'NV_20230328_DatasetSummary.xlsx'];
T = readtable(file_name);
fname_dict = T.FolderName;
fname_dict = cellfun(@(x) x(2:end-1), fname_dict, 'UniformOutput', false);
CTCCFlag = T.CTCCFlag;
BeadFlag = T.BeadFlag;
RatFlag = T.RatFlag;
inc = T.inc;
idx = find(RatFlag==1 & CTCCFlag==1 & inc==1);
Final_Fname = {fname_dict(idx)}';
Final_BeadFlag = BeadFlag(idx);

tempstorage = cell(1,length(folders));
bead_flag = Final_BeadFlag;
n = length(folders);
folders = Final_Fname{1};


threshold = [6,8,10,12,14];
performance = zeros(6,3);
performance(:,1) = [threshold,13]';
for j = 1:5
    thresh = threshold(j);
    disp(['Threshold ',num2str(j),' of 5'])
%     parfor i=1:length(folders)
%         disp([9 'Reprocessing Day ',num2str(i),' of 30'])
%         exp_name=folders{i};
%         cd(exp_name)
%         filepath=cd;
%         Fs=60e3;
%         Window_Low= 50;
%         Window_High= 10000;
%         date=exp_name(4:9);
%         date = textscan(date,'%2s');
%         date = date{1};
%         %% Peak Detection
%         outputfile = ['Threshold',num2str(thresh),'_070223_peak_values_',...
%                                         date{1},'_',date{2},'_',date{3}];
%         dirinfo = dir('*Cell*');
%         dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
%         dirinfo(ismember( {dirinfo.name}, {'.', '..'})) = [];  %remove . and ..
%         file_range = (1:length(dirinfo));
%         analysisvals = 1;
%         sample_type = 'Blood';
%         exp_num = round(6000/9);
%         std_threshold = thresh;
%         Spectralon_tail = '_1';
%         FWMH_threshold = 0;
%         intensity_threshold = 0.1;
%         bead_flag = Final_BeadFlag(i);
%         output = SamplePeakDetection_PCA_PN(filepath,outputfile,file_range,Window_Low,...
%             Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
%             Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
%         disp(output)
%         cd(mainpath)
%     end
    header = ['Threshold',num2str(thresh),'_070223'];
    [sens,pur] = ScatvFLRCalculator_DeepPeak_ThresholdTest(mainpath,header);
    performance(j,2:3) = [sens,pur];
    cd(mainpath)
end
header = 'NEW_peak_values';
[sens,pur] = ScatvFLRCalculator_DeepPeak_ThresholdTest(mainpath,header);
performance(6,2:3) = [sens,pur];
performance_final = sortrows(performance,1);

%% Ploting of data
figure;
performance_final(5,:) = [];
x = performance_final(:,1);
sens_final = performance_final(:,2).*100;
pur_final = performance_final(:,3).*100;
yyaxis('left')
plot(x,sens_final,'.','MarkerSize',23)
ylabel('Sensitivity Performance(%)')
yyaxis("right")
plot(x,pur_final,'.','MarkerSize',23)
legend('Sensitivity','Purity','NumColumns',2,'Location','southoutside')
xlabel('Threshold')
ylabel('Purity Performance (%)')
xticks([6:2:14])
xlim([5,15])
boldify
