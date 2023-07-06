%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  file: NV_032123_SensitivityStudy.m
%  authors: Nilay Vora
%  created: 03/21/23
%  modified: 03/28/23
%  purpose: This file will be used to create datasets 
%  for deep learning with a focus on specific number 
%  of clusters and FPs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
clear
clc
%% Load Data
cd('T:\Nilay\IVFC\Acquired Data\Blood Cell Data\DeepPeak2023\')
load("NV_060823_FormattedDataSet.mat","FP_Peak_ranges");
load("NV_060823_FormattedDataSet_locs.mat","Loc_Peak_ranges");
%% Main Loop
fields = fieldnames(FP_Peak_ranges);
CTCC_Count = zeros(length(fields),1);
for i = 1:length(fields)
CTCC_Count(i,1) = length(find(FP_Peak_ranges.(fields{i})(:,397)==0));
end
low_data = [];
totalPeaks = FP_Peak_ranges;
thresholds = [5,10,30,50,100];
for t = 1:length(thresholds)
    disp(['Applying threshold ',num2str(t), ' of ',...
                                            num2str(length(thresholds))])
    FP_Peak_ranges = struct();
    Loc_Peak_ranges = struct();
    threshold = thresholds(t);
    for i = 1:length(fields)
        disp([9 'Formatting field: ',fields{i}])
        data = totalPeaks.(fields{i});
        % Find the threshold CTCCs:
        allCTCCs = find(data(:,397)==0);
        if size(allCTCCs,1) > threshold
            thresholdCTCCs = allCTCCs(threshold);
        else
            thresholdCTCCs = allCTCCs(end);
            low_data{end+1} = {fields{i},thresholds(t)};
        end
        finaldata = data(1:thresholdCTCCs,:);
        finallocs = data(1:thresholdCTCCs,:);
        FP_Peak_ranges.(fields{i}) = finaldata;
        Loc_Peak_ranges.(fields{i}) = finallocs;
    end
    save_name = ['NV_060823_FormattedDataSet_thresh',num2str(threshold,'%03d')];
    save_name2 = ['NV_060823_FormattedDataSet_thresh',num2str(threshold,'%03d'),'_locs'];
    save(save_name,"FP_Peak_ranges")
    save(save_name2,"Loc_Peak_ranges")
end



