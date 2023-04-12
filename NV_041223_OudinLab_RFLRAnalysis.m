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
counts = struct;
for i=1:length(T)
    disp(['Analyzing Day ',num2str(i),' of ' num2str(length(T))])
    exp_name=T(i).name;
    cd(exp_name)
    outpath=cd;
    subfolders = dir();
    subfolders(~[subfolders.isdir]) = [];  %remove non-directories
    subfolders(ismember( {subfolders.name}, {'.', '..','Local Plots'})) = [];  %remove . and ..
    for j = 1:length(subfolders)
        data_path = subfolders(j).name;
        cd(data_path)
        data_folders = dir();
        data_folders(~[data_folders.isdir]) = [];  %remove non-directories
        data_folders(ismember( {data_folders.name}, {'.', '..','Local Plots'})) = [];  %remove . and ..
        data_name = data_folders(1).name;
        rFlr_data = dir('Corrected*');
        rFlr_data = rFlr_data(2);
        % Check
        if contains(rFlr_data.name,'RFLR') == 0
           disp('Error')
           break 
        end
        load(rFlr_data.name,'peak_values');
        counts.(data_name) = size(peak_values,1);
        cd(outpath)
    end
    cd(mainpath)
end
%% Cancer Mice
T = dir('*Cancer*');
for i=1:length(T)
    disp(['Analyzing Day ',num2str(i),' of ' num2str(length(T))])
    exp_name=T(i).name;
    cd(exp_name)
    outpath=cd;
    data_folders = dir();
    data_folders(~[data_folders.isdir]) = [];  %remove non-directories
    data_folders(ismember( {data_folders.name}, {'.', '..','Local Plots'})) = [];  %remove . and ..
    data_name = data_folders(1).name;
    rFlr_data = dir('Corrected*');
    rFlr_data = rFlr_data(2);
    % Check
    if contains(rFlr_data.name,'RFLR') == 0
       disp('Error')
       break 
    end
    load(rFlr_data.name,'peak_values');
    if isfield(counts,data_name)
        temp_val = counts.(data_name);
        new_val = [temp_val;size(peak_values,1)];
    else
        new_val = size(peak_values,1);
    end
    counts.(data_name) = new_val;
    cd(mainpath)
end
%% Analysis
cd('T:\Nilay\IVFC\Acquired Data\SydneyMice')
T = readtable("MetsNotebook.xlsx");
Obese = T(strcmp(T.DietType, 'Obese'),:);
Lean = T(strcmp(T.DietType, 'Lean'),:);
timePts = unique(T.TimePoint,'stable');
figure;
meanVal = zeros(length(timePts),1);
for i = 1:length(timePts)
    int = find(strcmp(Obese.TimePoint,timePts{i}) == 1);
    plot(repmat(i,length(int),1),Obese.Average(int),'r.','MarkerSize',23)
    hold on
    meanVal(i) = mean(Obese.Average(int));
end
plot(1:length(timePts),meanVal,'r--','LineWidth',3)

meanVal = zeros(length(timePts),1);
for i = 1:length(timePts)
    int = find(strcmp(Lean.TimePoint,timePts{i}) == 1);
    plot(repmat(i,length(int),1),Lean.Average(int),'b.','MarkerSize',23)
    hold on
    meanVal(i) = mean(Lean.Average(int));
end
plot(1:length(timePts),meanVal,'b--','LineWidth',3)

xticks(1:1:4)
xticklabels(timePts)
legend('Obese Counts','','','','Mean Obese Counts','Lean Counts','',...
    '','','Mean Lean Counts','Location','southoutside','NumColumns',2)
boldify

fields = fieldnames(counts);
Cancerfields = fields(5:end);
Average = zeros(23,1);
rownum = 1;
order = [1,2,3,4,6,5,8,7];
for n = 1:length(Cancerfields)
    i = order(n);
    data = counts.(Cancerfields{i});
    Average(rownum:length(data)+rownum-1,1) = data;
    rownum = rownum + length(data);
    clear data
end
Mouse = T.Mouse;
TimePoint = T.TimePoint;
DietType = T.DietType;

RFLRTable = table(Mouse,TimePoint,Average,DietType);
ObeseRFLR = RFLRTable(strcmp(RFLRTable.DietType, 'Obese'),:);
LeanRFLR = RFLRTable(strcmp(RFLRTable.DietType, 'Lean'),:);

%% Organize data
Obese_final = [Obese;ObeseRFLR];
Obese_final.Source = [repmat({'MetsCount'},12,1);repmat({'RFLR'},12,1)];

colors=[0.49,0.18,0.56;...
        0.30,0.75,0.93;...
        0.64,0.08,0.18;...
        0.93,0.69,0.13;...
        0.47,0.67,0.19];
figure;
for i = 1:length(timePts)
    int = find(strcmp(Obese.TimePoint,timePts{i}) == 1);
    int_DP = find(strcmp(ObeseRFLR.TimePoint,timePts{i}) == 1);
    plot(ObeseRFLR.Average(int_DP),Obese.Average(int),'.',...
            'Color',colors(i,:),'MarkerSize',23)
    hold on
end
legend('t = 4 weeks','t = 8 weeks','t = 12 weeks','t = 16 weeks',...
        'Location','southoutside','numColumns',4)
title('Obese Mice')
boldify
xlabel('RFLR Detected Events')
ylabel('Counted Mets')

Lean_final = [Lean;LeanRFLR];
Lean_final.Source = [repmat({'MetsCount'},11,1);repmat({'RFLR'},11,1)];

colors=[0.49,0.18,0.56;...
        0.30,0.75,0.93;...
        0.64,0.08,0.18;...
        0.93,0.69,0.13;...
        0.47,0.67,0.19];
figure;
for i = 1:length(timePts)
    int = find(strcmp(Lean.TimePoint,timePts{i}) == 1);
    int_DP = find(strcmp(LeanRFLR.TimePoint,timePts{i}) == 1);
    plot(LeanRFLR.Average(int_DP),Lean.Average(int),'.',...
            'Color',colors(i,:),'MarkerSize',23)
    hold on
end
legend('t = 4 weeks','t = 8 weeks','t = 12 weeks','t = 16 weeks',...
        'Location','southoutside','numColumns',4)
title('Lean Mice')
boldify
xlabel('RFLR Detected Events')
ylabel('Counted Mets')

figure;
meanVal = zeros(length(timePts),1);
for i = 1:length(timePts)
    int = find(strcmp(ObeseRFLR.TimePoint,timePts{i}) == 1);
    plot(repmat(i,length(int),1),ObeseRFLR.Average(int),'r.','MarkerSize',23)
    hold on
    meanVal(i) = mean(ObeseRFLR.Average(int));
end
plot(1:length(timePts),meanVal,'r--','LineWidth',3)

meanVal = zeros(length(timePts),1);
for i = 1:length(timePts)
    int = find(strcmp(LeanRFLR.TimePoint,timePts{i}) == 1);
    plot(repmat(i,length(int),1),LeanRFLR.Average(int),'b.','MarkerSize',23)
    hold on
    meanVal(i) = mean(LeanRFLR.Average(int));
end
plot(1:length(timePts),meanVal,'b--','LineWidth',3)

xticks(1:1:4)
xticklabels(timePts)
legend('Obese Counts','','','','Mean Obese Counts','Lean Counts','',...
    '','','Mean Lean Counts','Location','southoutside','NumColumns',2)
boldify

