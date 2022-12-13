%% NV_072222_ObesityStudy %%
%% Information about algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Low SNR in the Green FLR channel makes it difficult to analyze. A better
% way to analyze this data is either follow Joe's procedure for mahalnobis
% distance OR using my previously validated ML model. Here we will attempt
% the use of the ML model. In an alternative code, we will use Joe's
% Mahalnobis distance measurment to statistically determine the location of
% peaks.
%
% We plan to feed all the scattering data from all detected GFLR peaks, in
% essence, we should be cleaning all noise peaks from our data using this
% model... I hope
%
% Date Created: 07/28/2022
% Author: Nilay Vora
% Code Dependenacies: EBT Model saved in a folder (v large, monitor
% computer RAM to run!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
clear
clc
close all
mainpath = 'T:\Nilay\IVFC\Acquired Data\SydneyMice';
cd(mainpath)
%% Load Data and Format it for use in ML Algorithm
%% Load in all Lean Healthy Samples
dirinfo=dir("*HealthyControl");
k = 1;
r_store=cell(1,1);
label_store=cell(1,1);
peak_value_store=cell(1,1);
for i=1:length(dirinfo)
    disp(['Running Sample# ', num2str(i), ' of ' num2str(length(dirinfo))])
    cd([dirinfo(i).folder,'\',dirinfo(i).name])
    cd(dir("*Lean").name)
    mainfolder = cd;
    % Get list of of subfolders
    subdirinfo = dir();
    subdirinfo(~[subdirinfo.isdir]) = [];
    subdirinfo(ismember({subdirinfo.name}, {'.', '..'})) = [];  %remove . and ..
    % Generate peakvalue file name and folder names
    for j = 1:size(subdirinfo,1)
        disp(['     Running Timepoint# ', num2str(j), ' of ' num2str(size(subdirinfo,1))])
        peakfolder = subdirinfo(j).name;
        cd(peakfolder)
        date = peakfolder(4:9);
        date = [date(1:2),'_',date(3:4),'_',date(5:6)];
        fname = ['T',num2str(j),'-NEW_peak_values_',date,'_NoScatAll.mat'];
        [range_data,label,peak_values] = MLformating(fname,peakfolder);
        lbead=find(label==3); %Remove any beads from data
        label(lbead)=[];
        range_data(lbead,:)=[];
        peak_values(lbead,:)=[];
        peak_values(:,22)=k;
        X=[range_data,label];
        label_store{k,1}=X(:,109);
        r_store{k,1}=X(:,1:108);
        peak_value_store{k,1}=peak_values;
        k = k + 1;
        cd(mainfolder)
    end
end
LeanHealthyMice_data = struct('Timepoint0',r_store);
LeanHealthyMice_labels = struct('Timepoint0',label_store);
LeanHealthyMice_peaks = struct('Timepoint0',peak_value_store);
%% Load in all Obese Healthy Samples
cd(mainpath)
dirinfo=dir("*HealthyControl");
k = 1;
r_store=cell(1,1);
label_store=cell(1,1);
peak_value_store=cell(1,1);
for i=1:length(dirinfo)
    disp(['Running Sample# ', num2str(i), ' of ' num2str(length(dirinfo))])
    cd([dirinfo(i).folder,'\',dirinfo(i).name])
    cd(dir("*Obese").name)
    mainfolder = cd;
    % Get list of of subfolders
    subdirinfo = dir();
    subdirinfo(~[subdirinfo.isdir]) = [];
    subdirinfo(ismember({subdirinfo.name}, {'.', '..'})) = [];  %remove . and ..
    % Generate peakvalue file name and folder names
    for j = 1:size(subdirinfo,1)
        disp(['     Running Timepoint# ', num2str(j), ' of ' num2str(size(subdirinfo,1))])
        peakfolder = subdirinfo(j).name;
        cd(peakfolder)
        date = peakfolder(4:9);
        date = [date(1:2),'_',date(3:4),'_',date(5:6)];
        fname = ['T',num2str(j),'-NEW_peak_values_',date,'_NoScatAll.mat'];
        [range_data,label,peak_values] = MLformating(fname,peakfolder);
        lbead=find(label==3); %Remove any beads from data
        label(lbead)=[];
        range_data(lbead,:)=[];
        peak_values(lbead,:)=[];
        peak_values(:,22)=k;
        X=[range_data,label];
        label_store{k,1}=X(:,109);
        r_store{k,1}=X(:,1:108);
        peak_value_store{k,1}=peak_values;
        k = k + 1;
        cd(mainfolder)
    end
end
ObeseHealthyMice_data = struct('Timepoint0',r_store);
ObeseHealthyMice_labels = struct('Timepoint0',label_store);
ObeseHealthyMice_peaks = struct('Timepoint0',peak_value_store);
%% Experiments
cd(mainpath)
Obese_dates = {'041322';'051122';'060822';'070722'};
Lean_dates = {'041422';'051222';'060822';'070622'};
%% Lean Data
LeanCancerMice(4)=struct();
for i=1:length(Lean_dates)%Loop through Lean Mice Days
    label_store = [];
    r_store = [];
    peak_value_store = [];
    l = 1;
    cd(mainpath)
    disp(['Currently Processing Day # ',num2str(i),' of ',num2str(length(Lean_dates))])
    foldername = ['NV_',Lean_dates{i},'_LeanCancer'];
    LeanMice = dir([foldername,'*']);
    for j=1:length(LeanMice) %Loop through sample per day
        disp(['  Currently Processing Sample # ',num2str(j),' of ',num2str(length(LeanMice))])
        cd([LeanMice(j).folder,'\',LeanMice(j).name])
        mainfolder = cd;
        subdirinfo = dir();
        subdirinfo(~[subdirinfo.isdir]) = [];
        subdirinfo(ismember({subdirinfo.name}, {'.', '..'})) = [];  %remove . and ..
        % Generate peakvalue file name and folder names
        for k = 1:size(subdirinfo,1)
            disp(['    Running Timepoint# ', num2str(k), ' of ' num2str(size(subdirinfo,1))])
            peakfolder = subdirinfo(k).name;
            cd(peakfolder)
            date = peakfolder(4:9);
            date = [date(1:2),'_',date(3:4),'_',date(5:6)];
            fname = ['T',num2str(k),'-NEW_peak_values_',date,'_NoScatAll.mat'];
            [range_data,label,peak_values] = MLformating(fname,peakfolder);
            lbead=find(label==3); %Remove any beads from data
            label(lbead)=[];
            range_data(lbead,:)=[];
            peak_values(lbead,:)=[];
            peak_values(:,22)=l;
            X=[range_data,label];
            label_store{l,1}=X(:,109);
            r_store{l,1}=X(:,1:108);
            peak_value_store{l,1}=peak_values;
            l = l + 1;
            cd(mainfolder)
        end
    end
    LeanCancerMice(i).data = r_store;
    LeanCancerMice(i).labels = label_store;
    LeanCancerMice(i).peaks = peak_value_store;
end
%% Obese Data
ObeseCancerMice(4)=struct();
for i=1:length(Obese_dates)%Loop through Obese Mice Days
    label_store = [];
    r_store = [];
    peak_value_store = [];
    l = 1;
    cd(mainpath)
    disp(['Currently Processing Day # ',num2str(i),' of ',num2str(length(Lean_dates))])
    foldername = ['NV_',Obese_dates{i},'_ObeseCancer'];
    ObeseMice = dir([foldername,'*']);
    for j=1:length(ObeseMice) %Loop through sample per day
        disp(['  Currently Processing Sample # ',num2str(j),' of ',num2str(length(ObeseMice))])
        cd([ObeseMice(j).folder,'\',ObeseMice(j).name])
        mainfolder = cd;
        subdirinfo = dir();
        subdirinfo(~[subdirinfo.isdir]) = [];
        subdirinfo(ismember({subdirinfo.name}, {'.', '..'})) = [];  %remove . and ..
        % Generate peakvalue file name and folder names
        for k = 1:size(subdirinfo,1)
            disp(['    Running Timepoint# ', num2str(k), ' of ' num2str(size(subdirinfo,1))])
            peakfolder = subdirinfo(k).name;
            cd(peakfolder)
            date = peakfolder(4:9);
            date = [date(1:2),'_',date(3:4),'_',date(5:6)];
            fname = ['T',num2str(k),'-NEW_peak_values_',date,'_NoScatAll.mat'];
            [range_data,label,peak_values] = MLformating(fname,peakfolder);
            lbead=find(label==3); %Remove any beads from data
            label(lbead)=[];
            range_data(lbead,:)=[];
            peak_values(lbead,:)=[];
            peak_values(:,22)=l;
            X=[range_data,label];
            label_store{l,1}=X(:,109);
            r_store{l,1}=X(:,1:108);
            peak_value_store{l,1}=peak_values;
            l = l + 1;
            cd(mainfolder)
        end
    end
    ObeseCancerMice(i).data = r_store;
    ObeseCancerMice(i).labels = label_store;
    ObeseCancerMice(i).peaks = peak_value_store;
end
%% Load in ML Model #89
cd('T:\Nilay\IVFC\Acquired Data\Blood Cell Data\CNNData')
load('GBEnsembleModel.mat')
GBEnsemble = GBEnsemble_store{1};
clear GBEnsemble_store
%% Model Eval Script (Lean Healthy Mice)
count = 0;
healthy_lean_loc = [];
for i = 1:length(LeanHealthyMice_data)
    data = LeanHealthyMice_data(i).Timepoint0;
    testing_locs = LeanHealthyMice_peaks(i).Timepoint0;
    X = data(1:end,[1:81]);
    disp(['Working on Subfile ',num2str(i),' of ', num2str(length(LeanHealthyMice_data))])
    for iii=1:50
        disp(['     Testing: Running Model # ',num2str(iii),' of 50'])
        X_test=X(:,1:81);
        yfit = GBEnsemble(iii).predictFcn(X_test);
        X=[X_test(yfit==1,:)];
        testing_locs=testing_locs(yfit==1,:);
    end
    healthy_lean_loc = [healthy_lean_loc;testing_locs];
    count = count + size(X,1);
end
Healthy_Lean_count = count;

%% Model Eval Script (Obese Healthy Mice)
count = 0;
healthy_obese_loc = [];
for i = 1:length(ObeseHealthyMice_data)
    data = ObeseHealthyMice_data(i).Timepoint0;
    testing_locs = ObeseHealthyMice_peaks(i).Timepoint0;
    X = data(1:end,[1:81]);
    disp(['Working on Subfile ',num2str(i),' of ', num2str(length(LeanHealthyMice_data))])
    for iii=1:50
        disp(['     Testing: Running Model # ',num2str(iii),' of 50'])
        X_test=X(:,1:81);
        yfit = GBEnsemble(iii).predictFcn(X_test);
        X=[X_test(yfit==1,:)];
        testing_locs=testing_locs(yfit==1,:);
    end
    healthy_obese_loc = [healthy_obese_loc;testing_locs];
    count = count + size(X,1);
end
Healthy_Obese_count = count;
%% Model Eval Script (Lean Cancer Mice)
Cancer_Lean_tpt_loc = cell(1,4);
Cancer_Lean_count = zeros(1,4);
% Loop through timepts
for i = 1:length(LeanCancerMice)
    data = LeanCancerMice(i).data;
    testing_locs = LeanCancerMice(i).peaks;
    count = 0;
    Cancer_Lean_loc = [];
    % Loop through data per timept
    for ii = 1:length(data)
        t_data = data{ii};
        X = t_data(1:end,1:81);
        t_testing_locs = testing_locs{ii};
        disp(['     Working on Subfile ',num2str(ii),' of ', num2str(length(data))])
        %loop through model
        for iii=1:50
            disp(['          Testing: Running Model # ',num2str(iii),' of 50'])
            X_test=X(:,1:81);
            yfit = GBEnsemble(iii).predictFcn(X_test);
            X=[X_test(yfit==1,:)];
            t_testing_locs=t_testing_locs(yfit==1,:);
        end
        count = count + size(X,1);
        Cancer_Lean_loc = [Cancer_Lean_loc;t_testing_locs];
    end
    Cancer_Lean_tpt_loc{i}= Cancer_Lean_loc;
    Cancer_Lean_count(i) = count;
end
%% Model Eval Script (Obese Cancer Mice)
Cancer_Obese_tpt_loc = cell(1,4);
Cancer_Obese_count = zeros(1,4);
% Loop through timepts
for i = 1:length(ObeseCancerMice)
    data = ObeseCancerMice(i).data;
    testing_locs = ObeseCancerMice(i).peaks;
    count = 0;
    Cancer_Obese_loc = [];
    % Loop through data per timept
    for ii = 1:length(data)
        t_data = data{ii};
        t_testing_locs = testing_locs{ii};
        X = t_data(1:end,1:81);
        disp(['     Working on Subfile ',num2str(ii),' of ', num2str(length(data))])
        %loop through model
        for iii=1:50
            disp(['          Testing: Running Model # ',num2str(iii),' of 50'])
            X_test=X(:,1:81);
            yfit = GBEnsemble(iii).predictFcn(X_test);
            X=[X_test(yfit==1,:)];
            t_testing_locs=t_testing_locs(yfit==1,:);
        end
        count = count + size(X,1);
        Cancer_Obese_loc = [Cancer_Obese_loc;t_testing_locs];
    end
    Cancer_Obese_tpt_loc{i} = Cancer_Obese_loc;
    Cancer_Obese_count(i) = count;
end
%% Seperate Samples per Mouse in Healthy
pv = healthy_lean_loc;
Mouse_count = length(find(pv(1:end-1,21)-pv(2:end,21)>0))+1;
Mouse_idx = [1,find(pv(1:end-1,21)-pv(2:end,21)>0)+1];
Healthy_lean_locs_sorted = cell(1,Mouse_count);
Healthy_lean_count_sorted = zeros(1,Mouse_count);
for i = 1:Mouse_count
    if i+1 > Mouse_count
        range = Mouse_idx(i):size(pv,1);
    else
        range = Mouse_idx(i):(Mouse_idx(i+1)-1);
    end
    pv_temp = pv(range,:);
    pv_temp = pv_temp(pv_temp(:,19)>9,:);
    Healthy_lean_locs_sorted{1,i} = pv_temp;
    Healthy_lean_count_sorted(1,i) = size(pv_temp,1);
end

pv = healthy_obese_loc;
Mouse_count = length(find(pv(1:end-1,21)-pv(2:end,21)>0))+1;
Mouse_idx = [1,find(pv(1:end-1,21)-pv(2:end,21)>0)+1];
Healthy_obese_locs_sorted = cell(1,Mouse_count);
Healthy_obese_count_sorted = zeros(1,Mouse_count);
for i = 1:Mouse_count
    if i+1 > Mouse_count
        range = Mouse_idx(i):size(pv,1);
    else
        range = Mouse_idx(i):(Mouse_idx(i+1)-1);
    end
    pv_temp = pv(range,:);
    pv_temp = pv_temp(pv_temp(:,19)>9,:);
    Healthy_obese_locs_sorted{1,i} = pv_temp;
    Healthy_obese_count_sorted(1,i) = size(pv_temp,1);
end
%% Seperate Samples per Mouse in Cancer
Cancer_lean_locs_sorted = cell(length(Cancer_Lean_tpt_loc),3);
Cancer_lean_count_sorted = zeros(length(Cancer_Lean_count),3);
for j = 1:length(Cancer_Lean_tpt_loc)
    pv = Cancer_Lean_tpt_loc{j};
    Mouse_count = length(find(pv(1:end-1,21)-pv(2:end,21)>0))+1;
    Mouse_idx = [1;find(pv(1:end-1,21)-pv(2:end,21)>0)+1];
    for i = 1:Mouse_count
        if i+1 > Mouse_count
            range = Mouse_idx(i):size(pv,1);
        else
            range = Mouse_idx(i):(Mouse_idx(i+1)-1);
        end
        pv_temp = pv(range,:);
        pv_temp = pv_temp(pv_temp(:,19)>9,:);
        Cancer_lean_locs_sorted{j,i} = pv_temp;
        Cancer_lean_count_sorted(j,i) = size(pv_temp,1);
    end
end

Cancer_obese_locs_sorted = cell(length(Cancer_Obese_tpt_loc),3);
Cancer_obese_count_sorted = zeros(length(Cancer_Obese_count),3);
for j = 1:length(Cancer_Obese_tpt_loc)
    pv = Cancer_Obese_tpt_loc{j};
    Mouse_count = length(find(pv(1:end-1,21)-pv(2:end,21)>0))+1;
    Mouse_idx = [1;find(pv(1:end-1,21)-pv(2:end,21)>0)+1];
    for i = 1:Mouse_count
        if i+1 > Mouse_count
            range = Mouse_idx(i):size(pv,1);
        else
            range = Mouse_idx(i):(Mouse_idx(i+1)-1);
        end
        pv_temp = pv(range,:);
        pv_temp = pv_temp(pv_temp(:,19)>9,:);
        Cancer_obese_locs_sorted{j,i} = pv_temp;
        Cancer_obese_count_sorted(j,i) = size(pv_temp,1);
    end
end

%% Time Calculations for Healthy Lean
cd(mainpath)
dirinfo=dir("*HealthyControl");
Healthy_lean_time = zeros(1,length(dirinfo));
for i = 1:length(dirinfo)
    disp(['Running Mouse# ', num2str(i), ' of ' num2str(length(dirinfo))])
    cd([dirinfo(i).folder,'\',dirinfo(i).name])
    cd(dir("*Lean").name)
    mainfolder = cd;
    subdirinfo = dir();
    subdirinfo(~[subdirinfo.isdir]) = [];
    subdirinfo(ismember({subdirinfo.name}, {'.', '..'})) = [];  %remove . and ..
    % Generate peakvalue file name and folder names
    pts_mouse = 0;
    for j = 1:length(subdirinfo)
        disp(['     Running Sample# ', num2str(j), ' of ' num2str(size(subdirinfo,1))])
        peakfolder = subdirinfo(j).name;
        cd(peakfolder)
        files = dir('*raw.mat');
        for k = 1 : length(files)
            disp(['          Running Chunk# ', num2str(k), ' of ' num2str(size(files,1))])
            file = files(k).name;
            load(file,'M')
            npts = size(M,1);
            pts_mouse = pts_mouse + npts;
        end
        cd(mainfolder)
    end
    Healthy_lean_time(1,i) = pts_mouse./60e3/60; % Saves time in Minutes
end
%% Time Calculations for Healthy Obese
cd(mainpath)
dirinfo=dir("*HealthyControl");
Healthy_obese_time = zeros(1,length(dirinfo));
for i = 1:length(dirinfo)
    disp(['Running Mouse# ', num2str(i), ' of ' num2str(length(dirinfo))])
    cd([dirinfo(i).folder,'\',dirinfo(i).name])
    cd(dir("*Obese").name)
    mainfolder = cd;
    subdirinfo = dir();
    subdirinfo(~[subdirinfo.isdir]) = [];
    subdirinfo(ismember({subdirinfo.name}, {'.', '..'})) = [];  %remove . and ..
    % Generate peakvalue file name and folder names
    pts_mouse = 0;
    for j = 1:length(subdirinfo)
        disp(['     Running Sample# ', num2str(j), ' of ' num2str(size(subdirinfo,1))])
        peakfolder = subdirinfo(j).name;
        cd(peakfolder)
        files = dir('*raw.mat');
        for k = 1 : length(files)
            disp(['          Running Chunk# ', num2str(k), ' of ' num2str(size(files,1))])
            file = files(k).name;
            load(file,'M')
            npts = size(M,1);
            pts_mouse = pts_mouse + npts;
        end
        cd(mainfolder)
    end
    Healthy_obese_time(1,i) = pts_mouse./60e3/60; % Saves time in Minutes
end

%% Calculating Time for all cancer samples
Cancer_lean_time = zeros(4,3);
for i=1:length(Lean_dates)%Loop through Lean Mice Days
    cd(mainpath)
    disp(['Currently Processing Day # ',num2str(i),' of ',num2str(length(Lean_dates))])
    foldername = ['NV_',Lean_dates{i},'_LeanCancer'];
    LeanMice = dir([foldername,'*']);
    for j=1:length(LeanMice) %Loop through sample per day
        disp(['     Running Mouse# ', num2str(i), ' of ' num2str(length(LeanMice))])
        cd([LeanMice(j).folder,'\',LeanMice(j).name])
        mainfolder = cd;
        subdirinfo = dir();
        subdirinfo(~[subdirinfo.isdir]) = [];
        subdirinfo(ismember({subdirinfo.name}, {'.', '..'})) = [];  %remove . and ..
        pts_mouse = 0;
        for l = 1:length(subdirinfo)
            disp(['          Running Sample# ', num2str(l), ' of ' num2str(size(subdirinfo,1))])
            peakfolder = subdirinfo(l).name;
            cd(peakfolder)
            files = dir('*raw.mat');
            for k = 1 : length(files)
                disp(['               Running Chunk# ', num2str(k), ' of ' num2str(size(files,1))])
                file = files(k).name;
                load(file,'M')
                npts = size(M,1);
                pts_mouse = pts_mouse + npts;
            end
            cd(mainfolder)
        end
        Cancer_lean_time(i,j) = pts_mouse./60e3/60; % Saves time in Minutes
    end
end

Cancer_obese_time = zeros(4,3);
for i=1:length(Obese_dates)%Loop through Obese Mice Days
    cd(mainpath)
    disp(['Currently Processing Day # ',num2str(i),' of ',num2str(length(Obese_dates))])
    foldername = ['NV_',Obese_dates{i},'_ObeseCancer'];
    ObeseMice = dir([foldername,'*']);
    for j=1:length(ObeseMice) %Loop through sample per day
        disp(['     Running Mouse# ', num2str(i), ' of ' num2str(length(ObeseMice))])
        cd([ObeseMice(j).folder,'\',ObeseMice(j).name])
        mainfolder = cd;
        subdirinfo = dir();
        subdirinfo(~[subdirinfo.isdir]) = [];
        subdirinfo(ismember({subdirinfo.name}, {'.', '..'})) = [];  %remove . and ..
        pts_mouse = 0;
        for l = 1:length(subdirinfo)
            disp(['          Running Sample# ', num2str(l), ' of ' num2str(size(subdirinfo,1))])
            peakfolder = subdirinfo(l).name;
            cd(peakfolder)
            files = dir('*raw.mat');
            for k = 1 : length(files)
                disp(['               Running Chunk# ', num2str(k), ' of ' num2str(size(files,1))])
                file = files(k).name;
                load(file,'M')
                npts = size(M,1);
                pts_mouse = pts_mouse + npts;
            end
            cd(mainfolder)
        end
        Cancer_obese_time(i,j) = pts_mouse./60e3/60; % Saves time in Minutes
    end
end
%% Normalization of data
Lean_mat = [[Healthy_lean_count_sorted,NaN];Cancer_lean_count_sorted];
Lean_time_mat = [[Healthy_lean_time,0];Cancer_lean_time];
Obese_mat = [[Healthy_obese_count_sorted,NaN];Cancer_obese_count_sorted];
Obese_time_mat = [[Healthy_obese_time,0];Cancer_obese_time];

Lean_Norm = Lean_mat./Lean_time_mat;
Obese_Norm = Obese_mat./Obese_time_mat;
%% Unnormalized counts
colors = linspecer(2);
lean_norm = Lean_mat;
obese_norm = Obese_mat;
figure; plot(1:1:5,mean(lean_norm,2,'omitnan'),'.','Color',colors(1,:),'MarkerSize',30)
hold on
plot(1:1:5,mean(obese_norm,2,'omitnan'),'.','Color',colors(2,:),'MarkerSize',30)
title('Average Number of Events Per Timepoint')
xticks(1:1:5)
xticklabels({'4 Weeks', '8 Weeks', '12 Weeks', '16 Weeks', '20 Weeks'})
xlabel('Timepoint')
ylabel('Events')
boldify
errorbar(1:1:5,mean(lean_norm,2,'omitnan'),std(lean_norm,[],2,'omitnan'),'LineWidth',3,'CapSize',10,'Color',colors(1,:))
errorbar(1:1:5,mean(obese_norm,2,'omitnan'),std(obese_norm,[],2,'omitnan'),'LineWidth',2,'CapSize',10,'Color',colors(2,:))
legend('Lean Cancer Mice','Obese Cancer Mice','numcolumns',2,'Location','Southoutside')
xlim([0.5,5.5])
%% Normalized counts
colors = linspecer(2);
lean_norm = Lean_Norm;
obese_norm = Obese_Norm;
figure; plot(1:1:5,mean(lean_norm,2,'omitnan'),'.','Color',colors(1,:),'MarkerSize',30)
hold on
plot(1:1:5,mean(obese_norm,2,'omitnan'),'.','Color',colors(2,:),'MarkerSize',30)
title('Average Number of Events Per Timepoint')
xticks(1:1:5)
xticklabels({'4 Weeks','8 Weeks', '12 Weeks', '16 Weeks', '20 Weeks'})
xlabel('Timepoint')
ylabel('Events/Minute (min^{-1})')
boldify
errorbar(1:1:5,mean(lean_norm,2,'omitnan'),std(lean_norm,[],2,'omitnan'),'LineWidth',3,'CapSize',10,'Color',colors(1,:))
errorbar(1:1:5,mean(obese_norm,2,'omitnan'),std(obese_norm,[],2,'omitnan'),'LineWidth',2,'CapSize',10,'Color',colors(2,:))
legend('Lean Cancer Mice','Obese Cancer Mice','numcolumns',2,'Location','Southoutside')
xlim([0.5,5.5])

%% AHHHHHH MAHAL TIME CAUSE NOTHING ELSE WORKED
% Healthy sample peaks Lean
PV_Lean_WT = PvFormat(healthy_lean_loc);
PV_Lean_WT1 = PvFormat(Healthy_lean_locs_sorted{1});
PV_Lean_WT2 = PvFormat(Healthy_lean_locs_sorted{2});

mahal_Lean_WT1=mahal(PV_Lean_WT1,PV_Lean_WT);
WT1=repmat({'Lean_WT1'},[size(PV_Lean_WT1,1),1]);
mahal_Lean_WT2=mahal(PV_Lean_WT2,PV_Lean_WT);
WT2=repmat({'Lean_WT2'},[size(PV_Lean_WT2,1),1]);

M_dist=[mahal_Lean_WT1;mahal_Lean_WT2];
groups=[WT1;WT2];
threshold=max(M_dist);
outliers=[];
everything=[];
for i = 1:4
    for j = 1:3
        flag = 0;
        if isempty(Cancer_lean_locs_sorted{i,j})
            flag = 1;     
        else
            peak_values = Cancer_lean_locs_sorted{i,j};
            peakdata = PvFormat(peak_values);
            mahal_add=mahal(peakdata, PV_Lean_WT);
            l1=find(mahal_add>threshold);
            l2=find(mahal_add<=threshold);
            if ~isempty(l1)
                outliers=[outliers;[peak_values(l1,:),repmat(i*j,[length(l1),1])]];
            end
            everything=[everything;[peak_values(l2,:),repmat(i*j,[length(l2),1])]];
            M_dist=[M_dist;mahal_add];
            groups=[groups;repmat({['T',num2str(i),'Mouse#',num2str(j)]},[size(peakdata,1),1])];         
        end
    end
end

[Out_mat_Lean,In_mat_Lean,Sig_mat_Lean]=OutlierAnalysis(outliers,everything);

% Plotting Code
figure;
boxplot(M_dist,groups,'Colors','k','PlotStyle','Compact','Symbol','.','OutlierSize',10,'Jitter',0.25,'LabelOrientation','horizontal')
ax=gca;
ax.YAxis.Scale ="log";
threshold=max([mahal_Lean_WT1;mahal_Lean_WT2]);
yline(threshold,'b','LineWidth',3)
title('Two WT Lean Mice (Clusters Only)')
xlabel('Sample')
ylabel('Mahalanobis Distance')
boldify

% Healthy sample peaks obese
PV_Obese_WT = PvFormat(healthy_obese_loc);
PV_Obese_WT1 = PvFormat(Healthy_obese_locs_sorted{1});
PV_Obese_WT2 = PvFormat(Healthy_obese_locs_sorted{2});

mahal_Obese_WT1=mahal(PV_Obese_WT1,PV_Obese_WT);
WT1=repmat({'Obese_WT1'},[size(PV_Obese_WT1,1),1]);
mahal_Obese_WT2=mahal(PV_Obese_WT2,PV_Obese_WT);
WT2=repmat({'Obese_WT2'},[size(PV_Obese_WT2,1),1]);

M_dist=[mahal_Obese_WT1;mahal_Obese_WT2];
groups=[WT1;WT2];
threshold=max(M_dist);
outliers=[];
everything=[];
for i = 1:4
    for j = 1:3
        flag = 0;
        if isempty(Cancer_obese_locs_sorted{i,j})
            flag = 1;     
        else
            peak_values = Cancer_obese_locs_sorted{i,j};
            peakdata = PvFormat(peak_values);
            mahal_add=mahal(peakdata, PV_Obese_WT);
            l1=find(mahal_add>threshold);
            l2=find(mahal_add<=threshold);
            if ~isempty(l1)
                outliers=[outliers;[peak_values(l1,:),repmat(i*j,[length(l1),1])]];
            end
            everything=[everything;[peak_values(l2,:),repmat(i*j,[length(l2),1])]];
            M_dist=[M_dist;mahal_add];
            groups=[groups;repmat({['T',num2str(i),'Mouse#',num2str(j)]},[size(peakdata,1),1])];         
        end
    end
end

[Out_mat_Obese,In_mat_Obese,Sig_mat_Obese]=OutlierAnalysis(outliers,everything);

% Plotting Code
figure;
boxplot(M_dist,groups,'Colors','k','PlotStyle','Compact','Symbol','.','OutlierSize',10,'Jitter',0.25,'LabelOrientation','horizontal')
ax=gca;
ax.YAxis.Scale ="log";
threshold=max([mahal_Obese_WT1;mahal_Obese_WT2]);
yline(threshold,'b','LineWidth',3)
title('Two WT Obese Mice (Clusters Only)')
xlabel('Sample')
ylabel('Mahalanobis Distance')
boldify
%% AHHHHHH MAHAL TIME CAUSE NOTHING ELSE WORKED... okay it worked
% Healthy sample peaks Lean
PV_Lean_WT = PvFormat(healthy_lean_loc);
PV_Lean_WT1 = PvFormat(Healthy_lean_locs_sorted{1});
PV_Lean_WT2 = PvFormat(Healthy_lean_locs_sorted{2});

mahal_Lean_WT1=mahal(PV_Lean_WT1,PV_Lean_WT);
WT1=repmat({'Lean_WT1'},[size(PV_Lean_WT1,1),1]);
mahal_Lean_WT2=mahal(PV_Lean_WT2,PV_Lean_WT);
WT2=repmat({'Lean_WT2'},[size(PV_Lean_WT2,1),1]);

M_dist=[mahal_Lean_WT1;mahal_Lean_WT2];
groups=[WT1;WT2];
threshold=max(M_dist);
outliers=[];
everything=[];
for i = 1:4
    peak_values = [];
    for j = 1:3
        flag = 0;
        if isempty(Cancer_lean_locs_sorted{i,j})
            flag = 1;     
        else
            peak_values = [peak_values; Cancer_lean_locs_sorted{i,j}];
        end
    end
    peakdata = PvFormat(peak_values);
    mahal_add=mahal(peakdata, PV_Lean_WT);
    l1=find(mahal_add>threshold);
    l2=find(mahal_add<=threshold);
    if ~isempty(l1)
        outliers=[outliers;[peak_values(l1,:),repmat(i*j,[length(l1),1])]];
    end
    everything=[everything;[peak_values(l2,:),repmat(i*j,[length(l2),1])]];
    M_dist=[M_dist;mahal_add];
    groups=[groups;repmat({['T',num2str(i)]},[size(peakdata,1),1])];         
end

disp(['There are ',num2str(size(outliers,1)),...
    ' outliers and ',num2str(size(everything,1)), 'inliers in the Lean data'])
[Out_mat_Lean_tp,In_mat_Lean_tp,Sig_mat_Lean_tp]=OutlierAnalysis(outliers,everything);

% Plotting Code
figure;
boxplot(M_dist,groups,'Colors','k','PlotStyle','Compact','Symbol','.','OutlierSize',10,'Jitter',0.25,'LabelOrientation','horizontal')
ax=gca;
ax.YAxis.Scale ="log";
threshold=max([mahal_Lean_WT1;mahal_Lean_WT2]);
yline(threshold,'b','LineWidth',3)
title('Two WT Lean Mice (Clusters Only)')
xlabel('Sample')
ylabel('Mahalanobis Distance')
boldify

% Healthy sample peaks obese
PV_Obese_WT = PvFormat(healthy_obese_loc);
PV_Obese_WT1 = PvFormat(Healthy_obese_locs_sorted{1});
PV_Obese_WT2 = PvFormat(Healthy_obese_locs_sorted{2});

mahal_Obese_WT1=mahal(PV_Obese_WT1,PV_Obese_WT);
WT1=repmat({'Obese_WT1'},[size(PV_Obese_WT1,1),1]);
mahal_Obese_WT2=mahal(PV_Obese_WT2,PV_Obese_WT);
WT2=repmat({'Obese_WT2'},[size(PV_Obese_WT2,1),1]);

M_dist=[mahal_Obese_WT1;mahal_Obese_WT2];
groups=[WT1;WT2];
threshold=max(M_dist);
outliers=[];
everything=[];
for i = 1:4
    peak_values = [];
    for j = 1:3
        flag = 0;
        if isempty(Cancer_obese_locs_sorted{i,j})
            flag = 1;     
        else
            peak_values = [peak_values; Cancer_obese_locs_sorted{i,j}];
        end
    end
    peakdata = PvFormat(peak_values);
    mahal_add=mahal(peakdata, PV_Obese_WT);
    l1=find(mahal_add>threshold);
    l2=find(mahal_add<=threshold);
    if ~isempty(l1)
        outliers=[outliers;[peak_values(l1,:),repmat(i*j,[length(l1),1])]];
    end
    everything=[everything;[peak_values(l2,:),repmat(i*j,[length(l2),1])]];
    M_dist=[M_dist;mahal_add];
    groups=[groups;repmat({['T',num2str(i)]},[size(peakdata,1),1])];         
end
disp(['There are ',num2str(size(outliers,1)),...
    ' outliers and ',num2str(size(everything,1)), 'inliers in the Obese data'])
[Out_mat_Obese_tp,In_mat_Obese_tp,Sig_mat_Obese_tp]=OutlierAnalysis(outliers,everything);

% Plotting Code
figure;
boxplot(M_dist,groups,'Colors','k','PlotStyle','Compact','Symbol','.','OutlierSize',10,'Jitter',0.25,'LabelOrientation','horizontal')
ax=gca;
ax.YAxis.Scale ="log";
threshold=max([mahal_Obese_WT1;mahal_Obese_WT2]);
yline(threshold,'b','LineWidth',3)
title('Two WT Obese Mice (Clusters Only)')
xlabel('Sample')
ylabel('Mahalanobis Distance')
boldify
%% Formatting code
function [range_data,label,peak_values]=MLformating(fname,ftype)
Wn=[50 6E3]./(60E3/2);%Cutoff frequencies divided by Nyquist frequency
[b,a]=butter(2,Wn);
load(fname,'peak_values')
fileN=ftype;
peak_values(peak_values(:,9)<20,:)=[];
peak_values=sortrows(peak_values,[21,6,7]);
chunk=unique(peak_values(:,6));
range_data=zeros(size(peak_values,1),108);
label=2.*ones(size(peak_values,1),1);
count=1;
for i=1:max(chunk)
    l1=find(peak_values(:,6)==i);
    load([fileN,'_',num2str(i),'_raw.mat'],'M')
    M(:,1)=(M(:,1)-mean(M(:,1)))./std(M(:,1));
    M(:,2)=(M(:,2)-mean(M(:,2)))./std(M(:,2));
    M(:,3)=(M(:,3)-mean(M(:,3)))./std(M(:,3));
    M_filt(:,4)=filtfilt(b,a,M(:,4));
    M(:,4)=(M_filt(:,4)-mean(M_filt(:,4)))./std(M_filt(:,4));
    disp(['          Chunk # ',num2str(i),' of ',num2str(max(chunk))])
    for j=1:size(l1,1)
        loc=peak_values(l1(j),7);
        M_range=[M(loc-13:loc+13,1)',M(loc-13:loc+13,2)',...
            M(loc-13:loc+13,3)',M(loc-13:loc+13,4)'];
        if peak_values(l1(j),4)<0.3 && peak_values(l1(j),5)>0.15
            label(count)=1;
        elseif peak_values(l1(j),4)>0.3
            label(count)=3;
        end
        range_data(count,:)=M_range;
        count=count+1;
    end
    clear M_filt
end
end
%% Mahal Formatting All
function [range_data]=PvFormat(peak_values)
range_data=peak_values(:,[1:3,8:16]);
range_data(:,1:3)=abs([range_data(:,1)./range_data(:,2),range_data(:,1)./range_data(:,3),range_data(:,2)./range_data(:,3)]);
end
%% Mahal Formatting Joe
function [range_data]=PvFormat2(peak_values)
range_data=peak_values(:,[1:3,8,9]);
range_data(:,1:3)=abs([range_data(:,1)./range_data(:,2),range_data(:,1)./range_data(:,3),range_data(:,2)./range_data(:,3)]);
end
%% Analysis Code
function [Out_mat,In_mat,Sig_mat]=OutlierAnalysis(outliers,everything)
% Find oultiers greater than max wildtype
outliers=abs(outliers);
everything=abs(everything);
% Outliers
r405_488=[mean((outliers(:,1)./outliers(:,2))),std((outliers(:,1)./outliers(:,2)))];
r405_633=[mean((outliers(:,1)./outliers(:,3))),std((outliers(:,1)./outliers(:,3)))];
r488_633=[mean((outliers(:,2)./outliers(:,3))),std((outliers(:,2)./outliers(:,3)))];
Cum=[mean(outliers(:,8)),std(outliers(:,8))];
FWHM=[mean(outliers(:,9)),std(outliers(:,9))];
GFLR=[mean(outliers(:,5)),std(outliers(:,5))];
AUC=[mean(outliers(:,10)),std(outliers(:,10))];
Out_mat = [r405_488;r405_633;r488_633;Cum;FWHM;GFLR;AUC];

%Inliers
r405_4882=[mean((everything(:,1)./everything(:,2))),std((everything(:,1)./everything(:,2)))];
r405_6332=[mean((everything(:,1)./everything(:,3))),std((everything(:,1)./everything(:,3)))];
r488_6332=[mean((everything(:,2)./everything(:,3))),std((everything(:,2)./everything(:,3)))];
Cum2=[mean(everything(:,8)),std(everything(:,8))];
FWHM2=[mean(everything(:,9)),std(everything(:,9))];
GFLR2=[mean(everything(:,5)),std(everything(:,5))];
AUC2=[mean(everything(:,10)),std(everything(:,10))];
In_mat = [r405_4882;r405_6332;r488_6332;Cum2;FWHM2;GFLR2;AUC2];

% Stats
r405488=outliers(:,1)./outliers(:,2);
r405633=outliers(:,1)./outliers(:,3);
r488633=outliers(:,2)./outliers(:,3);

r4054882=everything(:,1)./everything(:,2);
r4056332=everything(:,1)./everything(:,3);
r4886332=everything(:,2)./everything(:,3);

p1=ranksum(r405488,r4054882);
p2=ranksum(r405633,r4056332);
p3=ranksum(r488633,r4886332);
p4=ranksum(outliers(:,8),everything(:,8));
p5=ranksum(outliers(:,9),everything(:,9));
p6=ranksum(outliers(:,5),everything(:,5));
p7=ranksum(outliers(:,10),everything(:,10));

Sig_mat = [p1;p2;p3;p4;p5;p6;p7];
end