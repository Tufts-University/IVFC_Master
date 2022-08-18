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
l = 1;
LeanCancerMice(4)=struct();
for i=1:length(Lean_dates)%Loop through Lean Mice Days
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
l = 1;
ObeseCancerMice(4)=struct();
for i=1:length(Obese_dates)%Loop through Lean Mice Days
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