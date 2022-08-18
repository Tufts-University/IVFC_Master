%% NV_072222_ObesityStudy %%
%% Initialization
clear 
clc
close all
mainpath = 'T:\Nilay\IVFC\Acquired Data\SydneyMice';
cd(mainpath)
%% Load in all Lean Healthy Samples
dirinfo=dir("*HealthyControl");
p_total=[];
for i=1:length(dirinfo)
    cd([dirinfo(i).folder,'\',dirinfo(i).name])
    cd(dir("*Lean").name)
    load(dir("*NoScatAll.mat").name)
    p_total=[p_total;peak_values];
end
%% Final characteristics of all peaks
pv=p_total(p_total(:,19)>19,:);
Lean_GFLR_thresh = max(pv(:,5));
Lean_GFLR_thresh2 = max(p_total(:,5));
%% Load in all Obese Healthy Samples
cd(mainpath)
dirinfo=dir("*HealthyControl");
p_total2=[];
for i=1:length(dirinfo)
    cd([dirinfo(i).folder,'\',dirinfo(i).name])
    cd(dir("*Obese").name)
    load(dir("*NoScatAll.mat").name)
    p_total2=[p_total2;peak_values];
end
%% Final characteristics of all peaks
pv2=p_total2(p_total2(:,19)>19,:);
Obese_GFLR_thresh = max(pv2(:,5));
Obese_GFLR_thresh2 = max(p_total2(:,5));
%% ThresholdingSet
Thresholds = Lean_GFLR_thresh;
cd(mainpath)
Obese_dates = {'041322';'051122';'060822';'070722'};
Lean_dates = {'041422';'051222';'060822';'070622'};
%% Lean Data
peaks_total = [];
counts = zeros([4,3]);
for i=1:length(Lean_dates)%Loop through Lean Mice Days
    cd(mainpath)
    disp(['Currently Processing Day # ',num2str(i),' of ',num2str(length(Lean_dates))])
    foldername = ['NV_',Lean_dates{i},'_LeanCancer'];
    LeanMice = dir([foldername,'*']);
    peaks = [];
    for j=1:length(LeanMice) %Loop through sample per day
        disp(['     Currently Processing Sample # ',num2str(j),' of ',num2str(length(LeanMice))])
        cd([LeanMice(j).folder,'\',LeanMice(j).name])
        load(dir("*NoScatAll.mat").name)
        peaks = [peaks; peak_values];
        counts(i,j) = size(peak_values,1)';
    end
    peaks_total = [peaks_total;peaks];
end
%% Obese Data
peaks_total2 = [];
counts2 = zeros([4,3]);
for i=1:length(Obese_dates)%Loop through Lean Mice Days
    cd(mainpath)
    disp(['Currently Processing Day # ',num2str(i),' of ',num2str(length(Obese_dates))])
    foldername = ['NV_',Obese_dates{i},'_ObeseCancer'];
    ObeseMice = dir([foldername,'*']);
    peaks2 = [];
    for j=1:length(ObeseMice) %Loop through sample per day
        disp(['     Currently Processing Sample # ',num2str(j),' of ',num2str(length(ObeseMice))])
        cd([ObeseMice(j).folder,'\',ObeseMice(j).name])
        load(dir("*NoScatAll.mat").name)
        peaks2 = [peaks2; peak_values];
        counts2(i,j) = size(peak_values,1)';
    end
    peaks_total2 = [peaks_total2;peaks2];
end

%% MAHAL Codes
cd(mainpath)
PV_Totals_lean_control = PvFormat(p_total);
dirinfo=dir("*HealthyControl");
mahalLean_store=[];
for i=1:length(dirinfo)
    cd([dirinfo(i).folder,'\',dirinfo(i).name])
    cd(dir("*Lean").name)
    load(dir("*NoScatAll.mat").name)
    mahalLeanControl = mahal(PvFormat(peak_values),PV_Totals_lean_control);
    mahalLean_store=[mahalLean_store;mahalLeanControl];
end
cd(mainpath)
PV_Totals_obese_control = PvFormat(p_total2);
dirinfo=dir("*HealthyControl");
mahalObese_store=[];
for i=1:length(dirinfo)
    cd([dirinfo(i).folder,'\',dirinfo(i).name])
    cd(dir("*Obese").name)
    load(dir("*NoScatAll.mat").name)
    mahalobeseControl = mahal(PvFormat(peak_values),PV_Totals_obese_control);
    mahalObese_store=[mahalObese_store;mahalobeseControl];
end
%% Running Codes for PV formating
function [range_data]=PvFormat(peak_values);
range_data=peak_values(:,[1:3,8:16]);
range_data(:,1:3)=abs([range_data(:,1)./range_data(:,2),range_data(:,1)./range_data(:,3),range_data(:,2)./range_data(:,3)]);
end

function [range_data]=PvFormat2(peak_values);
range_data=peak_values(:,[1:3,8,9]);
range_data(:,1:3)=abs([range_data(:,1)./range_data(:,2),range_data(:,1)./range_data(:,3),range_data(:,2)./range_data(:,3)]);
end