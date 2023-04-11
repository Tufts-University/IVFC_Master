%% NV_033123_CVTesting
%% Initialization
clear
clc
close all
mainpath = 'T:\Nilay\IVFC\Acquired Data\Blood Cell Data\DeepPeak2023';
cd(mainpath)
%% Main Loop
filecrit = 'NV_032823_FormattedDataSet_thresh*';
files = dir(filecrit);
files(1) = [];
CV_store = zeros(4,3);
thresh = [10,30,50,100];
events = zeros(136,2);
count = 1;
for i = 1:length(files)
    filename = files(i).name;
    load(filename)
    fields = fieldnames(FP_Peak_ranges);
    total_events = zeros(34,1);
    TP = zeros(34,1);
    n = length(TP);
    for j = 1:length(total_events)
        data = FP_Peak_ranges.(fields{j});
        label = data(:,end);
        total_events(j,1) = length(label);
        TP(j,1) = length(label(label==0));
    end
    events(count:count+n-1,1) = TP./total_events;
    events(count:count+n-1,2) = repmat(thresh(i),34,1);
    CV_store(i,:) = [thresh(i),mean(TP./total_events),...
                               std(TP./total_events)];
    count = count + n;
end
CV = CV_store(:,3)./CV_store(:,2)
%% Gather Metrics From DeepPeak Thresholding
mainpath = 'T:\Nilay\IVFC\DeepPeakResults\Trial_2';
cd(mainpath)
files = dir('*All_Thresh*');
xstore = zeros(136,1);
ystore = zeros(136,1);
count = 1;
savetable = zeros(4,3);
thresh = [10,100,30,50];
for i = 1:length(files)
    load(files(i).name)
    n = length(x)-1;
    xstore(count:count+n,1) = x';
    if length(unique(y))>1
       y = repmat(mode(y),1,n+1);
    end
    ystore(count:count+n,1) = y';
    count = count + n+1;
    savetable(i,:) = [thresh(i),mean(double(x)),std(double(x))];
end
savetable = sortrows(savetable,1);
values = [xstore,ystore];
sortedvalues = sortrows(values,2);
x = sortedvalues(:,1);
y = sortedvalues(:,2);
%% Sum of variance
rng(60)
Final_CV = zeros(4,2);
theoretical = [10, 0.3;30,0.19;50,0.15;100,0.11];
for k = 1:length(thresh)
    x_theoretical = poissrnd(thresh(k),[length(find(sortedvalues(:,2)...
                                                    ==thresh(k))),1]);
    x_theoretical = (x_theoretical-min(x_theoretical))./(max(x_theoretical)-min(x_theoretical));
    x_purity = events(events(:,2)==thresh(k),1);
    x_purity = (x_purity-min(x_purity))./(max(x_purity)-min(x_purity));
    covariance_mat = cov(x_theoretical,x_purity);
    stdev_sample = sqrt(covariance_mat(1,1)+covariance_mat(2,2) + ...
                    2.*covariance_mat(1,2));
    mean_sample = mean(x_theoretical) + mean(x_purity);
    Final_CV(k,:) = [thresh(k),stdev_sample./mean_sample];
end
sortedCV = sortrows(Final_CV,1);

%% TestSet
mainpath = 'T:\Nilay\IVFC\Acquired Data\Blood Cell Data\DeepPeak2023';
cd(mainpath)
load('NV_032823_FormattedDataSet_thresh001.mat')
fields = fieldnames(FP_Peak_ranges);
keys = [1,2,6,8,11,14,18,21,22,28];
fields_test = fields(keys);
filecrit = 'NV_032823_FormattedDataSet_thresh*';
files = dir(filecrit);
files(1) = [];
CV_store = zeros(4,3);
thresh = [10,30,50,100];
events = zeros(40,2);
count = 1;
for i = 1:length(files)
    filename = files(i).name;
    load(filename)
    fields = fields_test;
    total_events = zeros(10,1);
    TP = zeros(10,1);
    n = length(TP);
    for j = 1:length(total_events)
        data = FP_Peak_ranges.(fields{j});
        label = data(:,end);
        total_events(j,1) = length(label);
        TP(j,1) = length(label(label==0));
    end
    events(count:count+n-1,1) = TP./total_events;
    events(count:count+n-1,2) = repmat(thresh(i),10,1);
    CV_store(i,:) = [thresh(i),mean(TP./total_events),...
                               std(TP./total_events)];
    count = count + n;
end
CV = CV_store(:,3)./CV_store(:,2)
%% Test Set Sum of Variance
rng(60)
Final_CV2 = zeros(4,2);
theoretical = [10, 0.32;30,0.20;50,0.16;100,0.13];
for k = 1:length(thresh)
    x_theoretical = poissrnd(thresh(k),[10,1]);
    x_theoretical = (x_theoretical-min(x_theoretical))./(max(x_theoretical)-min(x_theoretical));
    x_purity = events(events(:,2)==thresh(k),1);
    x_purity = (x_purity-min(x_purity))./(max(x_purity)-min(x_purity));
    covariance_mat = cov(x_theoretical,x_purity);
    stdev_sample = sqrt(covariance_mat(1,1)+covariance_mat(2,2) + ...
                    2.*covariance_mat(1,2));
    mean_sample = mean(x_theoretical) + mean(x_purity);
    Final_CV2(k,:) = [thresh(k),stdev_sample./mean_sample];
end
sortedCV2 = sortrows(Final_CV2,1);
