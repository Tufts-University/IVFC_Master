cd%% NV_033123_CVTesting
%% Initialization
clear
clc
%close all
mainpath = 'T:\Nilay\IVFC\Acquired Data\Blood Cell Data\DeepPeak2023';
cd(mainpath)
%% Main Loop
filecrit = 'NV_060823_FormattedDataSet_thresh*';
files = dir(filecrit);
for n = length(files):-1:1
    if contains(files(n).name,'locs') == 1
        files(n) = [];
    end
end
% files(1) = [];
CV_store = zeros(5,3);
thresh = [5,10,30,50,100];
exclusion = [9,11,17];
total_days = 34;
total_events = (total_days-length(exclusion)).*5;
events = zeros(total_events,2);

count = 1;
for i = 1:length(files)
    filename = files(i).name;
    load(filename)
    fields = fieldnames(FP_Peak_ranges);
    fields(exclusion) = [];
    total_events = zeros(size(fields,1),1);
    TP = zeros(size(fields,1),1);
    n = length(TP);
    for j = 1:length(total_events)
        data = FP_Peak_ranges.(fields{j});
        label = data(:,end);
        total_events(j,1) = length(label);
        TP(j,1) = length(label(label==0));
    end
    events(count:count+n-1,1) = TP./total_events;
    events(count:count+n-1,2) = repmat(thresh(i),total_days...
                                                -length(exclusion),1);
    CV_store(i,:) = [thresh(i),mean(TP./total_events),...
                               std(TP./total_events)];
    count = count + n;
end
CV = CV_store(:,3)./CV_store(:,2);
%% Gather Metrics From DeepPeak Thresholding
mainpath = 'T:\Nilay\IVFC\DeepPeakResults\Trial_2';
cd(mainpath)
files = dir('*All_Thresh*');
total_events = (total_days-length(exclusion)).*5;
xstore = zeros(total_events,1);
ystore = zeros(total_events,1);
count = 1;
savetable = zeros(5,3);
thresh = [100,10,30,50,5];
for i = 1:length(files)
    load(files(i).name)
    x(exclusion) = [];
    y(exclusion) = [];
    n = length(x)-1;
    xstore(count:count+n,1) = x';
    if length(unique(y))>1
       y = repmat(thresh(i),1,n+1);
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
rng(120)
Final_CV = zeros(5,2);
theoretical = [5, sqrt(5)./5; 10, sqrt(10)./10;30,sqrt(30)./30;50,sqrt(50)./50;100,sqrt(100)./100];
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
sortedCV = sortrows(Final_CV,1)
%% TestSet
mainpath = 'T:\Nilay\IVFC\Acquired Data\Blood Cell Data\DeepPeak2023';
cd(mainpath)
load('NV_060823_FormattedDataSet_thresh005.mat')
fields = fieldnames(FP_Peak_ranges);
keys = [1,2,6,8,11,14,18,21,22,28];
exclusion = [5];
fields_test = fields(keys);
filecrit = 'NV_060823_FormattedDataSet_thresh*';
files = dir(filecrit);
for n = length(files):-1:1
    if contains(files(n).name,'locs') == 1
        files(n) = [];
    end
end
%files(1) = [];
CV_store = zeros(5,3);
thresh = [5,10,30,50,100];
total_days = 10;
total_events = (total_days-length(exclusion)).*5;
events = zeros(total_events,2);

count = 1;
for i = 1:length(files)
    filename = files(i).name;
    load(filename)
    fields = fields_test;
    fields(exclusion) = [];
    total_events = zeros(length(fields),1);
    TP = zeros(length(fields),1);
    n = length(TP);
    for j = 1:length(total_events)
        data = FP_Peak_ranges.(fields{j});
        label = data(:,end);
        total_events(j,1) = length(label);
        TP(j,1) = length(label(label==0));
    end
    events(count:count+n-1,1) = TP./total_events;
    events(count:count+n-1,2) = repmat(thresh(i),length(fields),1);
    CV_store(i,:) = [thresh(i),mean(TP./total_events),...
                               std(TP./total_events)];
    count = count + n;
end
CV = CV_store(:,3)./CV_store(:,2);
%% Test Set Sum of Variance
average = zeros(5,2);
seeds = randi(500,20,1);
for n = 1:length(seeds)
    seed = seeds(n);
    rng(seed)
    Final_CV2 = zeros(5,2);
    for k = 1:length(thresh)
        x_theoretical = poissrnd(thresh(k),[9,1]);
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
    average = average + sortedCV2;
end
average = average./20;
