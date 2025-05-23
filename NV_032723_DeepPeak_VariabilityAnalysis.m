%% NV_032723_DeepPeak_VariabilityAnalysis
%% Initialization
clear
clc
close all
mainpath = 'T:\Nilay\IVFC\DeepPeakResults\Trial_2';
cd(mainpath)
datapath = 'T:\Nilay\IVFC\Acquired Data\Blood Cell Data\DeepPeak2023';
cd(datapath)
exlusiondata = readtable('NV_060823_ExclusionTable_Threshold.xlsx','NumHeaderLines',1);
cd(mainpath);
%% Main Code
colors = colororder;
%% All Data Thresh
files = dir('*All_Thresh*');
alldata = table2array(exlusiondata(:,2:6));
total = sum(alldata,'all');
xstore = zeros(total,1);
ystore = zeros(total,1);
count = 1;
savetable = zeros(4,3);
alldata = alldata(:,[5,2,3,4,1]);
thresh = [100,10,30,50,5];
for i = 1:length(files)
    load(files(i).name)
    n = length(x)-1;
    x = x.*alldata(i);
    y = y.*alldata(i);
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
mdl = fitlm(savetable(:,1),savetable(:,2))
plot(mdl,'Marker','.','MarkerSize',23,'Color',colors(1,:))
xlabel('Expected Number of Events')
ylabel('Average DeepPeak Events')
title('Expected vs. DeepPeak Average Counts (Full Dataset)')
boldify
yticks(0:10:100)
str_r2 = num2str(mdl.Rsquared.Adjusted);
str_r2 = str_r2(3:end);
file_name = ['Seed2_SpikingRatio_FullDataset_LinFit_r_0p',str_r2];
 print(file_name,'-dsvg')

figure; 
boxchart(categorical(y),x,"MarkerStyle",".","MarkerSize",23)
title('Recovery of Spiked CTCCs (Full Dataset)')
xlabel('Expected Number of Events')
ylabel('Observed DeepPeak Events')
boldify
file_name = 'Seed2_SpikingRatio_FullDataset_BoxChart';
print(file_name,'-dsvg')
%% Test Data Thresh
files = dir('Pearson_Thresh*');
xstore = zeros(45,1);
ystore = zeros(45,1);
count = 1;
savetable = zeros(4,3);
thresh = [10,100,30,5,50];
for i = 1:length(files)
    load(files(i).name)
%     if i~=1 && i~=5
        x(5) = [];
        y(5) = [];
%     end
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
mdl = fitlm(savetable(:,1),savetable(:,2))
figure;
plot(mdl,'Marker','.','MarkerSize',23,'Color',colors(1,:))
xlabel('Expected Number of Events')
ylabel('Average DeepPeak Events')
title('Expected vs. DeepPeak Average Counts (Test Set)')
yticks(0:10:100)
boldify
str_r2 = num2str(mdl.Rsquared.Adjusted);
str_r2 = str_r2(3:end);
file_name = ['Seed2_SpikingRatio_TestDataset_LinFit_r_0p',str_r2];
print(file_name,'-dsvg')


figure; 
boxchart(categorical(y),x,"MarkerStyle",".","MarkerSize",23)
title('Recovery of Spiked CTCCs (Test Set)')
xlabel('Expected Number of Events')
ylabel('Observed DeepPeak Events')
boldify
file_name = 'Seed2_SpikingRatio_TestDataset_BoxChart';
print(file_name,'-dsvg')