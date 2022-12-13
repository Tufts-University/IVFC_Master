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
Thresholds = Lean_GFLR_thresh2;
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
        peaks = [peaks; peak_values(peak_values(:,5)>Thresholds,:)];
        counts(i,j) = size(peak_values(peak_values(:,5)>Thresholds,:),1)';
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
        peaks2 = [peaks2; peak_values(peak_values(:,5)>Thresholds,:)];
        counts2(i,j) = size(peak_values(peak_values(:,5)>Thresholds,:),1)';
    end
    peaks_total2 = [peaks_total2;peaks2];
end
%% Unnormalized counts
colors = linspecer(2);
lean_norm = counts;
obese_norm = counts2;
figure; plot(1:1:4,mean(lean_norm,2,'omitnan'),'.','Color',colors(1,:),'MarkerSize',30)
hold on
plot(1:1:4,mean(obese_norm,2,'omitnan'),'.','Color',colors(2,:),'MarkerSize',30)
title('Average Number of Events Per Timepoint')
xticks(1:1:4)
xticklabels({'8 Weeks', '12 Weeks', '16 Weeks', '20 Weeks'})
xlabel('Timepoint')
ylabel('Events')
boldify
errorbar(1:1:4,mean(lean_norm,2,'omitnan'),std(lean_norm,[],2,'omitnan'),'LineWidth',3,'CapSize',10,'Color',colors(1,:))
errorbar(1:1:4,mean(obese_norm,2,'omitnan'),std(obese_norm,[],2,'omitnan'),'LineWidth',2,'CapSize',10,'Color',colors(2,:))
legend('Lean Cancer Mice','Obese Cancer Mice','numcolumns',2,'Location','Southoutside')
xlim([0.5,4.5])
%% Unnormalized count Box Plot
timepts = repmat({'8 Weeks', '12 Weeks', '16 Weeks', '20 Weeks'},3,1);
timepts = [timepts(:);timepts(:)];
Diet = [repmat({'Lean'},12,1);repmat({'Obese'},12,1)];
counts(2,3)=nan;
Lean1 = counts';
Lean = Lean1(:);
Obese1 = counts2';
Obese = Obese1(:);
Events = [Lean;Obese];
T = table(timepts,Diet,Events);
dataOrder = {'8 Weeks', '12 Weeks', '16 Weeks', '20 Weeks'};
T.timepts = categorical(T.timepts,dataOrder);
colors = linspecer(2);
colororder(colors)
boxchart(T.timepts,T.Events,'GroupByColor',T.Diet)
title('Average Number of Events Per Timepoint')
xlabel('Timepoint')
xtickangle(45)
ylabel('Events')
boldify
ylim([-.1,inf])
legend('Lean Cancer Mice','Obese Cancer Mice','numcolumns',2,'Location','Southoutside')
%% Normalized count Box Plot
figure;
Lean_time =[60.1066666666667	55.1016666666667	60.1050000000000
28.6875766666667	60.1050000000000	0
60.1016666666667	60.6683333333333	61.9542433333333
60.1000000000000	60.1016666666667	60.1066666666667];
Obese_time =[45.0783333333333	60.1066666666667	60.2438888888889
80.1350000000000	60.1066666666667	60.1016666666667
60.1050000000000	60.1066666666667	60.1016666666667
46.2492433333333	60.1066666666667	60.1033333333333];

timepts = repmat({'8 Weeks', '12 Weeks', '16 Weeks', '20 Weeks'},3,1);
timepts = [timepts(:);timepts(:)];
Diet = [repmat({'Lean'},12,1);repmat({'Obese'},12,1)];
counts(2,3)=nan;
Lean1 = counts./Lean_time;
Lean1 = Lean1';
Lean = Lean1(:);
Obese1 = counts2./Obese_time;
Obese1 = Obese1';
Obese = Obese1(:);
Events = [Lean;Obese];
T = table(timepts,Diet,Events);
dataOrder = {'8 Weeks', '12 Weeks', '16 Weeks', '20 Weeks'};
T.timepts = categorical(T.timepts,dataOrder);
colors = linspecer(2);
colororder(colors)
boxchart(T.timepts,T.Events,'GroupByColor',T.Diet)
title('Average Number of Events Per Timepoint')
xlabel('Timepoint')
xtickangle(45)
ylabel('Events/Minute (min^{-1})')
boldify
legend('Lean Cancer Mice','Obese Cancer Mice','numcolumns',2,'Location','Southoutside')
%% Normalized counts
colors = linspecer(2);
lean_norm = counts./Lean_time;
obese_norm = counts2./Obese_time;
figure; plot(1:1:4,mean(lean_norm,2,'omitnan'),'.','Color',colors(1,:),'MarkerSize',30)
hold on
plot(1:1:4,mean(obese_norm,2,'omitnan'),'.','Color',colors(2,:),'MarkerSize',30)
title('Average Number of Events Per Timepoint')
xticks(1:1:4)
xticklabels({'8 Weeks', '12 Weeks', '16 Weeks', '20 Weeks'})
xlabel('Timepoint')
ylabel('Events/Minute (min^{-1})')
boldify
errorbar(1:1:4,mean(lean_norm,2,'omitnan'),std(lean_norm,[],2,'omitnan'),'LineWidth',3,'CapSize',10,'Color',colors(1,:))
errorbar(1:1:4,mean(obese_norm,2,'omitnan'),std(obese_norm,[],2,'omitnan'),'LineWidth',2,'CapSize',10,'Color',colors(2,:))
legend('Lean Cancer Mice','Obese Cancer Mice','numcolumns',2,'Location','Southoutside')
xlim([0.5,4.5])
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