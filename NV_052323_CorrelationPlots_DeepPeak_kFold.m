%% NV_052323_CorrelationPlots_DeepPeak_kfold
%% Initialization
clear
clc
mainpath = 'T:\Nilay\IVFC\DeepPeakResults\Trial_2\k-fold';
cd(mainpath)
folders = dir('Seed*');
[F_names] = natsortfiles({folders.name});
%% Plotting code
PearsonStorage = zeros(length(folders),7);
for i = 1:length(F_names)
    disp(['Plotting Seed # ', num2str(i),' of ',...
                        num2str(length(F_names))])
    cd(F_names{i})
    load(['Pearson_All_seed',num2str(i-1)])
    x = double(x);
    y = double(y);
    figure; plot(x,y,'b.','MarkerSize',23)
    xlabel('DeepPeak Detected Events')
    ylabel('True # of CTCCs')
    a = fitlm(x,y);
    hold on
    plot(a)
    r = corrcoef(x,y);
    disp([9 'Pearson Correlation = ', num2str(r(1,2)),...
                       ', R^2 = ',num2str(r(1,2).^2)])
    xlabel('DeepPeak Detected Events')
    ylabel('True # of CTCCs')
    title('DeepPeak vs. True Events')
    hold on
    upperbound = max(max(x),5000);
    plot((0:100:upperbound),(0:100:upperbound),'k--')
    ylim([-10,5000])
    xlim([-10,upperbound])
    boldify
    legend('Data','','Fit','Confidence Bound','','y=x')
    set(gcf,'color','w');
    FileName = 'LinFit_FullData';
    print(FileName,'-dpng')
    PearsonStorage(i,1:3) = [str2double(F_names{i}(find(F_names{i}==' ')+1:end)),r(1,2),r(1,2).^2];
    %% Lower Range Fit
    locs = find(x<max((max(x).*0.3),1500));
    if ~isempty(locs)
        x2 = x(locs);
        y2 = y(locs);
        
        figure; plot(x2,y2,'b.','MarkerSize',23)
        xlabel('DeepPeak Detected Events')
        ylabel('True # of CTCCs')
        a = fitlm(x2,y2);
        hold on
        plot(a)
        r = corrcoef(x2,y2);
        disp([9 'Pearson Correlation = ', num2str(r(1,2)),...
                            ', R^2 = ',num2str(r(1,2).^2)])
        xlabel('DeepPeak Detected Events')
        ylabel('True # of CTCCs')
        title('DeepPeak vs. True Events')
        hold on
        upperbound = max((max(x).*0.3),1500);
        plot((0:100:upperbound),(0:100:upperbound),'k--')
        ylim([-10,upperbound])
        xlim([-10,upperbound])
        boldify
        legend('Data','','Fit','Confidence Bound','','y=x')
        set(gcf,'color','w');
        FileName = 'LinFit_FullData_LowerRange';
        print(FileName,'-dpng')
        PearsonStorage(i,4:5) = [r(1,2),r(1,2).^2];
    end
    %% Test Set Fit
    load(['Pearson_seed',num2str(i-1)])
    x = double(x);
    y = double(y);
    figure; plot(x,y,'b.','MarkerSize',23)
    xlabel('DeepPeak Detected Events')
    ylabel('True # of CTCCs')
    a = fitlm(x,y);
    hold on
    plot(a)
    r = corrcoef(x,y);
    disp([9 'Pearson Correlation = ', num2str(r(1,2)),...
                        ', R^2 = ',num2str(r(1,2).^2)])
    xlabel('DeepPeak Detected Events')
    ylabel('True # of CTCCs')
    title('DeepPeak vs. True Events')
    hold on
    upperbound = max(max(x),5000);
    plot((0:100:upperbound),(0:100:upperbound),'k--')
    ylim([-10,5000])
    xlim([-10,upperbound])
    boldify
    legend('Data','','Fit','Confidence Bound','','y=x')
    set(gcf,'color','w');
    FileName = 'LinFit_TestSet';
    pause(1)
    print(FileName,'-dpng')
    PearsonStorage(i,6:7) = [r(1,2),r(1,2).^2];
    close all
    cd(mainpath)
end
%% FAR Calculation 
cd(mainpath)
FAR = zeros(7,1);
for i = 1:length(F_names)
    disp(['Calculating Seed # ', num2str(i),' of ',...
                        num2str(length(F_names))])
    cd(F_names{i})
    load(['Pearson_FP_seed',num2str(i-1)])
    x = double(x);
    y = double(y);
    exclusion = [4,6];
    x_final = x;
    x_final(exclusion) = [];
    total_FP = sum(x_final);
    total_time = 175; % mins
    flow_rate = 3; % uL/min
    total_volume = total_time.*flow_rate;
    FAR(i) = total_FP./total_volume;
    cd(mainpath)
end
%% Load k-fold Sheet
T = readtable("k-fold_Sheet.xlsx");
values = table2array(T(:,2:end));
% Mean Performance
avg_fold = mean(values);
% Variance
std_fold = std(values);
%% Split into Test and All
Test_avg_fold = avg_fold(:,1:5)./100;
All_avg_fold = avg_fold(:,6:10)./100;
Test_std_fold = std_fold(:,1:5)./100;
All_std_fold = std_fold(:,6:10)./100;

Test_avg_fold(:,5) = Test_avg_fold(:,5).*100;
All_avg_fold(:,5) = All_avg_fold(:,5).*100;
Test_std_fold(:,5) = Test_std_fold(:,5).*100;
All_std_fold(:,5) = All_std_fold(:,5).*100;

avg_FAR = avg_fold(:,11);
std_FAR = std_fold(:,11);

%% Plot k-fold performance
colors = linspecer(5);
Metrics = {'Accuracy','Sensitivity','Purity','Specificity',...
            'Pearson Corr.'};
figure('units','inch','position',[2,2,7.5,6]);
b = bar((1:5),[Test_avg_fold;All_avg_fold],'grouped','LineWidth',1.5,...
        'FaceAlpha',0.9);
colororder(colors)
boldify
set(gcf,'color','w');
set(gca,'FontSize',12,'FontName','Ariel','LineWidth',1.5);
yticks([0:0.1:1])
hold on

totalPerform_mean = [Test_avg_fold;All_avg_fold]';
totalPerform_err = [Test_std_fold;All_std_fold]';

% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(totalPerform_mean);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',totalPerform_mean,totalPerform_err,'k','linestyle','none',...
    'LineWidth',3,'CapSize',9);
xticklabels(Metrics)
legend('Test Set','All Data','location','southoutside','numcolumns',2,...
        'FontSize',12,'FontName','Ariel')
ylabel('Performance')
print(gcf, '-dmeta', 'k-fold_validation.emf');
%% Plot k-fold performance v2
colors = linspecer(5);
Metrics = {'Accuracy','Sensitivity','Purity','Specificity',...
            'Pearson Corr.'};
figure('units','inch','position',[2,2,7.5,6]);
data = [Test_avg_fold;All_avg_fold]';
data = data([3,4,2,1,5],:);
b = bar((1:2),data','grouped','LineWidth',1.5,...
        'FaceAlpha',0.9);
colororder(colors)
boldify
set(gcf,'color','w');
set(gca,'FontSize',12,'FontName','Ariel','LineWidth',1.5);
yticks([0:0.1:1])
hold on

totalPerform_mean = data';
totalPerform_err = [Test_std_fold;All_std_fold];

totalPerform_err = totalPerform_err(:,[3,4,2,1,5]);
% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(totalPerform_mean);
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
errorbar(x',totalPerform_mean,totalPerform_err,'k','linestyle','none',...
    'LineWidth',3,'CapSize',9);
xticklabels({'Test Data','All Data'})
legend('Purity','Specificity','Sensitivity','Accuracy','Pearson Corr.',...
        'Location','southoutside','numcolumns',5,...
        'FontSize',12,'FontName','Ariel')
ylabel('Performance')
print(gcf, '-dmeta', 'k-fold_validation_V2.emf');