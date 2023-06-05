%% NV_052323_CorrelationPlots_DeepPeak_Wavelength
%% Initialization
clear
clc
mainpath = 'T:\Nilay\IVFC\DeepPeakResults\Trial_2\k-fold\Seed 2';
cd(mainpath)
folders = dir('*combo');
[F_names] = natsortfiles({folders.name});
%% Plotting code
PearsonStorage = zeros(length(folders),7);
for i = 1:length(F_names)
    disp(['Plotting Combo # ', num2str(i),' of ',...
                        num2str(length(F_names))])
    cd(F_names{i})
    load(['Pearson_All_seed2'])
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
    load(['Pearson_seed2'])
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
    disp(['Calculating Combo # ', num2str(i),' of ',...
                        num2str(length(F_names))])
    cd(F_names{i})
    load(['Pearson_FP_seed2'])
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
%% Load Combo Sheet
T = readtable("Combo_Sheet.xlsx");
values = table2array(T(:,2:end));

%% Split into Test and All
Test_values = values(:,1:5)./100;
All_values = values(:,6:10)./100;
Test_values(:,5) = Test_values(:,5).*100;
All_values(:,5) = All_values(:,5).*100;
All_values = All_values(:,[3,4,2,1,5]);
Test_values = Test_values(:,[3,4,2,1,5]);
FAR_values = values(:,11);
%% Plot Test performance
colors = linspecer(5);
Metrics = {'Purity','Specificity','Sensitivity','Accuracy',...
            'Pearson Corr.'};
figure('units','inch','position',[2,2,7.5,6]);
b = bar((1:5),[Test_values],'grouped','LineWidth',1.5,...
        'FaceAlpha',0.9);
colororder(colors)
boldify
set(gcf,'color','w');
set(gca,'FontSize',12,'FontName','Ariel','LineWidth',1.5);
yticks([0:0.1:1])
combos = table2cell(T(:,1));
for k = 1:length(combos)
    temp = num2str(combos{k});
    if length(temp)>6
        temp = insertAfter(temp,6,'+');
        temp = insertAfter(temp,3,'+');
    elseif length(temp)>3
        temp = insertAfter(temp,3,'+');
    else
        temp = temp;
    end
    combos{k} = temp;
end
xticklabels(combos)
legend(Metrics,'location','southoutside','numcolumns',5,...
        'FontSize',12,'FontName','Ariel')
ylabel('Performance')
print(gcf, '-dmeta', 'Combos.emf');
%% Plot Combos All performance
colors = linspecer(5);
Metrics = {'Purity','Specificity','Sensitivity','Accuracy',...
            'Pearson Corr.'};
figure('units','inch','position',[2,2,7.5,6]);
b = bar((1:5),[All_values],'grouped','LineWidth',1.5,...
        'FaceAlpha',0.9);
colororder(colors)
boldify
set(gcf,'color','w');
set(gca,'FontSize',12,'FontName','Ariel','LineWidth',1.5);
yticks([0:0.1:1])
combos = table2cell(T(:,1));
for k = 1:length(combos)
    temp = num2str(combos{k});
    if length(temp)>6
        temp = insertAfter(temp,6,'+');
        temp = insertAfter(temp,3,'+');
    elseif length(temp)>3
        temp = insertAfter(temp,3,'+');
    else
        temp = temp;
    end
    combos{k} = temp;
end
xticklabels(combos)
legend(Metrics,'location','southoutside','numcolumns',5,...
        'FontSize',12,'FontName','Ariel')
ylabel('Performance')
print(gcf, '-dmeta', 'Combos_All.emf');
%% Plot Test performance v2
colors = linspecer(5);
Metrics = {'Purity','Specificity','Sensitivity','Accuracy',...
            'Pearson Corr.'};
figure('units','inch','position',[2,2,7.5,6]);
b = bar((1:5),[Test_values'],'grouped','LineWidth',1.5,...
        'FaceAlpha',0.9);
colororder(colors)
boldify
set(gcf,'color','w');
set(gca,'FontSize',12,'FontName','Ariel','LineWidth',1.5);
yticks([0:0.1:1])
combos = table2cell(T(:,1));
for k = 1:length(combos)
    temp = num2str(combos{k});
    if length(temp)>6
        temp = insertAfter(temp,6,'+');
        temp = insertAfter(temp,3,'+');
    elseif length(temp)>3
        temp = insertAfter(temp,3,'+');
    else
        temp = temp;
    end
    combos{k} = temp;
end
xticklabels(Metrics)
legend(combos,'location','southoutside','numcolumns',5,...
        'FontSize',12,'FontName','Ariel')
ylabel('Performance')
print(gcf, '-dmeta', 'Combos_v2.emf');
%% Plot Combos All performance
colors = linspecer(5);
Metrics = {'Purity','Specificity','Sensitivity','Accuracy',...
            'Pearson Corr.'};
figure('units','inch','position',[2,2,7.5,6]);
b = bar((1:5),[All_values'],'grouped','LineWidth',1.5,...
        'FaceAlpha',0.9);
colororder(colors)
boldify
set(gcf,'color','w');
set(gca,'FontSize',12,'FontName','Ariel','LineWidth',1.5);
yticks([0:0.1:1])
combos = table2cell(T(:,1));
for k = 1:length(combos)
    temp = num2str(combos{k});
    if length(temp)>6
        temp = insertAfter(temp,6,'+');
        temp = insertAfter(temp,3,'+');
    elseif length(temp)>3
        temp = insertAfter(temp,3,'+');
    else
        temp = temp;
    end
    combos{k} = temp;
end
xticklabels(Metrics)
legend(combos,'location','southoutside','numcolumns',5,...
        'FontSize',12,'FontName','Ariel')
ylabel('Performance')
print(gcf, '-dmeta', 'Combos_All_v2.emf');