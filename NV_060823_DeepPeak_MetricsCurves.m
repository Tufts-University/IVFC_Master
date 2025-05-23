%% NV_060823_MetricsCurves
%% Initialization
% clear
% clc
% close all
mainpath = 'T:\Nilay\IVFC\DeepPeakResults\Trial_2\k-fold\Seed 2\Metrics';
cd(mainpath)
folders = dir('*Cell_Metrics.mat');
[F_names] = natsortfiles({folders.name});
%% Load and restore values
Sensitivity_Final = zeros(10,10);
Purity_Final = zeros(10,10);
Specificity_Final = zeros(10,10);
for i = 1:10
    load(F_names{i});
    Sensitivity_Final(i,:) = Sensitivity;
    Purity_Final(i,:) = Purity;
    Specificity_Final(i,:) = Specificity;
end
%% Plotting
colors = linspecer(10);
figure;
plot((1:10),Sensitivity_Final')
boldify
ylabel('Sensitivity')
xlabel('Model #')
caption = cell(10,1);
for k = 1:10
    caption{k} = ['Day # ',num2str(k)];
end
colororder(colors);
xticks((0:1:10))
legend(caption,'NumColumns',5,'Location','southoutside')

figure;
plot((1:10),Purity_Final')
boldify
ylabel('Purity')
xlabel('Model #')
colororder(colors);
xticks((0:1:10))
legend(caption,'NumColumns',5,'Location','southoutside')

%% Pearson Correlation
folders = dir('*Pearson_Metrics.mat');
[F_names] = natsortfiles({folders.name});
load(F_names{1});
x = double(squeeze(Pearson(:,:,1)));
y = double(squeeze(Pearson(:,:,2)));
% keys = [1,2,6,8,11,14,18,21,22,28];
% x = x(keys,:);
% y = y(keys,:); 
Pearson_store = zeros(10,1);
for j = 1:10
    disp(['Model # ',num2str(j),' of 10'])
    x_temp = x(:,j);
    y_temp = y(:,j);
    r = corrcoef(x_temp,y_temp);
    disp([9 'Pearson Correlation = ', num2str(r(1,2)),...
                       ', R^2 = ',num2str(r(1,2).^2)])
    Pearson_store(j,1) = r(1,2);
end
figure;
plot((1:10),Pearson_store);
boldify
ylabel('Pearson Correlation')
xlabel('Model #')
xticks((0:1:10))
%% Weighted Sensitivity
totalevents = sum(y(:,1));
ratioevents = y(:,1)./totalevents.*100;
weigthedSensitivity = Sensitivity_Final .* ratioevents;
weigthedPurity = Purity_Final .* ratioevents;
meanSens = sum(weigthedSensitivity);
meanPur = sum(Purity_Final).*10;
figure;
plot((1:10),meanSens.*0.8509);
boldify
ylabel('Net Sensitivity (%)')
xlabel('Model #')
xticks((0:1:10))

figure;
plot((1:10),meanPur);
boldify
ylabel('Purity (%)')
xlabel('Model #')
xticks((0:1:10))
%% FAR
exclusion = [4,6];
folders = dir('*FAR_Metrics.mat');
[F_names] = natsortfiles({folders.name});
load(F_names{1});
FAR_events = double(FAR);
FAR_events(exclusion,:) = [];
FAR_store = zeros(10,1);
for j = 1:10
    disp(['Model # ',num2str(j),' of 10'])
    events = FAR_events(:,j);
    total_FP = sum(events);
    total_time = 175; % mins
    flow_rate = 3; % uL/min
    total_volume = total_time.*flow_rate;
    FAR_store(j) = total_FP./total_time;
end
figure;
plot((1:10),FAR_store);
boldify
ylabel('False Alarm Rate (\muL^{-1})')
xlabel('Model #')
xticks((0:1:10))
%% Joint Plot
colors = linspecer(4);
figure;
yyaxis left
plot((1:10),Pearson_store.*100,'-','Color',colors(1,:));
ylabel('Performance (%)')
xlabel('Model #')
xticks((0:1:10))
hold on
plot((1:10),meanSens.*0.8509,'-','Color',colors(2,:));
plot((1:10),meanPur,'-','Color',colors(3,:));
boldify
yyaxis("right")
plot((1:10),FAR_store,'Color',colors(4,:));
ylabel('False Alarm Rate (min^{-1})')
boldify
legend('Pearson Correlation','Net Sensitivity','Purity','FAR','Location',...
    'southoutside','NumColumns',2)
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';
