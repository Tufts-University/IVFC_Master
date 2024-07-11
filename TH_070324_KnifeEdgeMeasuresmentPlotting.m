%% TH_07032024_KnifeEdgeMeasuresmentPlotting
%% Initialization
clear
clc
close all
main_path = 'T:\Taras\Knife Edge Measurments';
%% Load in Data [Slit Length]
cd(main_path)
T = readtable("TH_07032024_KnifeEdge Measurments_L");
T = sortrows(T,"Var1");
Position = T(:,1);
w_405 = T(:,2);
w_488 = T(:,3);
w_633 = T(:,4);
Position = table2array(Position);
%% Power Normalization
signal = [w_405,w_488,w_633];
signal = table2array(signal);
signal = signal - min(signal);
signal = signal./ max(signal);
%% Plot
colors=[0.49,0.18,0.56;...
        0.30,0.75,0.93;...
        0.64,0.08,0.18];
figure;
colororder(colors)
plot(Position,signal)
ylabel('Normalized Power')
xlabel('Position')
legend('405','488','633','Numcolumns',3,'Location','southoutside',...
        'AutoUpdate','off')
boldify
xlim([20,80])
xticks(20:2:80)
yline(0.9,'k--','LineWidth',3)
yline(0.1,'k--','LineWidth',3)
% xline(219,'k--','LineWidth',3)
% xline(250,'k--','LineWidth',3)
set(gcf,'Color','w')
export_fig('TH_07032024SlitLength.png')
%% Load in Data [Slit Width]
cd(main_path)
T2 = readtable("TH_07032024_KnifeEdge Measurments_W");
Position = T2(2:end,1);
w_405 = T2(2:end,2);
w_488 = T2(2:end,3);
w_633 = T2(2:end,4);
Position = table2array(Position);
%% Power Normalization
signal = [w_405,w_488,w_633];
signal = table2array(signal);
signal = signal - min(signal);
signal = signal./ max(signal);
%% Plot
colors=[0.49,0.18,0.56;...
        0.30,0.75,0.93;...
        0.64,0.08,0.18];
figure;
colororder(colors)
plot(Position,signal)
ylabel('Normalized Power')
xlabel('Position')
legend('405','488','633','Numcolumns',3,'Location','southoutside',...
        'AutoUpdate','off')
boldify
xlim([18,28])
xticks(18:1:28)
yline(0.9,'k--','LineWidth',3)
yline(0.1,'k--','LineWidth',3)
%xline(7.5,'k--','LineWidth',3)
%xline(12,'k--','LineWidth',3)
set(gcf,'Color','w')
export_fig('TH_07032024SlitWidth.png')
%% Load in Data [Slit Z]
cd(main_path)
T2 = readtable("TH_07032024_KnifeEdge Measurments_Z");
Position = T2(2:end,1);
w_405 = T2(2:end,2);
w_488 = T2(2:end,3);
w_633 = T2(2:end,4);
Position = table2array(Position);
%% Power Normalization
signal = [w_405,w_488,w_633];
signal = table2array(signal);
signal = signal - min(signal);
signal = signal./ max(signal);
%% Plot
colors=[0.49,0.18,0.56;...
        0.30,0.75,0.93;...
        0.64,0.08,0.18];
figure;
colororder(colors)
plot(Position,signal)
ylabel('Normalized Power')
xlabel('Position')
legend('405','488','633','Numcolumns',3,'Location','southoutside',...
        'AutoUpdate','off')
boldify
xlim([-250,60])
xticks(-250:10:60)
yline(0.9,'k--','LineWidth',3)
yline(0.1,'k--','LineWidth',3)
%xline(7.5,'k--','LineWidth',3)
%xline(12,'k--','LineWidth',3)
set(gcf,'Color','w')
export_fig('TH_07032024SlitWidth.png')