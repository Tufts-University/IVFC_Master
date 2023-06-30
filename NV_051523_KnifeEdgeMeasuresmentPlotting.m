%% NV_051523_KnifeEdgeMeasuresmentPlotting
%% Initialization
clear
clc
close all
main_path = 'T:\Nilay\IVFC\KnifeEdge Measurments';
%% Load in Data [Slit Length]
cd(main_path)
T = readtable("NV_050123_KnifeEdge Measurments.xlsx");
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
figure("Units","inches","Position",[3,3,3,2.5]);
colororder(colors)
plot(Position,signal)
ylabel('Normalized Power')
xlabel('Position')
legend('405','488','633','Numcolumns',3,'Location','southoutside',...
        'AutoUpdate','off')
boldify
xlim([150,300])
xticks(150:10:300)
yline(0.9,'k--','LineWidth',3)
yline(0.1,'k--','LineWidth',3)
xline(219,'k--','LineWidth',3)
xline(250,'k--','LineWidth',3)
set(gcf,'Color','w')
set(gca,'FontName','Arial','FontSize',18)
print('SlitLengthPlot','-dmeta')
%% Load in Data [Slit Width]
cd(main_path)
T2 = readtable("NV_050323_KnifeEdgeMeasurements_slitwidth.xlsx");
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
figure("Units","inches","Position",[3,3,3,2.5]);
colororder(colors)
plot(Position,signal)
ylabel('Normalized Power')
xlabel('Position')
legend('405','488','633','Numcolumns',3,'Location','southoutside',...
        'AutoUpdate','off')
boldify
xlim([0,20])
xticks(0:2:20)
yline(0.9,'k--','LineWidth',3)
yline(0.1,'k--','LineWidth',3)
xline(7.5,'k--','LineWidth',3)
xline(12,'k--','LineWidth',3)
set(gcf,'Color','w')
set(gca,'FontName','Arial','FontSize',18)
print('SlitWidthPlot','-dmeta')