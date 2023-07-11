function [successmessage]=DailyCalibrationScript(files)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  file: DailyCalibrationScript.m
%  This file will read the outputted peak_values data files and
%  generate histograms and conduct ANOVA calculation for verification of
%  alignment from day to day
%  Written By: Nilay Vora (nvora01@tufts.edu)
%  Date Written: 01/12/2022
%  Modifying Author:
%  Date Modified:
%  Latest Revision: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function details
% Inputs:
%   files = Name of all peak value files to read
% Outputs:
%   successmessage = a string that indicates the completion of the code
%
% Usage DailyCalibrationScript(files)
% Example:
% files={'NEW_peak_values_11_22_21';...
%     'NEW_peak_values_11_29_21';...
%     'NEW_peak_values_11_30_21';...
%     'NEW_peak_values_12_14_21';...
%     'NEW_peak_values_12_20_21';...
%     'NEW_peak_values_01_12_22'};
% output= DailyCalibrationScript(files)
%% Initialization
close all hidden;

%% Check to see if all files are in Path
currentpath=pwd;
if ismac
    myFolders = split(currentpath,"/");
    myFolders(end)=[];
    myFolders(end)=[];
    myNewPath = join(myFolders,"/");
elseif ispc
    myFolders = split(currentpath,"\");
    myFolders(end)=[];
    myFolders(end)=[];
    myNewPath = join(myFolders,"\");
end
cd(myNewPath{1})
dirinfo = dir();
dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
dirinfo(ismember( {dirinfo.name}, {'.', '..'})) = [];  %remove . and ..
for k=1:size(dirinfo,1)
    if ispc
        folders=[dirinfo(k).folder,'\',dirinfo(k).name];
        addpath(genpath(folders))
    elseif ismac
        folders=[dirinfo(k).folder,'/',dirinfo(k).name];
        addpath(genpath(folders))
    end
end
cd(currentpath)
% List all the files you want to analyze and compare
num_files=size(files,1); % number of files listed above

% List the legend elements
legend_el=cell(1,size(files,1));
for j=1:size(files,1)
    v=str2double(strip(strsplit(files{j},'_'),'m'));
    v=v(~isnan(v));
    legend_el{j}=[num2str(v(1)),'/',num2str(v(2))];
end
%% Main Loop
colors={'405','488','633','Red FLR','Green FLR'};
for i=1:5
    figure;
    max_int=0;
    dataset=cell(1,num_files);
    for ii=1:num_files % loads each channel of data from all the files
        data= load(files{ii});
        dataset{1 ,ii}=data.peak_values(:,i);
        if max_int<max(dataset{1,ii})
            max_int=max(dataset{1,ii});
        end
    end
    %average(i,:)=[mean(dataset{1}),mean(dataset{2}),mean(dataset{3})];
    fits = zeros(100,2,numel(dataset));
    edges=linspace(0,max_int,20);
    hold on
    for k = 1:num_files
        total = numel(dataset{k}); % for normalizing
        %f = histfit(dataset{k},20,'kernel'); % draw the histogram and fit
        pd = fitdist(dataset{k},'kernel');
        q = icdf(pd,[0.0013499 0.99865]); % three-sigma range for normal distribution
        x = linspace(q(1),q(2)); 
        [bincounts,binedges] = histcounts(dataset{k}, edges);
        bincenters = binedges(1:end-1)+diff(binedges)/2;
        hh = bar(bincenters,bincounts,1);
        binwidth = binedges(2)-binedges(1); % Finds the width of each bin
        area = total * binwidth;
        y = area * pdf(pd,x);
        % Overlay the density
        hh1 = plot(x,y,'r-','LineWidth',2);
        f = [hh; hh1];

        % collect the curve data and normalize it:
        fits(:,:,k) = [f(2).XData; f(2).YData./total].';
        x = f(1).XData; % collect the bar positions
        n = f(1).YData; % collect the bar counts
        f.delete % delete the histogram and the fit
        bar(x,n./total,'histc'); % plot the bar
    end
    ax = gca; % get the axis handle
    % set all color and transparency for the bars:
    set(ax.Children,{'FaceColor'},mat2cell(lines(num_files),ones(num_files,1)))
    set(ax.Children,{'FaceAlpha'},repmat({0.5},num_files,1))
    % plot all the curves:
    plot(squeeze(fits(:,1,:)),squeeze(fits(:,2,:)),'LineWidth',3)
    hold off

    colorc=lines(num_files); % Establishes line colors to match bar colors.
    for iv=1:numel(dataset)
        ax.Children(iv).Color=colorc(iv,:);
    end


    title(strjoin([colors(i),'Data']))

    xlabel('Intensity')
    ylabel('Counts')
    legend(legend_el)

    fit_data=[];
    group=[];
    for iii=1:num_files % Sets up data for ANOVA
        fit_data=[fit_data;fits(:,2,iii).*fits(:,1,iii)]; %#ok<AGROW> 
        group=[group,repmat(legend_el(:,iii),1,100)]; %#ok<AGROW> 
    end

    group={group};
    [~,~,stats] = anovan(fit_data,group,'varnames',{'Date'});
    figure;
    [~,~,~,~] = multcompare(stats);
end

%% SAVE Figs
tempdir=cd;
cd(tempdir);
if ~exist([cd,'\Data_Plots'],'dir')
    mkdir Data_Plots
end
tempdir=[tempdir,'\Data_Plots'];
FolderName = tempdir;   % Your destination folder
for iFig = 1:15
    FigHandle = figure(iFig);
    FigName   = ['Stats_Plots(norm)_',num2str(iFig)];
    saveas(FigHandle, fullfile(FolderName, [FigName, '.png']));
end
successmessage='Daily Calibration Plots Saved';