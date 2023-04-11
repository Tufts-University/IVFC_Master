%% NV_033023_CorrelationPlots_DeepPeak
%% Initialization
clear
clc
mainpath = 'T:\Nilay\IVFC\DeepPeakResults';
cd(mainpath)
folders = dir('Trial_*');
[F_names] = natsortfiles({folders.name});
%% Plotting code
PearsonStorage = zeros(length(folders),7);
for i = 1:length(F_names)
    disp(['Plotting Trial # ', num2str(i),' of ',...
                        num2str(length(F_names))])
    cd(F_names{i})
    load('Pearson_All')
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
    PearsonStorage(i,1:3) = [str2double(F_names{i}(find(F_names{i}=='_')+1:end)),r(1,2),r(1,2).^2];
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
    load('Pearson')
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
    disp(['Calculating Trial # ', num2str(i),' of ',...
                        num2str(length(F_names))])
    cd(F_names{i})
    load('Pearson_FP')
    x = double(x);
    y = double(y);
    exclusion = [4,6];
    x_final = x;
%     x_final(exclusion) = [];
    total_FP = sum(x_final);
    total_time = 175; % mins
    flow_rate = 3; % uL/min
    total_volume = total_time.*flow_rate;
    FAR(i) = total_FP./total_volume;
    cd(mainpath)
end

