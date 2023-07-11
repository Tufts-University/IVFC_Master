%% 03/10/23_DataAnalysis
% Owner: Nilay Vora
% Data Type: Full Dataset Threshold Performance
% Run Date: 03/10/2023
%% Initialization
clear
clc
%% Calling Script
%% Data Files Names
mainpath = 'T:\Nilay\IVFC\Acquired Data\Blood Cell Data\CNNData';
cd(mainpath)
[~,T] = xlsread('dict_key.csv');
T(1)=[];
threshold = [6,8,10,12,14];
performance = zeros(6,3);
performance(:,1) = [threshold,13]';
for j = 1:5
    thresh = threshold(j);
    disp(['Threshold ',num2str(j),' of 5'])
    for i=1:30%length(T)
        disp([9 'Reprocessing Day ',num2str(i),' of 30'])
        exp_name=T{i};
        cd(exp_name(2:end))
        filepath=cd;
        Fs=60e3;
        Window_Low= 50;
        Window_High= 10000;
        date=exp_name(5:10);
        date = textscan(date,'%2s');
        date = date{1};
        %% Peak Detection
        outputfile= ['Threshold_031023_peak_values_',date{1},'_',date{2},'_',date{3}];
        dirinfo = dir('*Cell*');
        dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
        dirinfo(ismember( {dirinfo.name}, {'.', '..'})) = [];  %remove . and ..
        file_range= (1:length(dirinfo));
        analysisvals = 1;
        sample_type= 'Blood';
        exp_num=round(6000/9);
        std_threshold=thresh;
        if i==21 || i==23
             Spectralon_tail= '_2';
        elseif i==22
             Spectralon_tail= '_3';
        else
             Spectralon_tail= '_1';
        end
        FWMH_threshold=0;
        intensity_threshold= 0.1;
        if i>=6 && i<=13
            bead_flag=0;
        else
            bead_flag=1;
        end
        output=SamplePeakDetection_PCA_PN(filepath,outputfile,file_range,Window_Low,...
            Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
            Spectralon_tail,FWMH_threshold,intensity_threshold,bead_flag);
        disp(output)
        cd(mainpath)
    end
    header = 'Threshold_031023';
    [sens,pur] = ScatvFLRCalculator(mainpath,header);
    performance(j,2:3) = [sens,pur];
    cd(mainpath)
end
header = 'Corrected_012923';
[sens,pur] = ScatvFLRCalculator(mainpath,header);
performance(6,2:3) = [sens,pur];
performance_final = sortrows(performance,1);

%% Ploting of data
figure;
x = performance_final(:,1);
sens_final = performance_final(:,2);
pur_final = performance_final(:,3);
yyaxis('left')
plot(x,sens_final,'.','MarkerSize',23)
ylabel('Sensitivity Performance')
yyaxis("right")
plot(x,pur_final,'.','MarkerSize',23)
legend('Sensitivity','Purity','NumColumns',2,'Location','southoutside')
xlabel('Threshold')
ylabel('Purity Performance')
xticks([6:2:14])
xlim([5,15])

boldify