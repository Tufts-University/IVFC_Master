function successmessage=SamplePeakDetection_PCA_PN_TH_onlynormalization(filepath,...
    file_range,Window_Low,Window_High,Fs,sample_type,...
    exp_num,Spectralon_tail)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file:SamplePeakDetection.m
% ***Description***:
% This function serves automatically search for peaks in .mat files in 1.5
% minute chunks. The user must specify the sample type in order to ensure
% the correct thresholds are used. This algorithm is set up for mouse data,
% blood data, bead data, and cell data. Please read the input carefully as
% there are many ways to use this code!
% Written By: Nilay Vora (nvora01@tufts.edu)
% Date Written: 10/01/2021
% Modifying Author:Nilay Vora
% Date Modified: 01/13/2022
% Latest Revision: Added a new flag to seperate bead peaks from cell peaks
% in mixed blood samples.
% Modifying Author:Nilay Vora
% Date Modified: 06/06/2022
% Latest Revision: Changed the Green FLR FWHM measure script for peak
% equalization
% Date Modified: 08/21/2024 by Taras Hanulia (thanul01@tufts.edu)
% This function was change to remove the FWHM and peak area and added start
% and end of the peak extra sorting was added to have Scattering peak with
% FLR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function details
% Inputs:
%   filepath = Location of all main folder where all subfolders are
%   outputfile = Base name of all output files
%   file_range = Which files do you want to analyze in a specific day
%   Window_Low = Low frequency cutoff of Butterworth filter window
%   Window_High = High frequency cutoff of Butterworth filter window
%   Fs = Sample frequency used during acquisition
%   analysisvals = Which detection methods do you want to use
%   sample_type = (sting) Sample type (i.e. 'Cells','Blood','Beads',
%                  'Animal')
%   exp_num = Expected number of events in 1.5 minute chunk
%   std_threshold = Cluster cutoff value
%   Spectralon_tail= Tail format of Spectralon file
%   FWMH_threshold = Minimum peak width allowe do be detected
%   intensity_threshold= Minimum peak intensity for detection
%   bead_flag=Indicator that seperated beads from cell peaks in blood data
% Outputs:
%   successmessage = a string that indicates the completion of the code
%
% Usage SamplePeakDetection(filepath,outputfile,file_range,Window_Low,...
% Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
% Spectralon_tail,FWMH_threshold,intensity_threshold)
% Example:
% filepath = 'U:\Nilay\IVFC\Acquired Data\Blood Data\NV_092821_Blood_LNPs';
% outputfile= 'NEW_peak_values_09_28_21';
% file_range= [1:3];
% Window_Low= 50; % removes frequencies below 50/30000 Hz
% Window_High= 6000; % removes frequencies above 6000/30000 Hz
% Fs=60e3; %60,000 samples per second
% analysisvals=[1:2]; %For cells we may want [1:4] (Scat Only,FLR Only,...
%                       FLR+Scat,RFLR), for beads [1:3] (Scat Only,...
%                       FLR Only,FLR+Scat)
% sample_type= 'Beads';
% exp_num=[]; %used the default value of 6000/9
% std_threshold=3; 3*sigma(x) is used to detect a cluster event
% Spectralon_tail= '_1'; %in most cases there is a tail value, if there is
%                         none leave blank;
% FWMH_threshold=0; Usually kept at 0 and is inactive
% intensity_threshold= 0.1; Currently only used for FLR analysis!
% bead_flag=0;
% output=SamplePeakDetection(filepath,outputfile,file_range,Window_Low,...
% Window_High,Fs,analysisvals,sample_type,exp_num,std_threshold,...
% Spectralon_tail,FWMH_threshold,intensity_threshold,0)

%% Checking inputs
if isempty(sample_type)
    disp('Sample Type is not specified, please chose one of the options below');
    prompt = 'What Sample are you analyzing: Cells,Blood,Beads,Animal?';
    sample_type = input(prompt,'s');
end

if isempty(file_range)
    disp('File Range is not specified, please input range below');
    prompt = 'What is your File Range?';
    frange = input(prompt,'s');
    file_range=str2num(frange); %#ok<ST2NM>
end

if isempty(exp_num)
    exp_num=6000/9;
end

% if isempty(outputfile)
%     disp('Output File name is not specified, please input the name below');
%     prompt = 'What is your Output file name?';
%     outputfile = input(prompt,'s');
% end

if isempty(filepath)
    disp('Using Current Directory');
    filepath=pwd;
end
% if isempty(bead_flag) || strcmp(sample_type,'Blood')==0
%     bead_flag=0;
% else
%     bead_flag=1;
% end

%% Main Code
% Finding Folders with Files
mainFolder=filepath;
% fileName=outputfile;
cd(mainFolder)

dirinfo = dir();
dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories

dirinfo(ismember( {dirinfo.name}, {'.', '..'})) = [];  %remove . and ..
[~,c]=natsortfiles({dirinfo.name});
subdirinfo = cell(length(dirinfo));

for K = 1 : length(dirinfo)
    thisdir = dirinfo(c(K)).name;
    subdirinfo{K} = dir(fullfile(thisdir, '*_1_raw.mat'));
end

subdirinfo =  subdirinfo(~cellfun('isempty',subdirinfo));

Wn=[Window_Low Window_High]./(Fs/2);%Cutoff frequencies divided by Nyquist frequency
[b,a]=butter(2,Wn);
% Initialzing and pre-allocating

    for i=file_range
        disp(['Evaluating File # ',num2str(i),' of ',num2str(size(subdirinfo,1))]);
        tic
        if~isempty(subdirinfo{i})
            scatpath = subdirinfo{i}.folder;
            %spec_file =[subdirinfo{1}.name(1:end-4),'_Spectralon_Avg',Spectralon_tail,'.csv'];
            if subdirinfo{1}.name(end-10)=='l'
                spec_file =[subdirinfo{1}.name(1:end-10),'_Spectralon_Avg',Spectralon_tail,'.csv'];
            elseif subdirinfo{1}.name(end-11)=='l'
                spec_file =[subdirinfo{1}.name(1:end-11),'_Spectralon_Avg',Spectralon_tail,'.csv'];
            else 
                spec_file =[subdirinfo{1}.name(1:end-10),'_Spectralon_Avg',Spectralon_tail,'.csv'];
            end
            % Spectralon loading
            cd(mainFolder)
            folders = strsplit(mainFolder,'\');
            %folders(6) = [];
            folders = strjoin(folders,'\');
            spec_path = [folders,'\',spec_file];
            fid=fopen(spec_path);
            if fid<0
                folders = strsplit(mainFolder,'\');
                folders(6) = [];
                folders = strjoin(folders,'\');
                spec_path = [folders,'\',spec_file];
                fid=fopen(spec_path);
            end
            spec=textscan(fid,'%f %f %f %f %f %f %f %f','Delimiter',',');
            spec=cell2mat(spec);
            fclose(fid);

            cd(scatpath)
            data_type=subdirinfo{i}.name(1:end-10);
            cd(subdirinfo{i}.folder)
            dirinfo2 = dir('*_raw.mat');
            num_chunks=length(dirinfo2);

            % Load Standard Deviations
            load([folders,'\',data_type,'\',data_type,'_sigmas.mat'],'sigmas_final'); % variable is sigmas_final

           %% Determine expected cell count
            exp_num=floor(exp_num); %dimensionless number of cells.

            for ii=1:num_chunks
                %% Loading data
                disp([num2str(ii),' of ',num2str(num_chunks)])
                load([data_type,'_',num2str(ii),'_raw.mat'],'M') % variable =M %
                if isempty(M)

                else
                    %% Remove Clumps
                    cum=sum(M(:,1:3),2);
                    cum_std=std(cum);
                    stop = 0;
                    while 1
                        if cum_std<1.75 || stop == 3
                            break
                        end
                        for o=1:3
                            mean_clump=movmean(M(:,o),[250,250]);
                            %         clumps=zeros([length(M),1]);
                            %         clumps(mean_clump>5*std(mean_clump))=mean_clump(mean_clump>5*std(mean_clump));
                            M(:,o)=M(:,o)-mean_clump+mean(mean_clump);
                            %M(M(:,i)<-2*std(M(:,i)),i)=0;
                        end
                        disp('I ran')
                        cum=sum(M(:,1:3),2);
                        cum_std=std(cum);
                        stop = stop + 1;
                    end
                    %%%%%%%%%%%%%%%%%%%% 11/22 Edit %%%%%%%%%%%%%%%%%%%%%%%
                    %% Load Bead data if available:
                    % Check for bead file:
%                     bead_file = dir('*Beads.mat');
%                     if ~isempty(bead_file)
%                         M = bead_removal(ii, M, bead_file(1).name);
%                     end
                    M_filt = zeros(size(M)); 
                  
                    M_filt(:,1)=filtfilt(b,a,M(:,1))./abs(spec(1));
                    M_filt(:,2)=filtfilt(b,a,M(:,2))./abs(spec(2));
                    M_filt(:,3)=filtfilt(b,a,M(:,3))./abs(spec(3));
                    M_filt(:,4)=filtfilt(b,a,M(:,4));
                    M_filt(:,5)=filtfilt(b,a,M(:,5));

                    M_filt(:,1)=(M_filt(:,1)-mean(M_filt(:,1)));
                    M_filt(:,2)=(M_filt(:,2)-mean(M_filt(:,2)));
                    M_filt(:,3)=(M_filt(:,3)-mean(M_filt(:,3)));
                    M_filt(:,4)=(M_filt(:,4)-mean(M_filt(:,4)));
                    M_filt(:,5)=(M_filt(:,5)-mean(M_filt(:,5)));
                    clear M
                    %% Normalizing by standard deviation

                    SN_405=(M_filt(:,1)-mean(M_filt(:,1)))./sigmas_final(1);
                    SN_488=(M_filt(:,2)-mean(M_filt(:,2)))./sigmas_final(2);
                    SN_633=(M_filt(:,3)-mean(M_filt(:,3)))./sigmas_final(3);
                    SN_Red=(M_filt(:,4)-mean(M_filt(:,4)));
                    SN_Green=(M_filt(:,5)-mean(M_filt(:,5)));
                    SN_Green(SN_Green<0) = 0;
                    SN_Red(SN_Red<0) = 0;
                    M = [SN_405,SN_488,SN_633,SN_Red,SN_Green];
                    signal_max = 11;
                    %% Normalization
                    M2 = (M-min(M))./(max(M)-min(M));
                    %% Sensor cleaning
                    [~,~,~,tsquared,~,~] = pca(M2(:,1:3));
                    %% Normalization of the data set
                    M2(:,4:5) = [SN_Red,SN_Green]./signal_max;
                    norm = tsquared;
                    normalized_filename = fullfile(subdirinfo{i}.folder, [data_type, '_', num2str(ii), '_rawN.mat']);
                    save(normalized_filename, 'M2', '-v7.3');
                end
            end
        end
    end
successmessage='Completed Normalization';