function successmessage=SamplePeakDetection6PMTRF(filepath,outputfile,...
    file_range,Window_Low,Window_High,Fs,analysisvals,sample_type,...
    exp_num,std_threshold,Spectralon_tail,FWMH_threshold,...
    intensity_threshold,bead_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file:SamplePeakDetection6PMT.m
% ***Description***:
% This function serves automatically search for peaks in .mat files in 1.5
% minute chunks. The user must specify the sample type in order to ensure
% the correct thresholds are used. This algorithm is set up for mouse data,
% blood data, bead data, and cell data. Please read the input carefully as
% there are many ways to use this code! It is apdated for 6PMT(3 type of
% fluorescence (iRFP720, green, E2Crimson)
% Written By: Taras Hanulia (thanul01@tufts.edu)
% Date Written: 15/05/2023
% Modifying Author:Taras Hanulia
% Date Modified: 11/14/2024
% Latest Revision:peak sorting added for iRFP and GFP channel , Sorting for
% gingle iRFP and single GFP added
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
% analysisvals=[1:2]; %For cells we may want [1,2, 4,5,6,7] (Scat Only,FLR Only,...
%                       iRFP,E2C,RFLR,GFP), for beads [1,2,3] (Scat Only,...
%                       FLR Only,FLR+Scat)
%iRFP only signal from irfpchannel . E2C only signal from E2C channel, ...
%RFL from iRFP+E2C  , GFP only from GFP channel
% sample_type= 'Beads';
% exp_num=[]; %used the default value of 6000/9
% std_threshold=3; 3*sigma(x) is used to detect a cluster event
% Spectralon_tail= '_1'; %in most cases there is a tail value, if there is
%                         none leave blank;
% FWMH_threshold=0; Usually kept at 0 and is inactive
% intensity_threshold= 0.15; Currently only used for FLR analysis!
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

if isempty(outputfile)
    disp('Output File name is not specified, please input the name below');
    prompt = 'What is your Output file name?';
    outputfile = input(prompt,'s');
end

if isempty(filepath)
    disp('Using Current Directory');
    filepath=pwd;
end
if isempty(bead_flag) || strcmp(sample_type,'Blood')==0 
    bead_flag=0;
else
    bead_flag=1;
end
%% Modifying Evaluation Parameters
switch sample_type
    case 'Cells'
        if nargin==13
        else
            if nargin < 13
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 12
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 11
                Spectralon_tail= '_1'; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 10
                std_threshold=4; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= '_1'; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 9
                exp_num=6000/9; %Expected number of clusters
                std_threshold=4; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= '_1'; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 7
                analysisvals=(1:3); % Assumes Cells in blood with beads
                exp_num=6000/9; %Expected number of clusters
                std_threshold=4; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= '_1'; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 6
                Fs=60e3; %60,000 samples per second
                analysisvals=(1:3); % Assumes Cells in blood with beads
                exp_num=6000/9; %Expected number of clusters
                std_threshold=4; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= '_1'; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 5
                Window_High= 6000; % removes frequencies above 6000/30000 Hz
                Fs=60e3; %60,000 samples per second
                analysisvals=(1:3); % Assumes Cells in blood with beads
                exp_num=6000/9; %Expected number of clusters
                std_threshold=4; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= '_1'; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 4
                Window_Low= 50; % removes frequencies below 50/30000 Hz
                Window_High= 6000; % removes frequencies above 6000/30000 Hz
                Fs=60e3; %60,000 samples per second
                analysisvals=(1:3); % Assumes Cells in blood with beads
                exp_num=6000/9; %Expected number of clusters
                std_threshold=4; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= '_1'; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            end
        end
    case 'Blood'
        if nargin==13
        else
            if nargin < 13
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 12
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 11
                Spectralon_tail= '_1'; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 10
                std_threshold=3; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= '_1'; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 9
                exp_num=6000/9; %Expected number of clusters
                std_threshold=3; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= '_1'; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 7
                analysisvals=(1:4); % Assumes Cells in blood with beads
                exp_num=6000/9; %Expected number of clusters
                std_threshold=3; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= '_1'; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 6
                Fs=60e3; %60,000 samples per second
                analysisvals=(1:4); % Assumes Cells in blood with beads
                exp_num=6000/9; %Expected number of clusters
                std_threshold=3; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= '_1'; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 5
                Window_High= 6000; % removes frequencies above 6000/30000 Hz
                Fs=60e3; %60,000 samples per second
                analysisvals=(1:4); % Assumes Cells in blood with beads
                exp_num=6000/9; %Expected number of clusters
                std_threshold=3; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= '_1'; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 4
                Window_Low= 50; % removes frequencies below 50/30000 Hz
                Window_High= 6000; % removes frequencies above 6000/30000 Hz
                Fs=60e3; %60,000 samples per second
                analysisvals=(1:4); % Assumes Cells in blood with beads
                exp_num=6000/9; %Expected number of clusters
                std_threshold=3; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= '_1'; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            end
        end
    case 'Beads'
        if nargin==13
        else
            if nargin < 13
                intensity_threshold= 0.25; %Currently only used for FLR analysis!
            elseif nargin < 12
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.25; %Currently only used for FLR analysis!
            elseif nargin < 11
                Spectralon_tail= ''; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.25; %Currently only used for FLR analysis!
            elseif nargin < 10
                std_threshold=4; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= ''; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.25; %Currently only used for FLR analysis!
            elseif nargin < 9
                exp_num=6000/9; %Expected number of clusters
                std_threshold=4; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= ''; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.25; %Currently only used for FLR analysis!
            elseif nargin < 7
                analysisvals=(1:3); % Assumes Cells in blood with beads
                exp_num=6000/9; %Expected number of clusters
                std_threshold=4; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= ''; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.25; %Currently only used for FLR analysis!
            elseif nargin < 6
                Fs=60e3; %60,000 samples per second
                analysisvals=(1:3); % Assumes Cells in blood with beads
                exp_num=6000/9; %Expected number of clusters
                std_threshold=4; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= ''; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.25; %Currently only used for FLR analysis!
            elseif nargin < 5
                Window_High= 6000; % removes frequencies above 6000/30000 Hz
                Fs=60e3; %60,000 samples per second
                analysisvals=(1:3); % Assumes Cells in blood with beads
                exp_num=6000/9; %Expected number of clusters
                std_threshold=4; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= ''; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.25; %Currently only used for FLR analysis!
            elseif nargin < 4
                Window_Low= 50; % removes frequencies below 50/30000 Hz
                Window_High= 6000; % removes frequencies above 6000/30000 Hz
                Fs=60e3; %60,000 samples per second
                analysisvals=(1:3); % Assumes Cells in blood with beads
                exp_num=6000/9; %Expected number of clusters
                std_threshold=4; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= ''; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.25; %Currently only used for FLR analysis!
            end
        end
    case 'Animal'
        if nargin==13
        else
            if nargin < 13
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 12
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 11
                Spectralon_tail= '_1'; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 10
                std_threshold=5; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= '_1'; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 9
                exp_num=6000/9; %Expected number of clusters
                std_threshold=5; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= '_1'; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 7
                analysisvals=(1:4); % Assumes Cells in blood with beads
                exp_num=6000/9; %Expected number of clusters
                std_threshold=5; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= '_1'; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 6
                Fs=60e3; %60,000 samples per second
                analysisvals=(1:4); % Assumes Cells in blood with beads
                exp_num=6000/9; %Expected number of clusters
                std_threshold=5; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= '_1'; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 5
                Window_High= 6000; % removes frequencies above 6000/30000 Hz
                Fs=60e3; %60,000 samples per second
                analysisvals=(1:4); % Assumes Cells in blood with beads
                exp_num=6000/9; %Expected number of clusters
                std_threshold=5; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= '_1'; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            elseif nargin < 4
                Window_Low= 50; % removes frequencies below 50/30000 Hz
                Window_High= 6000; % removes frequencies above 6000/30000 Hz
                Fs=60e3; %60,000 samples per second
                analysisvals=(1:4); % Assumes Cells in blood with beads
                exp_num=6000/9; %Expected number of clusters
                std_threshold=5; %3*sigma(x) is used to detect a cluster event
                Spectralon_tail= '_1'; %in most cases there is a tail value, if there
                % is none leave blank;
                FWMH_threshold=0; %Usually kept at 0 and is inactive
                intensity_threshold= 0.1; %Currently only used for FLR analysis!
            end
        end
end
%% Main Code
% Finding Folders with Files
mainFolder=filepath;
fileName=outputfile;
cd(mainFolder)

dirinfo = dir();
dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories

dirinfo(ismember( {dirinfo.name}, {'.', '..'})) = [];  %remove . and ..
[~,c]=natsortfiles({dirinfo.name});
subdirinfo = cell(length(dirinfo));

for K = 1 : length(dirinfo)
    thisdir = dirinfo(c(K)).name;
    subdirinfo{K} = dir(fullfile(thisdir, '*.csv'));
end

subdirinfo =  subdirinfo(~cellfun('isempty',subdirinfo));

Wn=[Window_Low Window_High]./(Fs/2);%Cutoff frequencies divided by Nyquist frequency
[b,a]=butter(2,Wn);
types={'Scattering Only';'FLR Only';'Cumulative';'iRFP';'E2C';'RFLR';'GFLR'};
% Initialzing and pre-allocating
for f=analysisvals
    all_peaks=zeros(1,1);
    all_locs=zeros(1,1);
    widths_cum=zeros(1,1);
    peak_area_cum=zeros(1,1);
    all_chunks=zeros(1,1);
    peak_values=zeros(1,6);
    all_fwhm_405=zeros(1,1); all_fwhm_488=zeros(1,1); all_fwhm_633=zeros(1,1); all_fwhm_fl1=zeros(1,1); all_fwhm_fl2=zeros(1,1);all_fwhm_fl3=zeros(1,1);
    all_peak_area_405=zeros(1,1); all_peak_area_488=zeros(1,1); all_peak_area_633=zeros(1,1); all_peak_area_fl1=zeros(1,1); all_peak_area_fl2=zeros(1,1);all_peak_area_fl3=zeros(1,1);
    all_file_num=zeros(1,1);

    all_peaks_Store=zeros(1,1);
    all_locs_Store=zeros(1,1);
    widths_cum_Store=zeros(1,1);
    peak_area_cum_Store=zeros(1,1);
    all_chunks_Store=zeros(1,1);
    peak_values_Store=zeros(1,6);
    all_fwhm_405_Store=zeros(1,1); all_fwhm_488_Store=zeros(1,1); all_fwhm_633_Store=zeros(1,1); all_fwhm_fl1_Store=zeros(1,1); all_fwhm_fl2_Store=zeros(1,1);all_fwhm_fl3_Store=zeros(1,1);
    all_peak_area_405_Store=zeros(1,1); all_peak_area_488_Store=zeros(1,1); all_peak_area_633_Store=zeros(1,1); all_peak_area_fl1_Store=zeros(1,1); all_peak_area_fl2_Store=zeros(1,1);all_peak_area_fl3_Store=zeros(1,1);
    all_file_num_Store=zeros(1,1);

    disp(['Current evaluation based on ',types{f},' Data']);
    for i=file_range
        disp(['Evaluating File # ',num2str(i),' of ',num2str(size(subdirinfo,1))]);
        tic
        if~isempty(subdirinfo{i})
            scatpath=subdirinfo{i}.folder;
            spec_file =[subdirinfo{1}.name(1:end-4),'_Spectralon_Avg',Spectralon_tail,'.csv'];
            % Spectralon loading
            cd(mainFolder)
            fid=fopen(spec_file);
            spec=textscan(fid,'%f %f %f %f %f %f %f %f %f','Delimiter',',');
            spec=cell2mat(spec);
            fclose(fid);
            % Find all raw files
            cd(scatpath)
            data_type=subdirinfo{i}.name(1:end-4);
            dirinfo2 = dir();
            dirinfo2(~[dirinfo2.isdir]) = [];  %remove non-directories
            subdirinfo2 = cell(length(dirinfo2));
            for K = 1 : length(dirinfo2)
                thisdir2 = dirinfo2(K).name;
                subdirinfo2{K} = dir(fullfile(thisdir2, '*_raw.mat'));
            end
            num_chunks=length(subdirinfo2{1});
            % Select detection scheme
            if f==1
                flr_detect_1 = 0; % iRFP red Fluorescence
                flr_detect_2 = 0; % Green Fluorescence
                flr_detect_3 = 0; % E2C red Fluorescence
                fileN=[fileName,'_NoFLR'];
            elseif f==2
                flr_detect_1 = 0; % iRFP red Fluorescence
                flr_detect_2 = 1; % Green Fluorescence
                flr_detect_3 = 0; % E2C red Fluorescence
                if strcmp('Blood',sample_type) && bead_flag==1 
                    fileN=[fileName,'_NoScatAll'];
                else
                    fileN=[fileName,'_NoScat'];
                end
            elseif f==3
                flr_detect_1 = 1; % iRFP red Fluorescence
                flr_detect_2 = 1; % Green Fluorescence
                flr_detect_3 = 1; % E2C red Fluorescence
                fileN=[fileName]; %#ok<NBRAK>
            elseif f==4
                flr_detect_1 = 1; % iRFP red Fluorescence
                flr_detect_2 = 0; % Green Fluorescence
                flr_detect_3 = 0; % E2C red Fluorescence
                if strcmp('Blood',sample_type) && bead_flag==1 
                    fileN=[fileName,'_iRFPFLRAll'];
                else
                    fileN=[fileName,'_iRFPFLR']; 
                end
            elseif f==5
                flr_detect_1 = 0; % iRFP red Fluorescence
                flr_detect_2 = 0; % Green Fluorescence
                flr_detect_3 = 1; % E2C red Fluorescence
                fileN=[fileName,'_E2CPFLR']; 
            elseif f==7
                flr_detect_1 = 0; % iRFP red Fluorescence
                flr_detect_2 = 1; % Green Fluorescence
                flr_detect_3 = 1; % E2C red Fluorescence
                if strcmp('Blood',sample_type) && bead_flag==1 
                    fileN=[fileName,'_GFPFLRAll'];
                else
                    fileN=[fileName,'_GFPFLR']; 
                end
            else
                flr_detect_1 = 1; % iRFP red Fluorescence
                flr_detect_2 = 0; % Green Fluorescence
                flr_detect_3 = 1; % E2C red Fluorescence
                fileN=[fileName,'_RFLR'];
            end
            % Load Standard Deviations
            load([data_type,'_sigmas.mat'],'sigmas_final'); % variable is sigmas_final
            %% Determine expected cell count
            exp_num=floor(exp_num); %dimensionless number of cells.

            for ii=1:num_chunks
                %% Loading data
                disp([num2str(ii),' of ',num2str(num_chunks)])
                load([data_type,'_',num2str(ii),'_raw.mat'],'M') % variable =M %
                if isempty(M)

                else
                    %% Filtering %%
                    M_filt(:,1)=filtfilt(b,a,M(:,1))./abs(spec(1));
                    M_filt(:,2)=filtfilt(b,a,M(:,2))./abs(spec(2));
                    M_filt(:,3)=filtfilt(b,a,M(:,3))./abs(spec(3));
                    M_filt(:,4)=filtfilt(b,a,M(:,4));
                    M_filt(:,5)=filtfilt(b,a,M(:,5));
                    M_filt(:,6)=filtfilt(b,a,M(:,6));


                    M_filt(:,1)=(M_filt(:,1)-mean(M_filt(:,1)));
                    M_filt(:,2)=(M_filt(:,2)-mean(M_filt(:,2)));
                    M_filt(:,3)=(M_filt(:,3)-mean(M_filt(:,3)));
                    M_filt(:,4)=(M_filt(:,4)-mean(M_filt(:,4)));
                    M_filt(:,5)=(M_filt(:,5)-mean(M_filt(:,5)));
                    M_filt(:,6)=(M_filt(:,6)-mean(M_filt(:,6)));
                    clear M
                    %% Normalizing by standard deviation

                    SN_405=(M_filt(:,1)-mean(M_filt(:,1)))./sigmas_final(1);
                    SN_488=(M_filt(:,2)-mean(M_filt(:,2)))./sigmas_final(2);
                    SN_633=(M_filt(:,3)-mean(M_filt(:,3)))./sigmas_final(3);
                    SN_Red1=(M_filt(:,4)-mean(M_filt(:,4)));
                    SN_Green=(M_filt(:,5)-mean(M_filt(:,5)));
                    SN_Red2=(M_filt(:,6)-mean(M_filt(:,6)));
                    Norm=1;
                    SN_Green=SN_Green.*Norm;
                    SN_Red1=SN_Red1.*Norm;
                    SN_Red2=SN_Red2.*Norm;
                    
                    

                    if flr_detect_1==1 && flr_detect_2==1 && flr_detect_3==1 %ALL Channels
                        cumulative_det=SN_405+SN_488+SN_633+SN_Green;
                        cumulative=SN_405+SN_488+SN_633;
                    elseif flr_detect_1==1 && flr_detect_2==0 && flr_detect_3==1 % red FL (FL1 and FL3)
                        cumulative_det=SN_Red1+SN_Red2;
                        cumulative=SN_405+SN_488+SN_633;
                    elseif flr_detect_1==0 && flr_detect_2==1 && flr_detect_3==0 %FL1 +FL2 + FL3
                        cumulative_det=SN_Green+SN_Red1+SN_Red2;
                        cumulative=SN_405+SN_488+SN_633;
                    elseif flr_detect_1==1 && flr_detect_2==0 && flr_detect_3==0 % FL1 iRFP
                        cumulative_det=SN_Red1;
                        cumulative=SN_405+SN_488+SN_633;
                    elseif flr_detect_1==0 && flr_detect_2==0 && flr_detect_3==1 % FL3 E2C 
                        cumulative_det=SN_Red2;
                        cumulative=SN_405+SN_488+SN_633;
                    elseif flr_detect_1==0 && flr_detect_2==1 && flr_detect_3==1 %FL2 GFP
                        cumulative_det=SN_Green;
                        cumulative=SN_405+SN_488+SN_633;                       
                    else % ONLY SCATTERING
                        cumulative=SN_405+SN_488+SN_633;
                        cumulative_det=SN_405+SN_488+SN_633;
                    end
                    if f==1||f==3
                        for row=1:length(cumulative)
                            if cumulative(row)<0
                                cumulative(row)=0;
                                cumulative_det(row)=0;
                            end
                        end
                    end

                    cum_std(ii) = std(cumulative); %#ok<AGROW>
                    cum_avg(ii) = mean(cumulative); %#ok<AGROW>
                    cum_det_std(ii)=std(cumulative_det);%#ok<AGROW>
                    %% Peak Finding %%
                    %MINPEAKDISTANCE= specified how many points past the peak the next peak must be. Based on observations typical peak is 5 so choose 10
                    %NPEAKS= 3 times the expected number per chunk.  This is to protect against possibility that expected number varies from chunk to chunk
                    %
                    %[~,locs]=findpeaks(cumulative_det,'NPeaks',10*exp_num,'SortStr','descend');
                    if flr_detect_1==0 && flr_detect_2==1 && flr_detect_3==0 
                        [~,locs]=findpeaks(cumulative_det,'NPeaks',5*exp_num,'SortStr','descend');
                    else
                        [~,locs]=findpeaks(cumulative_det,'NPeaks',20*exp_num,'SortStr','descend');
                    end
                    locs=sort(locs,'ascend');
                    peaks=cumulative(locs);
                    peaks_det=cumulative_det(locs);
                    %% Clips peaks near edges
                    dist_left_edge=find(locs-100<101);
                    if length(dist_left_edge)>=1
                        peaks(dist_left_edge)=[];
                        peaks_det(dist_left_edge)=[];
                        locs(dist_left_edge)=[];
                    end

                    dist_right_edge=find(length(cumulative)-locs<101);
                    if length(dist_right_edge)>=1
                        peaks(dist_right_edge)=[];
                        peaks_det(dist_right_edge)=[];
                        locs(dist_right_edge)=[];
                    end
                    clear dist_left_edge dist_right_edge

                    %% Catch Clusters(round 1)
                    endpt1=find(cumulative_det(1:end-1)>1*cum_det_std(ii) & cumulative_det(2:end) < 1*cum_det_std(ii));
                    startpt1=find(cumulative_det(1:end-1)<1*cum_det_std(ii) & cumulative_det(2:end)> 1*cum_det_std(ii));
                    if length(startpt1)==length(endpt1)
                        ranges=[startpt1,endpt1];
                    elseif length(startpt1)>length(endpt1)
                        if startpt1(end)>endpt1(end)
                            startpt1(end)=[];
                        end
                        ranges=[startpt1,endpt1];
                    elseif length(startpt1)<length(endpt1)
                        if startpt1(1)>endpt1(1)
                            endpt1(1)=[];
                        end
                        ranges=[startpt1,endpt1];
                    end
                    save_locs=[];
                    for m=1:size(ranges,1)
                        if ranges(m,2)<ranges(m,1)
                            n=m+1;
                        else
                            n=m;
                        end
                        if n>size(ranges,1)

                        else
                            pt1=ranges(m,1);
                            pt2=ranges(n,2);
                            locs_clusters=locs(locs>=pt1 & locs<=pt2);
                            if ~isempty(locs_clusters)
                                save_locs=[save_locs;locs_clusters]; %#ok<AGROW>
                            end
                        end
                    end
                    p_store=[];
                    p_store2=[];
                    l_store=[];
                    for n=1:length(save_locs)
                        ind=find(locs==save_locs(n));
                        p_store=[p_store;peaks(ind)]; %#ok<AGROW>
                        p_store2=[p_store2;peaks_det(ind)]; %#ok<AGROW>
                        l_store=[l_store;locs(ind)]; %#ok<AGROW>
                    end
                    peaks=p_store;
                    peaks_det=p_store2;
                    locs=l_store;

                    clear p_store l_store save_lcos locs_clusters pt1 pt2 ranges pstore2
                    %% Catch Clusters
                    endpt=find(cumulative_det(1:end-1)>std_threshold*cum_det_std(ii) & cumulative_det(2:end) < std_threshold*cum_det_std(ii));
                    startpt=find(cumulative_det(1:end-1)<std_threshold*cum_det_std(ii) & cumulative_det(2:end)> std_threshold*cum_det_std(ii));

                    if length(startpt)==length(endpt)
                        ranges=[startpt,endpt];
                    elseif length(startpt)>length(endpt)
                        if startpt(end)>endpt(end)
                            startpt(end)=[];
                        end
                        ranges=[startpt,endpt];
                    elseif length(startpt)<length(endpt)
                        if startpt(1)>endpt(1)
                            endpt(1)=[];
                        end
                        ranges=[startpt,endpt];
                    end
                    tossed_peak=0;
                    clusterpts=[];
                    save_locs=[];
                    for m=1:size(ranges,1)
                        if ranges(m,2)<ranges(m,1)
                            n=m+1;
                        else
                            n=m;
                        end
                        if n>size(ranges,1)

                        else
                            pt1=ranges(m,1);
                            pt2=ranges(n,2);
                            locs_clusters=locs(locs>=pt1 & locs<=pt2);
                            if ~isempty(locs_clusters)
                                [~,I] = max(cumulative(locs_clusters));
                                save_locs=[save_locs;locs_clusters(I)]; %#ok<AGROW>
                                locs_clusters(I)=[];
                                tossed_peak=[tossed_peak;locs_clusters]; %#ok<AGROW>
                                clusterpts=[clusterpts;[pt1,pt2]]; %#ok<AGROW>
                            end
                        end
                    end
                    p_store=[];
                    l_store=[];
                    p_store2=[];
                    for n=1:length(save_locs)
                        ind=find(locs==save_locs(n));
                        p_store=[p_store;peaks(ind)]; %#ok<AGROW>
                        p_store2=[p_store2;peaks_det(ind)]; %#ok<AGROW>
                        l_store=[l_store;locs(ind)]; %#ok<AGROW>
                    end
                    peaks=p_store;
                    locs=l_store;

                    clear p_store l_store
                    %% Throws out peaks that are too close
                    % Peaks are sorted descending so that peaks that are thrown away start
                    % from the largest peak I.E. the most likely peaks
                    bad_peaks=0; %initializing
                    for k=1:length(locs)
                        peak_dist=abs(locs-locs(k)); %distances between all peaks and ith peak
                        close_peaks=find(peak_dist<=10); %distance is 5 pts between peaks
                        if length(close_peaks)>1 %peak will identify itself as less than 5 pts away
                            [~,idx]=sort(peaks(close_peaks),'descend');
                            cp_sorted=close_peaks(idx);
                            for j=2:length(close_peaks)
                                bad_peak=cp_sorted(j);
                                bad_peaks=[bad_peaks,bad_peak]; %#ok<AGROW>
                            end
                            clear bad_peak
                        end
                        clear close_peaks peak_dist
                    end
                    bad_peaks(1)=[]; % removing initial 0
                    peaks(bad_peaks)=[];
                    locs(bad_peaks)=[];
                    clusterpts(bad_peaks,:)=[]; %#ok<AGROW>
                    clear bad_peaks % reduces memory footprint

                    %% FWHM finding & Grabbing True Maximums
                    %FWHMs are found from the filtered peaks
                    fwhm=zeros(length(peaks),1);
                    peak_area_405=zeros(length(peaks),1); peak_area_488=zeros(length(peaks),1); peak_area_633=zeros(length(peaks),1); peak_area_fl1=zeros(length(peaks),1); peak_area_fl2=zeros(length(peaks),1); peak_area_fl3=zeros(length(peaks),1);
                    fwhm_405=zeros(length(peaks),1); fwhm_488=zeros(length(peaks),1); fwhm_633=zeros(length(peaks),1); fwhm_fl1=zeros(length(peaks),1); fwhm_fl2=zeros(length(peaks),1);fwhm_fl3=zeros(length(peaks),1);
                    peak_data=zeros(length(locs),6);

                    for m=1:length(peaks)
                        %% Grab True Maximums
                        data_range=M_filt(locs(m)-10:locs(m)+10,:); %grabs full cluster width by 5 matrix around peak
                        peak_data(m,:)=max(data_range); %grabs maximum which could be slightly different from channel to channel in time
                        clear data_range
                        %% Grab FWHMs
                        peak_height=M_filt(locs(m),:);% locs(m)
                        if clusterpts(m,2)+10>5400000
                            data_fwhm=M_filt(clusterpts(m,1)-10:end,:);
                            data_fwhm_cum=cumulative(clusterpts(m,1)-10:end);
                        else
                            data_fwhm=M_filt(clusterpts(m,1)-10:clusterpts(m,2)+10,:);
                            data_fwhm_cum=cumulative(clusterpts(m,1)-10:clusterpts(m,2)+10);
                        end
                        peak_height_cum=cumulative(locs(m)); % locs(m)
                        %405
                        [fwhm_405(m),peak_area_405(m)]=NV_101719_fwhm_measure(data_fwhm(:,1),peak_height(1));
                        %488
                        [fwhm_488(m),peak_area_488(m)]=NV_101719_fwhm_measure(data_fwhm(:,2),peak_height(2));
                        %633
                        [fwhm_633(m),peak_area_633(m)]=NV_101719_fwhm_measure(data_fwhm(:,3),peak_height(3));
                        %Fl1
                        [fwhm_fl1(m),peak_area_fl1(m)]=NV_101719_fwhm_measure(data_fwhm(:,4),peak_height(4));
                        %Fl2
                        [fwhm_fl2(m),peak_area_fl2(m)]=NV_052322_fwhm_measure(data_fwhm(:,5),peak_height(5));
                        %Fl3
                        [fwhm_fl3(m),peak_area_fl3(m)]=NV_101719_fwhm_measure(data_fwhm(:,6),peak_height(6)); % we dont have NV_052322_fwhm_measure for 6 
                        %cum
                        [fwhm(m),peak_area(m)]=NV_101719_fwhm_measure(data_fwhm_cum,peak_height_cum); %#ok<AGROW>
                        % Get rid of zeros
                    end
                    %% Preparing for next iteration of loop
                    if sample_type=='Beads' %#ok<BDSCA> 
                         idx=find(fwhm>0 & sum(peak_data,2)>1.5);
                    else
                        idx=find(fwhm>FWMH_threshold);% & peaks>0); % Was 4.0 on blood 2.5 for beads
                    end
                    if isempty(idx)==0
                        channel_map = containers.Map({2, 4, 5, 6, 7}, {[4, 5, 6], [4], [6], [4, 6], [5]});
                        if ismember(f, [2, 4, 5, 6, 7])% was if f==4||f==2 when we have GFP+ data and f=1-4 
                            channels = channel_map(f);
                            idx2 = find(fwhm(idx) > FWMH_threshold & any(peak_data(idx, channels) > intensity_threshold, 2));%idx2=find(fwhm(idx)>FWMH_threshold  & (peak_data(idx,5)>intensity_threshold | peak_data(idx,4)>intensity_threshold | peak_data(idx,6)>intensity_threshold));% was idx2=find(fwhm(idx)>FWMH_threshold  & peak_data(idx,5)>intensity_threshold); % 5 for green fluorescence 4 Irfp (peak_data(idx,5)>intensity_threshold | peak_data(idx,4)>intensity_threshold | peak_data(idx,6)>intensity_threshold)
                            if isempty(idx2)==0
                                peaks=peaks(idx(idx2(:)));
                                locs=locs(idx(idx2(:)));
                                fwhm=fwhm(idx(idx2(:)));
                                peak_data= peak_data(idx(idx2(:)),:);
                                peak_area=peak_area(idx(idx2(:)));
                                fwhm_405=fwhm_405(idx(idx2(:)));
                                fwhm_488=fwhm_488(idx(idx2(:)));
                                fwhm_633=fwhm_633(idx(idx2(:)));
                                fwhm_fl1=fwhm_fl1(idx(idx2(:)));
                                fwhm_fl2=fwhm_fl2(idx(idx2(:)));
                                fwhm_fl3=fwhm_fl3(idx(idx2(:)));
                                peak_area_405=peak_area_405(idx(idx2(:)));
                                peak_area_488=peak_area_488(idx(idx2(:)));
                                peak_area_633=peak_area_633(idx(idx2(:)));
                                peak_area_fl1=peak_area_fl1(idx(idx2(:)));
                                peak_area_fl2=peak_area_fl2(idx(idx2(:)));
                                peak_area_fl3=peak_area_fl3(idx(idx2(:)));
                                File_num=repmat(i,size(peak_data,1),1);

                                all_peaks=[all_peaks,peaks']; %#ok<AGROW>
                                all_locs=[all_locs,locs'];%#ok<AGROW>
                                all_chunks=[all_chunks,(ones(1,length(locs))*ii)];%#ok<AGROW>
                                widths_cum=[widths_cum,fwhm'];%#ok<AGROW>
                                peak_area_cum=[peak_area_cum,peak_area];%#ok<AGROW>
                                peak_values=[peak_values;peak_data];%#ok<AGROW>
                                all_fwhm_405=[all_fwhm_405;fwhm_405];%#ok<AGROW>
                                all_fwhm_488=[all_fwhm_488;fwhm_488];%#ok<AGROW>
                                all_fwhm_633=[all_fwhm_633;fwhm_633]; %#ok<AGROW>
                                all_fwhm_fl1=[all_fwhm_fl1;fwhm_fl1]; %#ok<AGROW>
                                all_fwhm_fl2=[all_fwhm_fl2;fwhm_fl2];%#ok<AGROW>
                                all_fwhm_fl3=[all_fwhm_fl3;fwhm_fl3];%#ok<AGROW>
                                all_peak_area_405=[all_peak_area_405;peak_area_405]; %#ok<AGROW>
                                all_peak_area_488=[all_peak_area_488;peak_area_488]; %#ok<AGROW>
                                all_peak_area_633=[all_peak_area_633;peak_area_633]; %#ok<AGROW>
                                all_peak_area_fl1=[all_peak_area_fl1;peak_area_fl1]; %#ok<AGROW>
                                all_peak_area_fl2=[all_peak_area_fl2;peak_area_fl2];%#ok<AGROW>
                                all_peak_area_fl3=[all_peak_area_fl3;peak_area_fl3];%#ok<AGROW>
                                all_file_num=[all_file_num;File_num];%#ok<AGROW>
                            end
                        else
                            peaks=peaks(idx(:));
                            locs=locs(idx(:));
                            fwhm=fwhm(idx(:));
                            peak_data= peak_data(idx(:),:);
                            peak_area=peak_area(idx(:));
                            fwhm_405=fwhm_405(idx(:));
                            fwhm_488=fwhm_488(idx(:));
                            fwhm_633=fwhm_633(idx(:));
                            fwhm_fl1=fwhm_fl1(idx(:));
                            fwhm_fl2=fwhm_fl2(idx(:));
                            fwhm_fl3=fwhm_fl3(idx(:));
                            peak_area_405=peak_area_405(idx(:));
                            peak_area_488=peak_area_488(idx(:));
                            peak_area_633=peak_area_633(idx(:));
                            peak_area_fl1=peak_area_fl1(idx(:));
                            peak_area_fl2=peak_area_fl2(idx(:));
                            peak_area_fl3=peak_area_fl3(idx(:));
                            File_num=repmat(i,size(peak_data,1),1);

                            all_peaks=[all_peaks,peaks'];%#ok<AGROW>
                            all_locs=[all_locs,locs'];%#ok<AGROW>
                            all_chunks=[all_chunks,(ones(1,length(locs))*ii)];%#ok<AGROW>
                            widths_cum=[widths_cum,fwhm'];%#ok<AGROW>
                            peak_area_cum=[peak_area_cum,peak_area];%#ok<AGROW>
                            peak_values=[peak_values;peak_data];%#ok<AGROW>
                            all_fwhm_405=[all_fwhm_405;fwhm_405]; %#ok<AGROW>
                            all_fwhm_488=[all_fwhm_488;fwhm_488]; %#ok<AGROW>
                            all_fwhm_633=[all_fwhm_633;fwhm_633]; %#ok<AGROW>
                            all_fwhm_fl1=[all_fwhm_fl1;fwhm_fl1]; %#ok<AGROW>
                            all_fwhm_fl2=[all_fwhm_fl2;fwhm_fl2];%#ok<AGROW>
                            all_fwhm_fl3=[all_fwhm_fl3;fwhm_fl3];%#ok<AGROW>
                            all_peak_area_405=[all_peak_area_405;peak_area_405]; %#ok<AGROW>
                            all_peak_area_488=[all_peak_area_488;peak_area_488]; %#ok<AGROW>
                            all_peak_area_633=[all_peak_area_633;peak_area_633]; %#ok<AGROW>
                            all_peak_area_fl1=[all_peak_area_fl1;peak_area_fl1]; %#ok<AGROW>
                            all_peak_area_fl2=[all_peak_area_fl2;peak_area_fl2];%#ok<AGROW>
                            all_peak_area_fl3=[all_peak_area_fl3;peak_area_fl3];%#ok<AGROW>
                            all_file_num=[all_file_num;File_num];%#ok<AGROW>
                        end
                    end
                    clear M_filt peaks locs cumulative peak_data
                    clear fwhm fwhm_405 fwhm_488 fwhm_633 fwhm_fl1 fwhm_fl2 fwhm_fl3
                    clear peak_area peak_area_405 peak_area_488 peak_area_633 peak_area_fl1 peak_area_fl2 peak_area_fl3 File_num
                    toc
                end
            end
            all_peaks_Store=[all_peaks_Store,all_peaks];%#ok<AGROW>
            all_locs_Store=[all_locs_Store,all_locs];%#ok<AGROW>
            all_chunks_Store=[all_chunks_Store,all_chunks];%#ok<AGROW>
            widths_cum_Store=[widths_cum_Store,widths_cum,];%#ok<AGROW>
            peak_area_cum_Store=[peak_area_cum_Store,peak_area_cum];%#ok<AGROW>
            peak_values_Store=[peak_values_Store;peak_values];%#ok<AGROW>
            all_fwhm_405_Store=[all_fwhm_405_Store;all_fwhm_405];%#ok<AGROW>
            all_fwhm_488_Store=[all_fwhm_488_Store;all_fwhm_488];%#ok<AGROW>
            all_fwhm_633_Store=[all_fwhm_633_Store;all_fwhm_633];%#ok<AGROW>
            all_fwhm_fl1_Store=[all_fwhm_fl1_Store;all_fwhm_fl1];%#ok<AGROW>
            all_fwhm_fl2_Store=[all_fwhm_fl2_Store;all_fwhm_fl2];%#ok<AGROW>
            all_fwhm_fl3_Store=[all_fwhm_fl3_Store;all_fwhm_fl3];%#ok<AGROW>
            all_peak_area_405_Store=[all_peak_area_405_Store;all_peak_area_405];%#ok<AGROW>
            all_peak_area_488_Store=[all_peak_area_488_Store;all_peak_area_488];%#ok<AGROW>
            all_peak_area_633_Store=[all_peak_area_633_Store;all_peak_area_633];%#ok<AGROW>
            all_peak_area_fl1_Store=[all_peak_area_fl1_Store;all_peak_area_fl1];%#ok<AGROW>
            all_peak_area_fl2_Store=[all_peak_area_fl2_Store;all_peak_area_fl2];%#ok<AGROW>
            all_peak_area_fl3_Store=[all_peak_area_fl3_Store;all_peak_area_fl3];%#ok<AGROW>
            all_file_num_Store=[all_file_num_Store;all_file_num];%#ok<AGROW>


            peak_values(:,7)=all_chunks; %chunck location
            peak_values(:,8)=all_locs;   %location within chunk
            peak_values(:,9)=all_peaks;  %Cumulative height
            peak_values(:,10)=widths_cum; %Cumulative Width
            peak_values(:,11)=peak_area_cum;
            peak_values(:,12)=all_fwhm_405; peak_values(:,13)=all_peak_area_405;
            peak_values(:,14)=all_fwhm_488; peak_values(:,15)=all_peak_area_488;
            peak_values(:,16)=all_fwhm_633; peak_values(:,17)=all_peak_area_633;
            peak_values(:,18)=all_fwhm_fl1; peak_values(:,19)=all_peak_area_fl1;
            peak_values(:,20)=all_fwhm_fl2; peak_values(:,21)=all_peak_area_fl2;
            peak_values(:,22)=all_fwhm_fl3; peak_values(:,23)=all_peak_area_fl3;
            peak_values(:,24)=all_file_num;
            peak_values(1,:)=[];
            peak_values=sortrows(peak_values,7);
            cd(subdirinfo{i}.folder)
            fileN2=['T',num2str(i),'-',fileN];
            save(fileN2,'peak_values');
            clear all_peaks all_locs all_chunks peak_values widths_cum
            clear all_fwhm_405 all_fwhm_488 all_fwhm_633 all_fwhm_fl1 all_fwhm_fl2 all_fwhm_fl3
            clear peak_area_cum all_peak_area_405 all_peak_area_488 all_peak_area_633 all_peak_area_fl1 all_peak_area_fl2 all_peak_area_fl3 all_file_num

            all_peaks=zeros(1,1);
            all_locs=zeros(1,1);
            widths_cum=zeros(1,1);
            peak_area_cum=zeros(1,1);
            all_chunks=zeros(1,1);
            peak_values=zeros(1,6);
            all_fwhm_405=zeros(1,1); all_fwhm_488=zeros(1,1); all_fwhm_633=zeros(1,1); all_fwhm_fl1=zeros(1,1); all_fwhm_fl2=zeros(1,1); all_fwhm_fl3=zeros(1,1);
            all_peak_area_405=zeros(1,1); all_peak_area_488=zeros(1,1); all_peak_area_633=zeros(1,1); all_peak_area_fl1=zeros(1,1); all_peak_area_fl2=zeros(1,1);all_peak_area_fl3=zeros(1,1);
            all_file_num=zeros(1,1);
        end
    end
    if ~isempty(peak_values_Store)
        peak_values=[]; %#ok<NASGU>
        peak_values=peak_values_Store;
        peak_values(:,7)=all_chunks_Store; %chunck location
        peak_values(:,8)=all_locs_Store;   %location within chunk
        peak_values(:,9)=all_peaks_Store;  %Cumulative height
        peak_values(:,10)=widths_cum_Store; %Cumulative Width
        peak_values(:,11)=peak_area_cum_Store;
        peak_values(:,12)=all_fwhm_405_Store; peak_values(:,13)=all_peak_area_405_Store;
        peak_values(:,14)=all_fwhm_488_Store; peak_values(:,15)=all_peak_area_488_Store;
        peak_values(:,16)=all_fwhm_633_Store; peak_values(:,17)=all_peak_area_633_Store;
        peak_values(:,18)=all_fwhm_fl1_Store; peak_values(:,19)=all_peak_area_fl1_Store;
        peak_values(:,20)=all_fwhm_fl2_Store; peak_values(:,21)=all_peak_area_fl2_Store;
        peak_values(:,22)=all_fwhm_fl3_Store; peak_values(:,23)=all_peak_area_fl3_Store;
        peak_values(:,24)=all_file_num_Store;
        % Sorting peaks
        peak_values=sortrows(peak_values,7);
        delrow=peak_values(:,24)==0;
        peak_values(delrow,:)=[];

        clear cumulative

        % - - - - - - %
        %   S A V E   %
        % - - - - - - %

        cd(mainFolder)

        save(fileN,'peak_values');

        % Save standard deviation and mean of cumulative channel
        cumulative.avg = mean(cum_avg);
        cumulative.stdev = mean(cum_std);
        save(strcat(data_type, '_', 'cumulative_peak_data'), 'cumulative');

    else
        disp('Skipped')
    end
end
if bead_flag==1
    successmessage=BeadSorting6PMT(filepath,file_range);
    disp(successmessage)
end
if bead_flag==1
    successmessage=BeadSorting_GFP6PMT(filepath,file_range);
    disp(successmessage)
end
if bead_flag==1
    successmessage=BeadSorting_iRFP6PMT(filepath,file_range);
    disp(successmessage)
end
%%
successmessage=BeadSorting_GFPwithiRFP6PMT(filepath,file_range,bead_flag);
disp(successmessage)
successmessage=BeadSorting_iRFPwithiGFP6PMT(filepath,file_range,bead_flag);
disp(successmessage)
successmessage='Completed Peak Detection';