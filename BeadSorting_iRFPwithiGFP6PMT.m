function successmessage=BeadSorting_iRFPwithiGFP6PMT(filepath,file_range,bead_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file:BeadSorting_GFP6PMT.m
% ***Description***:
% This function serves automatically search for iRFP Peaks with and without GFP.
% Currently this function is run from within SamplePeakDetection6PMTRF.m function
% Written By: Taras Hanulia (thanul01@tufts.edu) for 6PMT
% Date Written: 11/14/2024

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function details
% Inputs:
%   filepath = Location of all main folder where all subfolders are
%   file_range = Which files do you want to analyze in a specific day
% Outputs:
%   successmessage = a string that indicates the completion of the code
%
% Usage SamplePeakDetection(filepath,file_range)
% Example:
% filepath = 'U:\Nilay\IVFC\Acquired Data\Blood Data\NV_092821_Blood_LNPs';
% file_range= [1:3];
% 
% output=SamplePeakDetection(filepath,file_range)
%% Initialization
clc
%% Find Folders
mainFolder=filepath;
cd(mainFolder)


dirinfo = dir();
dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories

dirinfo(ismember( {dirinfo.name}, {'.', '..'})) = [];  %remove . and ..
[~,c]=natsortfiles({dirinfo.name});
subdirinfo = cell(length(dirinfo));

for K = 1 : length(dirinfo)
    thisdir = dirinfo(c(K)).name;
    if bead_flag == 0
        subdirinfo{K} = dir(fullfile(thisdir, '*iRFPFLR.mat'));
    else
        subdirinfo{K} = dir(fullfile(thisdir, '*iRFPFLRCell.mat'));
    end
end

subdirinfo =  subdirinfo(~cellfun('isempty',subdirinfo));

% if file_range(1)~=1
%     file_range=file_range-(file_range(1)-1);
% end
%% Main Loop: Split iRFP Peaks with and without GFP
for i=file_range
    disp(['Evaluating File # ',num2str(i),' of ',num2str(size(subdirinfo,1))]);
    del=[]; % Matrix of iRFP peaks withGFP
    truepeaks=[]; % Matrix of iRFP peaks without GFP
    if~isempty(subdirinfo{i})
        cd(subdirinfo{i}(1).folder)
        load(subdirinfo{i}(1).name)
        pv1=peak_values;
        if isempty(pv1)
            disp(['Skipping File # ', num2str(i), ' as it contains no peaks.']);
            continue;
        end
        chunks=unique(peak_values(:,7));
        for j=min(chunks):max(chunks)
            disp(['Chunk # ',num2str(j),' of ',num2str(max(chunks))]);
            a = find(pv1(:, 7) == j);
            loc1=pv1(a,8);
            wid1=pv1(a,10);
            if bead_flag == 0
                comparisonFile = [subdirinfo{i}(1).name(1:end-15), 'GFPFLR.mat'];
            else
                comparisonFile = [subdirinfo{i}(1).name(1:end-15), 'GFPFLRCell.mat'];
            end
            load(comparisonFile, 'peak_values')
            pv2=peak_values;
            if isempty(pv2)
                disp(['Skipping File # ', num2str(i), ' as it contains no peaks.']);
                continue;
            end
            a2 = pv2(:, 7) == j;
            loc2=pv2(a2,8);
            wid2=pv2(a2,10);
            for k=1:length(loc1)
                l1=loc1(k);
                match=find(loc2>=l1-wid1(k).*1.5 & loc2<=l1+wid1(k).*1.5,1);
                if ~isempty(match)
                    del=[del;pv1(a(k),:)]; %#ok<AGROW> 
                else
                    truepeaks=[truepeaks;pv1(a(k),:)]; %#ok<AGROW> 
                end
            end
        end
        peak_values=truepeaks;
        save([subdirinfo{i}(1).name(1:end-15),'iRFPFLR_Single.mat'],'peak_values')
        peak_values=del;
        save([subdirinfo{i}(1).name(1:end-15),'iRFPFLR_with_GFP.mat'],'peak_values')
    end 
end

% %% Save all sorted files
% pv1=[];
% pv2=[];
% for i=file_range
%     disp(['Evaluating File # ',num2str(i),' of ',num2str(size(subdirinfo,1))]);
%     if~isempty(subdirinfo{i}(1))
%         cd(subdirinfo{i}(1).folder)
%         load([subdirinfo{i}(1).name(1:end-15),'iRFPFLR_with_GFP.mat'],'peak_values')
%         pv1=[pv1;peak_values]; %#ok<AGROW> 
%         load([subdirinfo{i}(1).name(1:end-15),'iRFPFLR_Single.mat'],'peak_values')
%         pv2=[pv2;peak_values]; %#ok<AGROW> 
%     end
% end
% cd(mainFolder)
%  peak_values=pv1;
%  save([subdirinfo{i}(1).name(4:end-15),'iRFPFLR_with_GFP.mat'],'peak_values')
%  peak_values=pv2;
%  save([subdirinfo{i}(1).name(4:end-15),'iRFPFLR_Single.mat'],'peak_values')
% successmessage='Completed Sorting iRFP Peaks With and Without GFP FLR';
%%
% pv1 = [];
% pv2 = [];
% for i = file_range
%     disp(['Evaluating File # ', num2str(i), ' of ', num2str(size(subdirinfo,1))]);
%     
%     if ~isempty(subdirinfo{i}(1))
%         cd(subdirinfo{i}(1).folder);
% 
%         % Define file names
%         file_with_GFP = [subdirinfo{i}(1).name(1:end-15), 'iRFPFLR_with_GFP.mat'];
%         file_single = [subdirinfo{i}(1).name(1:end-15), 'iRFPFLR_Single.mat'];
% 
%         % Check if 'iRFPFLR_with_GFP.mat' exists before loading
%         if exist(file_with_GFP, 'file')
%             load(file_with_GFP, 'peak_values');
%             pv1 = [pv1; peak_values]; %#ok<AGROW> 
%         else
%             disp(['Skipping: ', file_with_GFP, ' not found.']);
%         end
% 
%         % Check if 'iRFPFLR_Single.mat' exists before loading
%         if exist(file_single, 'file')
%             load(file_single, 'peak_values');
%             pv2 = [pv2; peak_values]; %#ok<AGROW> 
%         else
%             disp(['Skipping: ', file_single, ' not found.']);
%         end
%     end
% end
% 
% cd(mainFolder);
% 
% % Save only if data exists
% if ~isempty(pv1)
%     peak_values = pv1;
%     save([subdirinfo{i}(1).name(4:end-15), 'iRFPFLR_with_GFP.mat'], 'peak_values');
% else
%     disp('No data found for iRFPFLR_with_GFP.mat, skipping save.');
% end
% 
% if ~isempty(pv2)
%     peak_values = pv2;
%     save([subdirinfo{i}(1).name(4:end-15), 'iRFPFLR_Single.mat'], 'peak_values');
% else
%     disp('No data found for iRFPFLR_Single.mat, skipping save.');
% end
% 
% successmessage = 'Completed Sorting iRFP Peaks With and Without GFP FLR';
%%
pv1 = [];
pv2 = [];

for i = file_range
    disp(['Evaluating File # ', num2str(i), ' of ', num2str(size(subdirinfo, 1))]);

    if ~isempty(subdirinfo{i}(1))
        cd(subdirinfo{i}(1).folder);

        % Ensure filename is long enough before extracting substrings
        base_name = subdirinfo{i}(1).name;
        if length(base_name) > 18  % Adjusted length check for safety
            save_filename1 = [base_name(4:end-15), 'iRFPFLR_with_GFP.mat'];
            save_filename2 = [base_name(4:end-15), 'iRFPFLR_Single.mat'];
        else
            warning(['Filename too short: ', base_name, ' Skipping file.']);
            continue;
        end

        % Define file names for loading
        file_with_GFP = [base_name(1:end-15), 'iRFPFLR_with_GFP.mat'];
        file_single = [base_name(1:end-15), 'iRFPFLR_Single.mat'];

        % Load data if the files exist
        if exist(file_with_GFP, 'file')
            load(file_with_GFP, 'peak_values');
            pv1 = [pv1; peak_values]; %#ok<AGROW>
        else
            disp(['Skipping: ', file_with_GFP, ' not found.']);
        end

        if exist(file_single, 'file')
            load(file_single, 'peak_values');
            pv2 = [pv2; peak_values]; %#ok<AGROW>
        else
            disp(['Skipping: ', file_single, ' not found.']);
        end
    end
end

cd(mainFolder);

% Ensure subdirinfo{i} is within bounds before using it
if ~isempty(subdirinfo) && i <= length(subdirinfo)
    base_name = subdirinfo{i}(1).name;
    if length(base_name) > 18  % Ensure the filename is long enough
        save_filename1 = [base_name(4:end-15), 'iRFPFLR_with_GFP.mat'];
        save_filename2 = [base_name(4:end-15), 'iRFPFLR_Single.mat'];
    else
        warning('Filename too short, skipping save.');
        return;
    end
else
    warning('subdirinfo index out of bounds, skipping save.');
    return;
end

% Save only if data exists
if ~isempty(pv1)
    peak_values = pv1;
    save(fullfile(mainFolder, save_filename1), 'peak_values');
    disp(['Saved: ', fullfile(mainFolder, save_filename1)]);
else
    disp('No data found for iRFPFLR_with_GFP.mat, skipping save.');
end

if ~isempty(pv2)
    peak_values = pv2;
    save(fullfile(mainFolder, save_filename2), 'peak_values');
    disp(['Saved: ', fullfile(mainFolder, save_filename2)]);
else
    disp('No data found for iRFPFLR_Single.mat, skipping save.');
end

successmessage = 'Completed Sorting iRFP Peaks With and Without GFP FLR';