function successmessage=BeadSorting_GFPwithiRFP6PMT(filepath,file_range,bead_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file:BeadSorting_GFP6PMT.m
% ***Description***:
% This function serves automatically search for GFP Peaks with and without iRFP.
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
        subdirinfo{K} = dir(fullfile(thisdir, '*GFPFLR.mat'));
    else
        subdirinfo{K} = dir(fullfile(thisdir, '*GFPFLRCell.mat'));
    end
end

subdirinfo =  subdirinfo(~cellfun('isempty',subdirinfo));

% if file_range(1)~=1
%     file_range=file_range-(file_range(1)-1);
% end
%% Main Loop: Split the GFP Peaks with and without iRFP
for i=file_range
    disp(['Evaluating File # ',num2str(i),' of ',num2str(size(subdirinfo,1))]);
    del=[]; % Matrix the GFP Peaks with iRFP
    truepeaks=[]; % Matrix of the GFP Peaks without iRFP
    if~isempty(subdirinfo{i})
        cd(subdirinfo{i}(1).folder)
        load(subdirinfo{i}(1).name)
        pv1=peak_values;
        if isempty(pv1)
            disp(['Skipping chunk ', num2str(i), ' as no peaks are found in the file.']);
            continue;
        end
        chunks=unique(peak_values(:,7));
        for j=min(chunks):max(chunks)
            disp(['Chunk # ',num2str(j),' of ',num2str(max(chunks))]);
            a = find(pv1(:, 7) == j);
            loc1=pv1(a,8);
            wid1=pv1(a,10);
            if bead_flag == 0
                comparisonFile = [subdirinfo{i}(1).name(1:end-14), 'iRFPFLR.mat'];
            else
                comparisonFile = [subdirinfo{i}(1).name(1:end-14), 'iRFPFLRCell.mat'];
            end
            load(comparisonFile, 'peak_values')
            pv2=peak_values;
            if isempty(pv2)
                disp(['Skipping chunk ', num2str(j), ' as no peaks are found in the file.']);
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
        save([subdirinfo{i}(1).name(1:end-14),'GFPFLR_Single.mat'],'peak_values')
        peak_values=del;
        save([subdirinfo{i}(1).name(1:end-14),'GFPFLR_with_iRFP.mat'],'peak_values')
    end 
end

% %% Save all sorted files
% pv1=[];
% pv2=[];
% for i=file_range
%     disp(['Evaluating File # ',num2str(i),' of ',num2str(size(subdirinfo,1))]);
%     if~isempty(subdirinfo{i}(1))
%         cd(subdirinfo{i}(1).folder)
%         load([subdirinfo{i}(1).name(1:end-14),'GFPFLR_with_iRFP.mat'],'peak_values')
%         pv1=[pv1;peak_values]; %#ok<AGROW> 
%         load([subdirinfo{i}(1).name(1:end-14),'GFPFLR_Single.mat'],'peak_values')
%         pv2=[pv2;peak_values]; %#ok<AGROW> 
%     end
% end
% cd(mainFolder)
%  peak_values=pv1;
%  save([subdirinfo{i}(1).name(4:end-14),'GFPFLR_with_iRFP.mat'],'peak_values')
%  peak_values=pv2;
%  save([subdirinfo{i}(1).name(4:end-14),'GFPFLR_Single.mat'],'peak_values')
% successmessage='Completed Sorting GFP Peaks With and Without iRFP FLR';
%%
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
%         file_with_GFP = [subdirinfo{i}(1).name(1:end-14), 'GFPFLR_with_iRFP.mat'];
%         file_single = [subdirinfo{i}(1).name(1:end-14), 'GFPFLR_Single.mat'];
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
%     save([subdirinfo{i}(1).name(4:end-14), 'GFPFLR_with_iRFP.mat'], 'peak_values');
% else
%     disp('No data found for GFPFLR_with_iRFP.mat, skipping save.');
% end
% 
% if ~isempty(pv2)
%     peak_values = pv2;
%     save([subdirinfo{i}(1).name(4:end-14), 'GFPFLR_Single.mat'], 'peak_values');
% else
%     disp('No data found for GFPFLR_Single.mat, skipping save.');
% end
% 
% successmessage = 'Completed Sorting GFP Peaks With and Without iRFP FLR';
%%
pv1 = [];
pv2 = [];

for i = file_range
    disp(['Evaluating File # ', num2str(i), ' of ', num2str(size(subdirinfo,1))]);
    
    if ~isempty(subdirinfo{i}(1))
        cd(subdirinfo{i}(1).folder);

        % Extract base filename dynamically
        base_filename = subdirinfo{i}(1).name;
        
        % Find where 'cell' starts to handle varying lengths
        cell_idx = strfind(base_filename, 'cell');
        if ~isempty(cell_idx)
            base_name = base_filename(1:cell_idx(1)-1); % Extract up to "cellX"
        else
            base_name = base_filename(1:end-14); % Fallback if 'cell' is missing
        end

        % Define file names
        file_with_GFP = [base_name, 'GFPFLR_with_iRFP.mat'];
        file_single = [base_name, 'GFPFLR_Single.mat'];

        % Load 'GFPFLR_with_iRFP.mat' if it exists
        if exist(file_with_GFP, 'file')
            load(file_with_GFP, 'peak_values');
            pv1 = [pv1; peak_values]; %#ok<AGROW>
        else
            disp(['Skipping: ', file_with_GFP, ' not found.']);
        end

        % Load 'GFPFLR_Single.mat' if it exists
        if exist(file_single, 'file')
            load(file_single, 'peak_values');
            pv2 = [pv2; peak_values]; %#ok<AGROW>
        else
            disp(['Skipping: ', file_single, ' not found.']);
        end
    end
end

cd(mainFolder);

% Save only if data exists, using fullfile for safety
if ~isempty(pv1)
    peak_values = pv1;
    save_filename = fullfile(mainFolder, [base_name, 'GFPFLR_with_iRFP.mat']);
    disp(['Saving: ', save_filename]);
    save(save_filename, 'peak_values');
else
    disp('No data found for GFPFLR_with_iRFP.mat, skipping save.');
end

if ~isempty(pv2)
    peak_values = pv2;
    save_filename = fullfile(mainFolder, [base_name, 'GFPFLR_Single.mat']);
    disp(['Saving: ', save_filename]);
    save(save_filename, 'peak_values');
else
    disp('No data found for GFPFLR_Single.mat, skipping save.');
end

successmessage = 'Completed Sorting GFP Peaks With and Without iRFP FLR';
