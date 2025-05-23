function successmessage=BeadSorting_iRFPwithiGFP_withScat6PMT(filepath,file_range)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file:BeadSorting_GFP6PMT.m
% ***Description***:
% This function serves automatically search for iRFP Peaks with GFP for with Scat and without .
% Currently this function is run from within SamplePeakDetection6PMTRF.m function
% Written By: Taras Hanulia (thanul01@tufts.edu) for 6PMT
% Date Written: 01/14/2025

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
    
    subdirinfo{K} = dir(fullfile(thisdir, '*iRFPFLR_with_GFP.mat'));
    
end

subdirinfo =  subdirinfo(~cellfun('isempty',subdirinfo));

% if file_range(1)~=1
%     file_range=file_range-(file_range(1)-1);
% end
%% Main Loop: Split iRFP Peaks with and without GFP
for i=file_range
    if i > length(subdirinfo) || isempty(subdirinfo{i}) % new021825
        warning(['Skipping file #', num2str(i), ' as it does not exist.']);% new021825
        continue;% new021825
    end% new021825
    disp(['Evaluating File # ',num2str(i),' of ',num2str(size(subdirinfo,1))]);
    del=[]; % Matrix of iRFP peaks withGFP
    truepeaks=[]; % Matrix of iRFP peaks without GFP
    if~isempty(subdirinfo{i})
        cd(subdirinfo{i}(1).folder)
        load(subdirinfo{i}(1).name)
        if ~exist('peak_values', 'var') || isempty(peak_values)
            warning(['File ', subdirinfo{i}(1).name, ' does not contain valid peak_values. Skipping.']);
            continue;
        end
        
        if size(peak_values, 2) < 7
            warning(['File ', subdirinfo{i}(1).name, ' has invalid peak_values structure. Skipping.']);
            continue;
        end
        pv1=peak_values;
        chunks=unique(peak_values(:,7));
        for j=min(chunks):max(chunks)
            disp(['Chunk # ',num2str(j),' of ',num2str(max(chunks))]);
            a = find(pv1(:, 7) == j);
            loc1=pv1(a,8);
            wid1=pv1(a,10);
            comparisonFile = [subdirinfo{i}(1).name(1:end-20), 'NoFLR.mat'];
            load(comparisonFile, 'peak_values')
            pv2=peak_values;
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
        save([subdirinfo{i}(1).name(1:end-20),'iRFPFLR_with_GFP_noScat.mat'],'peak_values')
        peak_values=del;
        save([subdirinfo{i}(1).name(1:end-20),'iRFPFLR_with_GFP_withScat.mat'],'peak_values')
    end 
end

% %% Save all sorted files
% pv1 = [];
% pv2 = [];
% for i = file_range
%     if i > length(subdirinfo) || isempty(subdirinfo{i}) % new021825
%         warning(['Skipping file #', num2str(i), ' as it does not exist.']);% new021825
%         continue;% new021825
%     end% new021825
%     disp(['Evaluating File # ', num2str(i), ' of ', num2str(size(subdirinfo, 1))]);
%     if ~isempty(subdirinfo{i}(1))
%         cd(subdirinfo{i}(1).folder);
%         
%         % Construct file paths
%         file_with_scat = [subdirinfo{i}(1).name(1:end-20), 'iRFPFLR_with_GFP_withScat.mat'];
%         file_no_scat = [subdirinfo{i}(1).name(1:end-20), 'iRFPFLR_with_GFP_noScat.mat'];
%         
%         % Load file_with_scat
%         if isfile(file_with_scat)
%             load(file_with_scat, 'peak_values');
%             if exist('peak_values', 'var') && ~isempty(peak_values) && size(peak_values, 2) >= 7
%                 pv1 = [pv1; peak_values]; %#ok<AGROW>
%             else
%                 warning(['File ', file_with_scat, ' does not contain valid peak_values. Skipping.']);
%             end
%         else
%             warning(['File ', file_with_scat, ' not found. Skipping.']);
%         end
%         
%         % Load file_no_scat
%         if isfile(file_no_scat)
%             load(file_no_scat, 'peak_values');
%             if exist('peak_values', 'var') && ~isempty(peak_values) && size(peak_values, 2) >= 7
%                 pv2 = [pv2; peak_values]; %#ok<AGROW>
%             else
%                 warning(['File ', file_no_scat, ' does not contain valid peak_values. Skipping.']);
%             end
%         else
%             warning(['File ', file_no_scat, ' not found. Skipping.']);
%         end
%     end
% end
% cd(mainFolder)
% % peak_values=pv1;
% % save([subdirinfo{i}(1).name(4:end-20),'iRFPFLR_with_GFP_withScat.mat'],'peak_values')
% % peak_values=pv2;
% % save([subdirinfo{i}(1).name(4:end-20),'iRFPFLR_with_GFP_noScat.mat'],'peak_values')
% if ~isempty(pv1)
%     save([subdirinfo{i}(1).name(4:end-20),'iRFPFLR_with_GFP_withScat.mat'], 'pv1');
% end
% if ~isempty(pv2)
%     save([subdirinfo{i}(1).name(4:end-20),'iRFPFLR_with_GFP_noScat.mat'], 'pv2');
% end
% successmessage='Completed Sorting iRFP Peaks With GFP FLR With and Without Scattering';
%%
mainFolder = pwd;

% Initialize the variables
pv1 = [];
pv2 = [];

% Process each file in the range
for i = file_range
    % Ensure that i is within bounds and subdirinfo{i} is not empty
    if i > length(subdirinfo) || isempty(subdirinfo{i})
        warning(['Skipping file #', num2str(i), ' as it does not exist.']);
        continue;
    end
    
    disp(['Evaluating File # ', num2str(i), ' of ', num2str(length(subdirinfo))]);
    
    % Check for files in subdirinfo and process
    if ~isempty(subdirinfo{i}) && ~isempty(subdirinfo{i}(1).name)
        cd(subdirinfo{i}(1).folder);  % Change directory to the file location
        
        % Set last folder name for saving later
        lastFolderName = subdirinfo{i}(1).name(1:end-20);  % Extract name prefix to use when saving
        
        % Define file names
        file_with_scat = [lastFolderName, 'iRFPFLR_with_GFP_withScat.mat'];
        file_no_scat = [lastFolderName, 'iRFPFLR_with_GFP_noScat.mat'];
        
        % Load and process the file with scattering
        if isfile(file_with_scat)
            load(file_with_scat, 'peak_values');
            if exist('peak_values', 'var') && ~isempty(peak_values) && size(peak_values, 2) >= 7
                pv1 = [pv1; peak_values];  %#ok<AGROW>
            else
                warning(['File ', file_with_scat, ' does not contain valid peak_values. Skipping.']);
            end
        else
            warning(['File ', file_with_scat, ' not found. Skipping.']);
        end
        
        % Load and process the file without scattering
        if isfile(file_no_scat)
            load(file_no_scat, 'peak_values');
            if exist('peak_values', 'var') && ~isempty(peak_values) && size(peak_values, 2) >= 7
                pv2 = [pv2; peak_values];  %#ok<AGROW>
            else
                warning(['File ', file_no_scat, ' does not contain valid peak_values. Skipping.']);
            end
        else
            warning(['File ', file_no_scat, ' not found. Skipping.']);
        end
    else
        warning(['Subdirectory info for file #', num2str(i), ' is empty. Skipping.']);
    end
end

% Extract the parent folder (one level up from the last subfolder)
parentFolder = fileparts(mainFolder);  % This gives the parent folder of mainFolder

disp(['Parent Folder: ', parentFolder]);  % Debug message to check the parent path

% Ensure saving is done in the parent folder
% Check if pv1 (with scattering) is not empty and save as peak_values
if ~isempty(pv1)
    peak_values = pv1;  % Assign pv1 to peak_values
    fullPath = fullfile(parentFolder, [lastFolderName, 'iRFPFLR_with_GFP_withScat.mat']);
    disp(['Saving to: ', fullPath]);  % Debug message to check the save location
    save(fullPath, 'peak_values');
end

% Check if pv2 (without scattering) is not empty and save as peak_values
if ~isempty(pv2)
    peak_values = pv2;  % Assign pv2 to peak_values
    fullPath = fullfile(parentFolder, [lastFolderName, 'iRFPFLR_with_GFP_noScat.mat']);
    disp(['Saving to: ', fullPath]);  % Debug message to check the save location
    save(fullPath, 'peak_values');
end

% Final success message
successmessage = 'Completed Sorting iRFP Peaks With GFP FLR With and Without Scattering';

% function successmessage = BeadSorting_iRFPwithiGFP_withScat6PMT(filepath, file_range)
% % Function to sort iRFP Peaks with GFP, with and without Scattering
% % Written By: Taras Hanulia (thanul01@tufts.edu)
% % Date Written: 01/14/2025
% 
% clc;
% mainFolder = filepath;
% cd(mainFolder);
% 
% dirinfo = dir();
% dirinfo(~[dirinfo.isdir]) = [];  % Remove non-directories
% 
% dirinfo(ismember({dirinfo.name}, {'.', '..'})) = [];  % Remove . and ..
% [~, c] = natsortfiles({dirinfo.name});
% subdirinfo = cell(length(dirinfo), 1);
% 
% for K = 1:length(dirinfo)
%     thisdir = dirinfo(c(K)).name;
%     subdirinfo{K} = dir(fullfile(thisdir, '*iRFPFLR_with_GFP.mat'));
% end
% 
% subdirinfo = subdirinfo(~cellfun('isempty', subdirinfo));
% 
% for i = file_range
%     if i > length(subdirinfo) || isempty(subdirinfo{i})
%         warning(['Skipping file index ', num2str(i), ' due to out-of-bounds or empty entry.']);
%         continue;
%     end
%     disp(['Evaluating File # ', num2str(i), ' of ', num2str(size(subdirinfo, 1))]);
%     del = [];
%     truepeaks = [];
%     
%     if ~isempty(subdirinfo{i})
%         cd(subdirinfo{i}(1).folder);
%         load(subdirinfo{i}(1).name, 'peak_values');
%         
%         if ~exist('peak_values', 'var') || isempty(peak_values) || size(peak_values, 2) < 7
%             warning(['Skipping file: ', subdirinfo{i}(1).name]);
%             continue;
%         end
%         
%         pv1 = peak_values;
%         chunks = unique(pv1(:, 7));
%         
%         for j = min(chunks):max(chunks)
%             disp(['Chunk # ', num2str(j), ' of ', num2str(max(chunks))]);
%             a = find(pv1(:, 7) == j);
%             loc1 = pv1(a, 8);
%             wid1 = pv1(a, 10);
%             
%             comparisonFile = [subdirinfo{i}(1).name(1:end-20), 'NoFLR.mat'];
%             load(comparisonFile, 'peak_values');
%             
%             if ~exist('peak_values', 'var') || isempty(peak_values)
%                 warning(['Skipping comparison file: ', comparisonFile]);
%                 continue;
%             end
%             
%             pv2 = peak_values;
%             loc2 = pv2(pv2(:, 7) == j, 8);
%             wid2 = pv2(pv2(:, 7) == j, 10);
%             
%             for k = 1:length(loc1)
%                 l1 = loc1(k);
%                 match = find(loc2 >= l1 - wid1(k) * 1.5 & loc2 <= l1 + wid1(k) * 1.5, 1);
%                 if ~isempty(match)
%                     del = [del; pv1(a(k), :)];
%                 else
%                     truepeaks = [truepeaks; pv1(a(k), :)];
%                 end
%             end
%         end
%         
%         save([subdirinfo{i}(1).name(1:end-20), 'iRFPFLR_with_GFP_noScat.mat'], 'peak_values', '-v7.3');
%         save([subdirinfo{i}(1).name(1:end-20), 'iRFPFLR_with_GFP_withScat.mat'], 'peak_values', '-v7.3');
%     end
% end
% %%
% pv1 = [];
% pv2 = [];
% for i = file_range
%     if i > length(subdirinfo)
%         warning(['Skipping file index ', num2str(i), ' as it exceeds available subdirectories.']);
%         continue;
%     end
%     if isempty(subdirinfo{i}), continue;
%     end
%     cd(subdirinfo{i}(1).folder);
%     
%     file_with_scat = [subdirinfo{i}(1).name(1:end-20), 'iRFPFLR_with_GFP_withScat.mat'];
%     file_no_scat = [subdirinfo{i}(1).name(1:end-20), 'iRFPFLR_with_GFP_noScat.mat'];
%     
%     if isfile(file_with_scat)
%         load(file_with_scat, 'peak_values');
%         if exist('peak_values', 'var') && ~isempty(peak_values)
%             pv1 = [pv1; peak_values];
%         end
%     end
%     
%     if isfile(file_no_scat)
%         load(file_no_scat, 'peak_values');
%         if exist('peak_values', 'var') && ~isempty(peak_values)
%             pv2 = [pv2; peak_values];
%         end
%     end
% end
% 
% cd(mainFolder);
% save([subdirinfo{1}(1).name(4:end-20), 'iRFPFLR_with_GFP_withScat.mat'], 'pv1', '-v7.3');
% save([subdirinfo{1}(1).name(4:end-20), 'iRFPFLR_with_GFP_noScat.mat'], 'pv2', '-v7.3');
% 
% successmessage = 'Completed Sorting iRFP Peaks With GFP FLR With and Without Scattering';
% end
