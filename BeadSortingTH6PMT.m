function successmessage=BeadSortingTH6PMT(filepath,file_range)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file:BeadSorting.m
% ***Description***:
% This function serves automatically search for Bead Peaks and sort them
% from the Cell peaks in Blood samples where Beads were also spiked.
% Currently this function is run from within SamplePeakDetection.m function
% when a flag indicated Beads were added to blood sample.
% Written By: Nilay Vora (nvora01@tufts.edu)
% Date Written: 01/13/2022
% Modifying Author:Nilay Vora
% Date Modified: 01/19/2022
% Latest Revision: Manually modified indexing due to error in lines 49-52,
% not sure if this is the right correction.
% Date Modified: 08/21/2024
% Modifying Author:Taras Hanulia
% modification of peak width to peak start and end
%Date Modified: 04/09/2025
% Modifying Author:Taras Hanulia
% modification for *NTH* data
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
    subdirinfo{K} = dir(fullfile(thisdir, '*NTH*RFLR.mat'));
end

subdirinfo =  subdirinfo(~cellfun('isempty',subdirinfo));

if file_range(1)~=1
    file_range=file_range-(file_range(1)-1);
end
%% Main Loop: Split the Cells from Beads
for i=file_range
    disp(['Evaluating File # ',num2str(i),' of ',num2str(size(subdirinfo,1))]);
    del=[]; % Matrix of bead peaks
    truepeaks=[]; % Matrix of cell peaks
    if~isempty(subdirinfo{i})
        cd(subdirinfo{i}(1).folder)
        load(subdirinfo{i}(1).name)
        pv1=peak_values;
        chunks=unique(peak_values(:,7));
        for j=min(chunks):max(chunks)
            disp(['Chunk # ',num2str(j),' of ',num2str(max(chunks))]);
            a=pv1(:,7)==j;
            loc=pv1(a,8);
            load([subdirinfo{i}(1).name(1:end-8),'NoScatAll.mat'],'peak_values')
            pv2=peak_values;
            a2=find(pv2(:,7)==j);
            loc2=pv2(a2,8);
%             start2 = pv2(a2, 9);
%             end2 = pv2(a2, 10);
            wid2=pv2(a2,10);
            for k=1:length(loc2)
                l2=loc2(k);
                match=find(loc>=l2-wid2(k).*1.5 & loc<=l2+wid2(k).*1.5,1);
%                 match = find(loc >= start2(k) & loc <= end2(k), 1);
                if ~isempty(match)
                    del=[del;pv2(a2(k),:)]; %#ok<AGROW> 
                else
                    truepeaks=[truepeaks;pv2(a2(k),:)]; %#ok<AGROW> 
                end
            end
        end
        peak_values=truepeaks;
        save([subdirinfo{i}(1).name(1:end-8),'NoScatCell.mat'],'peak_values')
        peak_values=del;
        save([subdirinfo{i}(1).name(1:end-8),'NoScatBeads.mat'],'peak_values')
    end 
end

% %% Save all sorted files
% pv1=[];
% pv2=[];
% for i=file_range
%     disp(['Evaluating File # ',num2str(i),' of ',num2str(size(subdirinfo,1))]);
%     if~isempty(subdirinfo{i}(1))
%         cd(subdirinfo{i}(1).folder)
%         load([subdirinfo{i}(1).name(1:end-8),'NoScatBeads.mat'],'peak_values')
%         pv1=[pv1;peak_values]; %#ok<AGROW> 
%         load([subdirinfo{i}(1).name(1:end-8),'NoScatCell.mat'],'peak_values')
%         pv2=[pv2;peak_values]; %#ok<AGROW> 
%     end
% end
% cd(mainFolder)
%  peak_values=pv1;
%  save([subdirinfo{i}(1).name(4:end-8),'NoScatBeads.mat'],'peak_values')
%  peak_values=pv2;
%  save([subdirinfo{i}(1).name(4:end-8),'NoScatCell.mat'],'peak_values')
% 
%  
% successmessage='Completed Seperating Bead Peaks';
%% Save all sorted files (skip files with no data)
pv1 = [];
pv2 = [];

for i = file_range
    disp(['Evaluating File # ', num2str(i), ' of ', num2str(size(subdirinfo, 1))]);

    if ~isempty(subdirinfo{i}(1))
        cd(subdirinfo{i}(1).folder);
        base_name = subdirinfo{i}(1).name(1:end-8);

        % Try to load Beads data
        beads_file = [base_name, 'NoScatBeads.mat'];
        if isfile(beads_file)
            temp = load(beads_file, 'peak_values');
            if isfield(temp, 'peak_values') && ~isempty(temp.peak_values)
                pv1 = [pv1; temp.peak_values]; %#ok<AGROW>
            end
        end

        % Try to load Cell data
        cell_file = [base_name, 'NoScatCell.mat'];
        if isfile(cell_file)
            temp = load(cell_file, 'peak_values');
            if isfield(temp, 'peak_values') && ~isempty(temp.peak_values)
                pv2 = [pv2; temp.peak_values]; %#ok<AGROW>
            end
        end
    end
end


cd(mainFolder)

% Get full original name
original_name = subdirinfo{i}(1).name;

% Remove prefix like 'T19-' or 'T7-'
name_no_prefix = regexprep(original_name, '^T\d{1,2}-', '');

% Remove suffix like '_RFLR.mat'
base_name = regexprep(name_no_prefix, '_[^_]+\.mat$', '');

% Save combined Beads
if ~isempty(pv1)
    peak_values = pv1;
    save([base_name, '_NoScatBeads.mat'], 'peak_values');
end

% Save combined Cells
if ~isempty(pv2)
    peak_values = pv2;
    save([base_name, '_NoScatCell.mat'], 'peak_values');
end
successmessage = 'Completed Separating Bead Peaks';
