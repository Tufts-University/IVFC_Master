function successmessage = BeadSorting_iRFP6PMT(filepath, file_range)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file: BeadSorting_iRFP6PMT.m
%
% *** Description ***
% This function automatically searches for Bead Peaks and sorts them
% from the Cell peaks in Blood samples where Beads were also spiked.
% It is executed from within SamplePeakDetection6PMTRF.m when a flag
% indicates that Beads were added to the blood sample.
%
% Written By: Nilay Vora (nvora01@tufts.edu) for 5PMT
% Date Written: 01/13/2022
%
% Modifying Author: Nilay Vora
% Date Modified: 01/19/2022
% Latest Revision: Manually modified indexing due to an error in lines 49-52.
%
% Modifying Author: Taras Hanulia
% Date Modified: 06/30/2022
% Latest Revision: 6PMT modification
%
% Modifying Author: Taras Hanulia
% Date Modified: 11/14/2022
% Latest Revision: 6PMT modification and iRFP
% Beads are in channel E2C.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Function details
% Inputs:
%   filepath   - Location of the main folder where all subfolders are stored.
%   file_range - Specifies which files to analyze for a given day.
%
% Outputs:
%   successmessage - A string indicating the completion of the process.
%
% Usage:
%   BeadSorting_iRFP6PMT(filepath, file_range)
%
% Example:
%   filepath = 'U:\Nilay\IVFC\Acquired Data\Blood Data\NV_092821_Blood_LNPs';
%   file_range = 1:3;
%   output = BeadSorting_iRFP6PMT(filepath, file_range);

%% Initialization
clc;

%% Find Folders
mainFolder = filepath;
cd(mainFolder);

dirinfo = dir();
dirinfo(~[dirinfo.isdir]) = [];  % Remove non-directory entries

dirinfo(ismember({dirinfo.name}, {'.', '..'})) = [];  % Remove '.' and '..'
[~, c] = natsortfiles({dirinfo.name});
subdirinfo = cell(length(dirinfo), 1);

for K = 1:length(dirinfo)
    thisdir = dirinfo(c(K)).name;
    subdirinfo{K} = dir(fullfile(thisdir, '*E2CPFLR.mat'));
end

subdirinfo = subdirinfo(~cellfun('isempty', subdirinfo));

%% Main Loop: Split the Cells from Beads
for i = file_range
    disp(['Evaluating File # ', num2str(i), ' of ', num2str(size(subdirinfo, 1))]);
    del = [];       % Matrix of bead peaks
    truepeaks = []; % Matrix of cell peaks
    
    if ~isempty(subdirinfo{i})
        cd(subdirinfo{i}(1).folder);
        load(subdirinfo{i}(1).name, 'peak_values');
        pv1 = peak_values;
        chunks = unique(pv1(:, 7));
        
        for j = min(chunks):max(chunks)
            disp(['Chunk # ', num2str(j), ' of ', num2str(max(chunks))]);
            a = pv1(:, 7) == j;
            loc = pv1(a, 8);
            
            load([subdirinfo{i}(1).name(1:end-11), 'iRFPFLRAll.mat'], 'peak_values');
            pv2 = peak_values;
            a2 = find(pv2(:, 7) == j);
            loc2 = pv2(a2, 8);
            wid2 = pv2(a2, 10);
            
            for k = 1:length(loc2)
                l2 = loc2(k);
                match = find(loc >= l2 - wid2(k) * 1.5 & loc <= l2 + wid2(k) * 1.5, 1);
                
                if ~isempty(match)
                    del = [del; pv2(a2(k), :)]; %#ok<AGROW>
                else
                    truepeaks = [truepeaks; pv2(a2(k), :)]; %#ok<AGROW>
                end
            end
        end
        
        peak_values = truepeaks;
        save([subdirinfo{i}(1).name(1:end-11), 'iRFPFLRCell.mat'], 'peak_values');
        peak_values = del;
        save([subdirinfo{i}(1).name(1:end-11), 'iRFPFLRBeads.mat'], 'peak_values');
    end 
end

%% Save all sorted files
pv1 = [];
pv2 = [];

for i = file_range
    disp(['Evaluating File # ', num2str(i), ' of ', num2str(size(subdirinfo, 1))]);
    
    if ~isempty(subdirinfo{i})
        cd(subdirinfo{i}(1).folder);
        
        load([subdirinfo{i}(1).name(1:end-11), 'iRFPFLRBeads.mat'], 'peak_values');
        pv1 = [pv1; peak_values]; %#ok<AGROW>
        
        load([subdirinfo{i}(1).name(1:end-11), 'iRFPFLRCell.mat'], 'peak_values');
        pv2 = [pv2; peak_values]; %#ok<AGROW>
    end
end

cd(mainFolder);
peak_values = pv1;
save([subdirinfo{i}(1).name(5:end-11), 'iRFPFLRBeads.mat'], 'peak_values');

peak_values = pv2;
save([subdirinfo{i}(1).name(5:end-11), 'iRFPFLRCell.mat'], 'peak_values');

successmessage = 'Completed Separating iRFP Bead Peaks';