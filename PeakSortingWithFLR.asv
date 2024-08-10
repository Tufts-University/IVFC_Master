function successmessage = PeakSortingWithFLR(filepath, file_range, bead_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file: PeakSortingWithFLR.m
% ***Description***:
% This function searches for peaks in files ending with NoFLR.mat, compares
% them to peaks in files ending with NoScatAll.mat or NoScatCell.mat, and 
% then sorts the peaks into those with FLR and without FLR.
% Inputs:
%   filepath = Location of all main folder where all subfolders are
%   file_range = Which files do you want to analyze in a specific day
%   bead_flag = Flag to determine comparison file 
%               (0 for NoScatAll, 1 for NoScatCell)
% Outputs:
%   successmessage = a string that indicates the completion of the code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialization
clc
mainFolder = filepath;
cd(mainFolder)

% Get the list of directories
dirinfo = dir();
dirinfo(~[dirinfo.isdir]) = [];  % Remove non-directories
dirinfo(ismember({dirinfo.name}, {'.', '..'})) = [];  % Remove . and ..
[~, c] = natsortfiles({dirinfo.name});
subdirinfo = cell(length(dirinfo));

% Find files ending with NoFLR.mat in each directory
for K = 1:length(dirinfo)
    thisdir = dirinfo(c(K)).name;
    subdirinfo{K} = dir(fullfile(thisdir, '*NoFLR.mat'));
end

subdirinfo = subdirinfo(~cellfun('isempty', subdirinfo));

%% Main Loop: Sort peaks with FLR and without FLR
for i = file_range
    disp(['Evaluating File # ', num2str(i), ' of ', num2str(size(subdirinfo, 1))]);
    peaksWithFLR = []; % Matrix of peaks with FLR
    peaksWithoutFLR = []; % Matrix of peaks without FLR
    
    if ~isempty(subdirinfo{i})
        cd(subdirinfo{i}(1).folder)
        load(subdirinfo{i}(1).name)
        pv1 = peak_values;
        chunks = unique(peak_values(:, 7));
        
        % Determine the comparison file based on bead_flag
        if bead_flag == 0
            comparisonFile = [subdirinfo{i}(1).name(1:end-9), 'NoScat.mat'];
        else
            comparisonFile = [subdirinfo{i}(1).name(1:end-9), 'NoScatCell.mat'];
        end
        load(comparisonFile, 'peak_values')
        pv2 = peak_values; 
        
        % Sort peaks
        for j = min(chunks):max(chunks)
            disp(['Chunk # ', num2str(j), ' of ', num2str(max(chunks))]);
            a = find(pv1(:, 7) == j);%a = pv1(:, 7) == j;
            loc1 = pv1(a, 8);
            start1 = pv1(a, 10);
            end1 = pv1(a, 11);
            a2 = pv2(:, 7) == j;%  a2 = find(pv2(:, 7) == j);
            loc2 = pv2(a2, 8);
            start2 = pv2(a2, 10);
            end2 = pv2(a2, 11);
            
            for k = 1:length(loc1) %for k = 1:length(loc2)
                l1 = loc1(k);      % l2 = loc2(k);
                % Check if any peak in loc1 falls within the range defined by start2 and end2
                match = find(loc2 >= start1(k) & loc2 <= end1(k), 1);%    match = find(loc1 >= start2(k) & loc1 <= end2(k), 1);
                if ~isempty(match)
                    peaksWithFLR = [peaksWithFLR; pv1(a(k), :)]; %#ok<AGROW> peaksWithFLR = [peaksWithFLR; pv2(a2(k), :)];
                else
                    peaksWithoutFLR = [peaksWithoutFLR; pv1(a(k), :)]; %#ok<AGROW>peaksWithoutFLR = [peaksWithoutFLR; pv2(a2(k), :)];
                end
            end
        end
        
        % Save sorted peaks
        peak_values = peaksWithFLR;
        save([subdirinfo{i}(1).name(1:end-9), 'NoFLRWithFLR.mat'], 'peak_values')
        peak_values = peaksWithoutFLR;
        save([subdirinfo{i}(1).name(1:end-9), 'NoFLRWithoutFLR.mat'], 'peak_values')
    end 
end

%% Save all sorted files
pv1 = [];
pv2 = [];
for i = file_range
    disp(['Saving sorted peaks for File # ', num2str(i), ' of ', num2str(size(subdirinfo, 1))]);
    if ~isempty(subdirinfo{i}(1))
        cd(subdirinfo{i}(1).folder)
        load([subdirinfo{i}(1).name(1:end-9), 'NoFLRWithFLR.mat'], 'peak_values')
        pv1 = [pv1; peak_values]; %#ok<AGROW>
        load([subdirinfo{i}(1).name(1:end-9), 'NoFLRWithoutFLR.mat'], 'peak_values')
        pv2 = [pv2; peak_values]; %#ok<AGROW>
    end
end
cd(mainFolder)
peak_values = pv1;
save([subdirinfo{file_range(1)}(1).name(4:end-9), 'NoFLRWithFLR.mat'], 'peak_values')
peak_values = pv2;
save([subdirinfo{file_range(1)}(1).name(4:end-9), 'NoFLRWithoutFLR.mat'], 'peak_values')

successmessage = 'Completed Sorting Peaks With and Without FLR';
end