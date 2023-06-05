%% DeepPeak Dataset Copier
%% Initialization
clear
clc
close all
alldata_path = 'T:\Nilay\IVFC\Acquired Data\Blood Cell Data';
deeppeak_path = 'T:\Nilay\IVFC\Acquired Data\Blood Cell Data\DeepPeak2023';
cd(deeppeak_path)
fname_dict = readtable("dict_key.csv",'Delimiter','''');
fname_dict = fname_dict(:,2);
%% Main Script
parfor i = 1:size(fname_dict,1)
    % Make a duplicate folder in DeepPeak2023 Path
    disp(['Working on copying ',fname_dict{i,1}{1}])
    mydir = fullfile([deeppeak_path,'\',fname_dict{i,1}{1}]);
    cd([deeppeak_path,'\',fname_dict{i,1}{1}])
    new_loc = cd;
    % Grab all folder names in orginal location
    folders = dir('*Cell*');
    folders(~[folders.isdir]) = []; 
    cd([alldata_path,'\',fname_dict{i,1}{1}])
    noFLR = dir('NEW*');

    for nflr = length(noFLR):-1:1
        if contains(noFLR(nflr).name,'_NoFLR.mat')==1
        else
            noFLR(nflr) = [];
        end
    end
    % Check that we only copy _noFLR files:
    if length(noFLR)>1
        peak_values_store = [];
        for nflr = 1:length(noFLR)
            pv = load(noFLR(nflr).name,'-mat','peak_values');
            pv = pv.peak_values;
            peak_values_store = [peak_values_store;pv];
        end
        filename = strrep(noFLR(1).name,'NEW','Old');
        temp_idx = strfind(filename,'Old')+3;
        if ~strcmp(filename(temp_idx),'_')
            filename(temp_idx) = [];
        end
        destination = fullfile(new_loc,filename);
        peak_values = peak_values_store;
        parsave(destination,peak_values)
    else
        for nflr = 1:length(noFLR)
            if contains(noFLR(nflr).name,'_NoFLR')==1
                source = fullfile(noFLR(nflr).folder,noFLR(nflr).name);
                filename = strrep(noFLR(nflr).name,'NEW','Old');
                destination = fullfile(new_loc,filename);
                copyfile(source,destination)
            end
        end
    end

    currentDir = cd;
    for n = 1:length(folders)
        cd(folders(n).name)
        noFLR = dir('*NEW*');
        for nflr = length(noFLR):-1:1
            if contains(noFLR(nflr).name,'_NoFLR')==1
            else
                noFLR(nflr) = [];
            end
        end

        tempvar = dir('*.csv');
        currentFoldername = tempvar(1).name(1:end-4);
        checkPath = fullfile([deeppeak_path,'\',fname_dict{i,1}{1},'\',...
                                                       currentFoldername]);
        if not(isfolder(checkPath))
            mkdir(checkPath)
        end
        % Check that we only copy _noFLR files:
        for nflr = 1:length(noFLR)
            if contains(noFLR(nflr).name,'_NoFLR.mat')==1 && ~contains(noFLR(nflr).name,'Thresh')
                source = fullfile(noFLR(nflr).folder,noFLR(nflr).name);
                filename = strrep(noFLR(nflr).name,'NEW','Old');
                temp_idx = strfind(filename,'Old')+3;
                if ~strcmp(filename(temp_idx),'_')
                    filename(temp_idx) = [];
                end
                destination = fullfile(checkPath,filename);
                copyfile(source,destination)
            end
        end
        cd(currentDir)
    end
    cd(deeppeak_path)
end
%% NoScat
parfor i = 1:size(fname_dict,1)
    % Make a duplicate folder in DeepPeak2023 Path
    disp(['Working on copying ',fname_dict{i,1}{1}])
    mydir = fullfile([deeppeak_path,'\',fname_dict{i,1}{1}]);
    cd([deeppeak_path,'\',fname_dict{i,1}{1}])
    new_loc = cd;
    % Grab all folder names in orginal location
    folders = dir('*Cell*');
    folders(~[folders.isdir]) = []; 
    cd([alldata_path,'\',fname_dict{i,1}{1}])
    NoScat = dir('NEW*');

    for nscat = length(NoScat):-1:1
        if contains(NoScat(nscat).name,'_NoScatCell')==1 || contains(NoScat(nscat).name,'_NoScat.mat')==1
        else
            NoScat(nscat) = [];
        end
    end
    % Check that we only copy _NoScat files:
    if length(NoScat)>1
        peak_values_store = [];
        for nscat = 1:length(NoScat)
            pv = load(NoScat(nscat).name,'-mat','peak_values');
            pv = pv.peak_values;
            peak_values_store = [peak_values_store;pv];
        end
        filename = strrep(NoScat(1).name,'NEW','Old');
        temp_idx = strfind(filename,'Old')+3;
        if ~strcmp(filename(temp_idx),'_')
            filename(temp_idx) = [];
        end
        destination = fullfile(new_loc,filename);
        peak_values = peak_values_store;
        parsave(destination,peak_values)
    else
        for nscat = 1:length(NoScat)
            if contains(NoScat(nscat).name,'_NoScat')==1
                source = fullfile(NoScat(nscat).folder,NoScat(nscat).name);
                filename = strrep(NoScat(nscat).name,'NEW','Old');
                destination = fullfile(new_loc,filename);
                copyfile(source,destination)
            end
        end
    end

    currentDir = cd;
    for n = 1:length(folders)
        cd(folders(n).name)
        NoScat = dir('*NEW*');
        for nscat = length(NoScat):-1:1
            if contains(NoScat(nscat).name,'_NoScatCell')==1 || contains(NoScat(nscat).name,'_NoScat.mat')==1
            else
                NoScat(nscat) = [];
            end
        end

        tempvar = dir('*.csv');
        currentFoldername = tempvar(1).name(1:end-4);
        checkPath = fullfile([deeppeak_path,'\',fname_dict{i,1}{1},'\',...
                                                       currentFoldername]);
        if not(isfolder(checkPath))
            mkdir(checkPath)
        end
        % Check that we only copy _NoScat files:
        for nscat = 1:length(NoScat)
            if contains(NoScat(nscat).name,'_NoScat')==1 && ~contains(NoScat(nscat).name,'Thresh')
                source = fullfile(NoScat(nscat).folder,NoScat(nscat).name);
                filename = strrep(NoScat(nscat).name,'NEW','Old');
                temp_idx = strfind(filename,'Old')+3;
                if ~strcmp(filename(temp_idx),'_')
                    filename(temp_idx) = [];
                end
                destination = fullfile(checkPath,filename);
                copyfile(source,destination)
            end
        end
        cd(currentDir)
    end
    cd(deeppeak_path)
end
%% Saving fcn
function parsave(fname, peak_values)
  save(fname, 'peak_values')
end