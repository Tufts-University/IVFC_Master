%% 063023_sort
% Owner: Taras Hanulia
% Data Type: Blood Cell Data from Human plus Cell
% Flow Date:  
%% Notes
% Medford 1
%% Initialization
clear
clc
addpath 'C:\Users\thanul01\Documents\MATLAB\ivfc_master'
%% Calling Script

    filepath = 'T:\Taras\IVFC\Acquired Data\Cell\TH_051524_Cell';
    % Labview Conversion
    Fs=60e3;
    Window_Low= 50;
    Window_High= 10000;
   % output=Labview_convert_rawdata_batch_6PMT(filepath,Fs,Window_Low,Window_High);
   % disp(output)
    close all
    
    cd(filepath)
    dirinfo = dir();
    dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
    dirinfo(ismember( {dirinfo.name}, {'.', '..'})) = [];  %remove . and ..
    exp_name = dirinfo(1).name;
    date = exp_name(4:9);
    date = strread(date,'%2s');
    %%Peak Detection
    outputfile= ['NEW_peak_values_',date{1},'_',date{2},'_',date{3}];
    dirinfo = dir('*Cell*');
    dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
    dirinfo(ismember( {dirinfo.name}, {'.', '..'})) = [];  %remove . and ..
    file_range= (1:length(dirinfo));
    sample_type= 'Cell';
    exp_num=[];
    std_threshold=4;
    Spectralon_tail= '';
    FWMH_threshold=0;
    intensity_threshold= 0.1;
    bead_flag=0;
    output=FLRSorting6PMT(filepath,file_range);
    disp(output)