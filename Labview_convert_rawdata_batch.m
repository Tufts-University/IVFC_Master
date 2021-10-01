function [successmessage] = Labview_convert_rawdata_batch(filepath,Fs,Window_Low,Window_High)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file:Labview_convert_rawdata_batch.m
% ***Description***:
% This function serves to convert any .csv file found in the parent folder
% into .mat files. This script automatically loops through the .csv file
% and save the chunks in their respective folders. 
% Written By: Nilay Vora (nvora01@tufts.edu)
% Date Written: 10/01/2021
% Modifying Author:
% Date Modified:
% Latest Revision: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function details
% Inputs:
%   filepath = Location of all main folder where all subfolders are
%   Fs = Sample frequency used during acquisition
%   Window_Low = Low frequency cutoff of Butterworth filter window
%   Window_High= High frequency cutoff of Butterworth filter window
% Outputs:
%   successmessage = a string that indicates the completion of the code
%
% Usage Labview_convert_rawdata_batch(filepath,Fs,Window)
% Example:
% filepath = 'U:\Nilay\IVFC\Acquired Data\Blood Data\NV_092821_Blood_LNPs';
% Fs=60e3 %60,000 samples per second
% Window_Low= 50; % removes frequencies below 50/30000 Hz
% Window_High= 6000; % removes frequencies above 6000/30000 Hz
% output=Labview_convert_rawdata_batch(filepath,Fs,Window_Low,Window_High)
%% Checking inputs
if nargin < 2
    Fs=60e3;
    Window_Low=50;
    Window_High=6000;
end

if nargin < 3
    Window_Low=50;
    Window_High=6000;
end

if nargin < 4
    Window_High=6000;
end
%% Finding Folders with Files
mainFolder=filepath;

cd(mainFolder)

dirinfo = dir();
dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories

dirinfo(ismember( {dirinfo.name}, {'.', '..'})) = [];  %remove . and ..

subdirinfo = cell(length(dirinfo));

for K = 1 : length(dirinfo)
    thisdir = dirinfo(K).name;
    subdirinfo{K} = dir(fullfile(thisdir, '*.csv'));
end

subdirinfo =  subdirinfo(~cellfun('isempty',subdirinfo));
%% Main Loop
Window=[Window_Low,Window_High];
for i=1:length(subdirinfo)
    disp(['Evaluating File # ',num2str(i),' of ',num2str(size(subdirinfo,1))]);
    if~isempty(subdirinfo{i})
        scatpath=fullfile(subdirinfo{i}.folder,subdirinfo{i}.name);
        cd(subdirinfo{i}.folder)
        data_type=subdirinfo{i}.name(1:end-4);
        
        %% Converts Labview processed CSV text file of RAW Data to Mat files
        feature accel on
        format shortEng
        chunk_time=1.5;
        smpl_rate=Fs; %per channel
        chunk_size=smpl_rate*60*chunk_time; 
        
        Wn=Window./(Fs/2);%Cutoff frequencies divided by Nyquist frequency
        [b,a]=butter(2,Wn);
        error=0;
        fid=fopen(scatpath);
        ii=1;
        while 1
            disp(['Chunk # ',num2str(ii)])
            block=chunk_size;
            [M,~]=textscan(fid,'%f %f %f %f %f',block,'Delimiter',',');
            M=cell2mat(M);
            if isempty(M)
                %Nothing happens if that chunk is empty
            else
            M_filt=zeros(length(M),5);
            avg_scat(ii,:)=mean(M,'omitnan');
            std_scat(ii,:)=std(M,'omitnan');
            scat_file=[data_type,'_',num2str(ii),'_raw'];

            M_filt=zeros(length(M),5);
            save(scat_file,'M')
            M_filt(:,1)=filtfilt(b,a,M(:,1));
            M_filt(:,2)=filtfilt(b,a,M(:,2));
            M_filt(:,3)=filtfilt(b,a,M(:,3));
            M_filt(:,4)=filtfilt(b,a,M(:,4));
            M_filt(:,5)=filtfilt(b,a,M(:,5));
            sigmas(ii,:)=std(M_filt,0,1);
            if size(M,1)<chunk_size
                break
            end
            ii=ii+1;
            clear M M_filt
            end
        end
        num_chunks=ii;
        sigmas_final=mean(sigmas,1);
        save([data_type,'_sigmas'],'sigmas_final');
        
        h=figure;
        subplot(2,3,1)
        errorbar(1:num_chunks-error,avg_scat(:,1),std_scat(:,1));
        title('Assesing Channel Drift','FontSize',12)
        xlabel('Time(1.5min chunks)','FontSize',10)
        ylabel('Average 405 Intensity','FontSize',10)
        
        subplot(2,3,2)
        errorbar(1:num_chunks-error,avg_scat(:,2),std_scat(:,2));
        xlabel('Time(1.5min chunks)','FontSize',10)
        ylabel('Average 488 Intensity','FontSize',10)
        
        subplot(2,3,3)
        errorbar(1:num_chunks-error,avg_scat(:,3),std_scat(:,3));
        xlabel('Time(1.5min chunks)','FontSize',10)
        ylabel('Average 633 Intensity','FontSize',10)
        
        subplot(2,3,4)
        errorbar(1:num_chunks-error,avg_scat(:,4),std_scat(:,4));
        xlabel('Time(1.5min chunks)','FontSize',10)
        ylabel('Average Red Flr Intensity','FontSize',10)
        
        subplot(2,3,5)
        errorbar(1:num_chunks-error,avg_scat(:,5),std_scat(:,5));
        xlabel('Time(1.5min chunks)','FontSize',10)
        ylabel('Average Green Flr Intensity','FontSize',10)
        
        set(gcf,'Position',[680.0000e+000   559.0000e+000   798.0000e+000   419.0000e+000]);
        filename=[data_type,'_Channel_Drift.jpeg'];
        saveas(h,filename)
    end
    fclose(fid);
    clear avg_scat std_scat
end

cd(mainFolder)
successmessage='Completed coversion of all .csv files to .mat files';
end