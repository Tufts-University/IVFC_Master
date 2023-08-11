%% SciRep 2022 Paper Data Formatting
%% Initialization
clear
clc
%% Formatting
mainpath='T:\Nilay\IVFC\Acquired Data\Blood Cell Data\DeepPeak2023';
cd(mainpath)
file_folders= '*Blood*';
folders = dir(file_folders);
folders(~[folders.isdir]) = [];  %remove non-directories
% TP_Peak_ranges = struct();
% for i=1:length(folders)
%     cd([mainpath,'\',folders(i).name])
%     disp(['Day # ',num2str(i),' of ',num2str(length(folders))])
%     scat_file = dir('*ScatCell.mat');
%     if isempty(scat_file)
%         scat_file = dir('*Scat.mat');
%     end
%     load(scat_file.name)
%     [range_data,label,peak_values] = ...
%         MLformating(scat_file.name);
%     TP_Peak_ranges.(folders(i).name) = [range_data,label];
% end
%% Scattering datasets Corrected
FP_Peak_ranges = struct();
for i = 1:length(folders)
    FP_Peak_ranges.(folders(i).name) = [];
end
tempstorage = cell(1,length(folders));
bead_flag = ones(length(folders),1);
bead_flag(10:17) = 0;

parfor i=1:length(folders)
    cd([mainpath,'\',folders(i).name])
    disp(['Day # ',num2str(i),' of ',num2str(length(folders))])
    scat_file = dir('*NoFLR.mat');
    [range_data,label,peak_values] = ...
        MLformating(scat_file(1).name,bead_flag(i));
    lbead=find(label==3); %Remove any beads from data
    label(lbead)=[];
    range_data(lbead,:)=[];
    peak_values(lbead,:)=[];
    tempstorage{i} = [range_data,label];
end
%% Save with keys
for i=1:length(folders)
   FP_Peak_ranges.(folders(i).name) = tempstorage{i};
end
%% Scattering datasets
% I deleted the NOFLR data from CNNData so we're going to try something
% fancy using a interative search algorithm... I hope(?)
% all_data_path = 'T:\Nilay\IVFC\Acquired Data\Blood Cell Data';
% FP_Peak_ranges = struct();
% for i=1:length(folders)
%     cd([mainpath,'\',folders(i).name])
%     subfolder = dir([mainpath,'\',folders(i).name]);
%     subfolder(~[subfolder.isdir]) = [];  %remove non-directories
%     subfolder(ismember( {subfolder.name}, {'.', '..','Local Plots'})) = [];
%     cd(subfolder(1).name)
%     root_name = dir('*-Corrected*');
%     idx_ = find(root_name(1).name=='_');
%     fname = [root_name(1).name(4:idx_(end)-1),'_NoFLR.mat'];
%     %%%%%%%%%%%%%%% Added for the new scat files %%%%%%%%%%%%%%%%%%%%%%%%%%
%     cd([mainpath,'\',folders(i).name])
%     scat_file = dir(fname);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     cd([all_data_path,'\',folders(i).name])
%     disp(['Day # ',num2str(i),' of ',num2str(length(folders))])
%     %%%%%%%%%%%%%%%Removed for the new scat files %%%%%%%%%%%%%%%%%%%%%%%%%
%     % scat_file = dir(fname);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     load([scat_file(end).folder,'\',scat_file(end).name])
%     [range_data,label,peak_values] = ...
%         MLformating([scat_file(end).folder,'\',scat_file(end).name]);
%     lbead=find(label==3); %Remove any beads from data
%     label(lbead)=[];
%     range_data(lbead,:)=[];
%     peak_values(lbead,:)=[];
%
%     lTP=find(label==1); %Remove TPs from data
%     label(lTP)=[];
%     range_data(lTP,:)=[];
%     peak_values(lTP,:)=[];
%     FP_Peak_ranges.(folders(i).name) = [range_data,label];
% end
%% Data Formatting
function [range_data,label,peak_values]=MLformating(fname,bead_flag)
load(fname,'-mat','peak_values')
peak_values(peak_values(:,9)<17,:) = [];
peak_values = sortrows(peak_values,[21,6,7]);
range_data = zeros(size(peak_values,1),396);
label = 1.*ones(size(peak_values,1),1);
count = 1;

dirinfo = dir();
dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
dirinfo(ismember( {dirinfo.name}, {'.', '..','Local Plots'})) = [];  %remove . and ..
[~,c]=natsortfiles({dirinfo.name});
dirinfo = dirinfo(c);

for j=1:length(unique(peak_values(:,21)))
    disp(['     Subfile # ',num2str(j),' of ',...
        num2str(length(unique(peak_values(:,21))))])
    unique_num = unique(peak_values(:,21));
    pv_temp = peak_values(peak_values(:,21)==unique_num(j),:);
    cd([dirinfo(unique_num(j)).folder,'\',dirinfo(unique_num(j)).name])
    chunks = dir('*raw.mat');
    [~,c]=natsortfiles({chunks.name});
    for k=1:size(chunks,1)
        disp(['          Chunk # ',num2str(k),' of ',...
            num2str(size(chunks,1))])
        load(chunks(c(k)).name)
        M(:,1)=(M(:,1)-mean(M(:,1)))./std(M(:,1));
        M(:,2)=(M(:,2)-mean(M(:,2)))./std(M(:,2));
        M(:,3)=(M(:,3)-mean(M(:,3)))./std(M(:,3));
        chunk_loc=pv_temp(pv_temp(:,6)==k,:);
        for l=1:size(chunk_loc,1)
            loc2=chunk_loc(l,7);
            M_range=[M(loc2-49:loc2+49,1)',M(loc2-49:loc2+49,2)',...
                M(loc2-49:loc2+49,3)',M(loc2-49:loc2+49,5)'];
            if bead_flag == 1
                if chunk_loc(l,4)<0.05 & chunk_loc(l,5)>0.05
                    if chunk_loc(l,19)>20
                        label(count)=0; % CTCCs
                    else
                        label(count)= 0; %CTC
                    end
                elseif chunk_loc(l,4)>0.05
                    label(count)=3; % Bead
                end
            else
                if chunk_loc(l,5)>0.05
                    if chunk_loc(l,19)>30
                        label(count)=0; % CTCC
                    else
                        label(count)= 0; %CTC
                    end
                end
            end
            range_data(count,:)=M_range;
            count=count+1;
        end
    end
end
end