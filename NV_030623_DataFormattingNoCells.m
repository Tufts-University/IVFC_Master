%% SciRep 2022 Paper Data Formatting
%% Initialization
clear
clc
%% Formatting
mainpath='T:\Nilay\IVFC\Acquired Data\Blood Cell Data';
cd(mainpath)
file_folders= '*Blood_noCell*';
folders = dir(file_folders);
folders(~[folders.isdir]) = [];  %remove non-directories
folders(5) = [];

%% Scattering datasets Corrected
FP_Peak_ranges = struct();
bead_flag = ones(length(folders),1).*0;
tempstorage = cell(1,length(folders));
% bead_flag(6:14) = 0;
for i=1:length(folders)
    cd([mainpath,'\',folders(i).name])
    disp(['Day # ',num2str(i),' of ',num2str(length(folders))])
    scat_file = dir('*NoFLR.mat');
    load(scat_file(1).name)
    [range_data,label,peak_values] = ...
        MLformating(scat_file(1).name,bead_flag(i));
    lbead=find(label==3); %Remove any beads from data
    label(lbead)=[];
    range_data(lbead,:)=[];
    peak_values(lbead,:)=[];
    label = ones(size(peak_values,1),1);
    FP_Peak_ranges.(folders(i).name) = [range_data,label];
end
%% Data Formatting
function [range_data,label,peak_values]=MLformating(fname,bead_flag)
load(fname,'-mat','peak_values')
peak_values(peak_values(:,9)<20,:) = [];
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
                if chunk_loc(l,4)<0.05 & chunk_loc(l,5)>0.15
                    if chunk_loc(l,19)>20
                        label(count)=0; % CTCC
                    else
                        label(count)= 0; %s
                    end
                elseif chunk_loc(l,4)>0.05
                    label(count)=3; % Bead
                end
            else
                if chunk_loc(l,5)>0.05
                    if chunk_loc(l,19)>20
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