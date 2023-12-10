%% DeepPeak Paper Data Formatting CAL27EGFP+ Cells
%% Initialization
clear
clc
%% Formatting
mainpath=['T:\Taras\IVFC\Acquired Data\Blood Cell Data'];
cd(mainpath)
file_folders= '*Blood*';
folders = dir(file_folders);
folders(~[folders.isdir]) = [];  %remove non-directories
%% Scattering datasets Corrected
FP_Peak_ranges = struct();
for i = 1:length(folders)
    FP_Peak_ranges.(folders(i).name) = [];
end
file_name = ['T:\Nilay\IVFC\Acquired Data\Blood Cell Data\',...
                                'NV_20230328_DatasetSummary.xlsx'];
T = readtable(file_name);
fname_dict = T.FolderName;
fname_dict = cellfun(@(x) x(2:end-1), fname_dict, 'UniformOutput', false);
CTCCFlag = T.CTCCFlag;
BeadFlag = T.BeadFlag;
RatFlag = T.RatFlag;
inc = T.inc;
cellType = T.Cell231;
idx = find(RatFlag==1 & CTCCFlag==1 & inc==1 & cellType==0);
Final_Fname = {fname_dict(idx)}';
Final_BeadFlag = BeadFlag(idx);

tempstorage = cell(1,length(folders));
bead_flag = Final_BeadFlag;
n = length(folders);
folders = Final_Fname{1};
parfor i=1:length(folders)
    cd([mainpath,'\',folders{i}])
    disp(['Day # ',num2str(i),' of ',num2str(n)])
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
   FP_Peak_ranges.(folders{i}) = tempstorage{i};
end
%% Data Formatting
function [range_data,label,peak_values]=MLformating(fname,bead_flag)
load(fname,'-mat','peak_values')
peak_values(peak_values(:,10)<17,:) = [];
peak_values = sortrows(peak_values,[24,7,8]);
range_data = zeros(size(peak_values,1),396);
label = 1.*ones(size(peak_values,1),1);
count = 1;

dirinfo = dir();
dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories
dirinfo(ismember( {dirinfo.name}, {'.', '..','Local Plots'})) = [];  %remove . and ..
[~,c]=natsortfiles({dirinfo.name});
dirinfo = dirinfo(c);

for j=1:length(unique(peak_values(:,24)))
    disp(['     Subfile # ',num2str(j),' of ',...
        num2str(length(unique(peak_values(:,24))))])
    unique_num = unique(peak_values(:,24));
    pv_temp = peak_values(peak_values(:,24)==unique_num(j),:);
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
        chunk_loc=pv_temp(pv_temp(:,7)==k,:);
        for l=1:size(chunk_loc,1)
            loc2=chunk_loc(l,8);
            M_range=[M(loc2-49:loc2+49,1)',M(loc2-49:loc2+49,2)',...
                M(loc2-49:loc2+49,3)',M(loc2-49:loc2+49,5)'];
            if bead_flag == 1
                if chunk_loc(l,4)<0.05 && chunk_loc(l,5)>0.05
                    label(count)=0; % CTCCs
                elseif chunk_loc(l,4)>0.05
                    label(count)=3; % Bead
                end
            else
                if chunk_loc(l,5)>0.05
                    label(count)=0; % CTCC
                end
            end
            range_data(count,:)=M_range;
            count=count+1;
        end
    end
end
end