%% DeepPeak Paper Data Formatting
%% Initialization
clear
clc
%% Formatting
mainpath=['T:\Nilay\IVFC\Acquired Data\Blood Cell Data\DeepPeak2023\'];
cd(mainpath)
file_folders= '*Blood*';
folders = dir(file_folders);
folders(~[folders.isdir]) = [];  %remove non-directories
%% Scattering datasets Corrected
Loc_Peak_ranges = struct();
for i = 1:length(folders)
    Loc_Peak_ranges.(folders(i).name) = [];
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
idx = find(RatFlag==1 & CTCCFlag==1 & inc==1);
Final_Fname = {fname_dict(idx)}';
Final_BeadFlag = BeadFlag(idx);

tempstorage = cell(1,length(folders));
bead_flag = Final_BeadFlag;
n = length(folders);
folders = Final_Fname{1};
for threshold = 100%  [1,10,30,50,100]
    disp(['Running Threshold: ',num2str(threshold)])
    parfor i=1:length(folders)
        cd([mainpath,'\',folders{i}])
        disp(['Day # ',num2str(i),' of ',num2str(n)])
        scat_file = dir('*NoFLR.mat');
        day = i;
        [locs_save,label] = ...
            Locformating(scat_file(1).name,bead_flag(i),day);
        lbead=find(label==3); %Remove any beads from data
        label(lbead)=[];
        locs_save(lbead,:)=[];
        allCTCCs = find(label==0);
        try
            locs_save = locs_save(1:allCTCCs(threshold),:);
            label = label(1:allCTCCs(threshold),:);
        catch
            locs_save = locs_save(1:allCTCCs(end),:);
            label = label(1:allCTCCs(end),:);
            disp(folders{i})
        end

        tempstorage{i} = [locs_save,label];
    end
    %% Save with keys
    for i=1:length(folders)
       Loc_Peak_ranges.(folders{i}) = tempstorage{i};
    end
    filename = ['NV_032823_FormattedDataSet_thresh',...
        num2str(threshold,'%03.0f'),'_locs.mat'];
    save(filename,"Loc_Peak_ranges")
end
%% Data Formatting
function [locs_save,label]=Locformating(fname,bead_flag,day)
load(fname,'-mat','peak_values')
peak_values(peak_values(:,9)<17,:) = [];
peak_values = sortrows(peak_values,[21,6,7]);
label = ones(size(peak_values,1),1);
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
            count=count+1;
        end
    end
end
locs_save = [repmat(day,size(peak_values,1),1),peak_values(:,21),...
    peak_values(:,6),peak_values(:,7)];
end