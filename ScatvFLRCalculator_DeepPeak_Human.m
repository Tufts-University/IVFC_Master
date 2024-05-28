function [finalSens, finalPur] = ScatvFLRCalculator_DeepPeak_Human(mainFolder,header)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  file: ScatvFLRCalculator_DeepPeak_Human.m
%  authors: Nilay Vora
%  created: 06/04/2023
% modified:
%  purpose: This file will read peak data and correlate it against
%  detection schemes for peaks. Namely, Scat vs. FLR vs. Scat+FLR. This
%  data is then returned as a table to be saved per folder as a .xlsx.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
datasheet_path = 'T:\Taras\IVFC\Acquired Data\Cell\Test';% path was changed
cd(mainFolder)
T = readtable([datasheet_path,'\TH_20230705_DatasetSummary.xlsx']);
fname_dict = T.FolderName;
fname_dict = cellfun(@(x) x(2:end-1), fname_dict, 'UniformOutput', false);
CTCCFlag = T.CTCCFlag;
BeadFlag = T.BeadFlag;
idx = find(CTCCFlag==1 & T.inc==1);
Final_Fname = {fname_dict{idx}}';
Final_BeadFlag = BeadFlag(idx);

total_tot = [];
count_tot = [];
SC_tot = [];
SC_tot2 = [];
total_tot2 = [];
count_tot2 = [];
peaks_FLR = [];
peaks_Scat = [];
peaks_FLR2 = [];
peaks_Scat2 = [];
tempnames = {};
%% Interogation
bead_flags = Final_BeadFlag;
for i=1:length(Final_Fname)
    disp(['Day ',num2str(i), ' of ',num2str(length(Final_Fname))])
    newFolder = Final_Fname{i};
    cd(newFolder)
    newFolder = cd;
    dirinfo2 = dir();
    dirinfo2(~[dirinfo2.isdir]) = [];  %remove non-directories

    dirinfo2(ismember( {dirinfo2.name}, {'.', '..','Local Plots'})) = [];  %remove . and ..

    bead_flag = bead_flags(i);
    subdirinfo2 =  dirinfo2;
    if bead_flag == 1
        search_tag = '_NoScatCell';
    else
        search_tag = '_NoScat.mat';
    end
    for j=1:length(subdirinfo2)
        disp(['     File ',num2str(j),' of ',num2str(length(subdirinfo2))])
        newFolder2 = [subdirinfo2(j).folder,'\',subdirinfo2(j).name];
        cd(newFolder2)
        base_name = [];
        files = [dir('*21*'); dir('*22*'); dir('*23*'); dir('*24*')];
        if length(files)>1
            for k=1:length(files)
                if contains(files(k).name,header) && contains(files(k).name,search_tag)
                    base_name=files(k).name;
                    break
                end
            end
            for k=1:length(files)
                if contains(files(k).name,header) && contains(files(k).name,'_NoFLR')
                    base_name2=files(k).name;
                    break
                end
            end
            if isempty(base_name) == 1
                for k=1:length(files)
                    if contains(files(k).name,header) && contains(files(k).name,'_NoScat')
                        base_name=files(k).name;
                        break
                    end
                end
            end
        end
        [count,total,SC,peaks1]=NV_070823_ScatvFLRML_6PMT(newFolder2,base_name,base_name2,bead_flag);
        [count2,total2,SC2,peaks2]=NV_070823_ScatvFLRML_6PMT(newFolder2,base_name2,base_name,bead_flag);
        SC_tot=[SC_tot;[SC,SC2]];
        total_tot=[total_tot;[total,total2]];
        count_tot=[count_tot;[count,count2]];
        peaks_FLR = [peaks_FLR;[repmat(i,size(peaks1,1),1),repmat(j,size(peaks1,1),1),peaks1]];
        peaks_Scat = [peaks_Scat;[repmat(i,size(peaks2,1),1),repmat(j,size(peaks2,1),1),peaks2]];
        tempnames = [tempnames;subdirinfo2(j).name];
    end
    cd(newFolder)
    files = [dir('*21*'); dir('*22*'); dir('*23*');dir('*24*')];
    base_name = [];
    if ~isempty(files)
        for k=1:length(files)
            if contains(files(k).name,header) && contains(files(k).name,search_tag)
                base_name=files(k).name;
                break
            end
        end
        for k=1:length(files)
            if contains(files(k).name,header) && contains(files(k).name,'_NoFLR')
                base_name2=files(k).name;
                break
            end
        end
        if isempty(base_name) == 1
            for k=1:length(files)
                if contains(files(k).name,header) && contains(files(k).name,'_NoScat')
                    base_name=files(k).name;
                    break
                end
            end
        end

    end
    [count3,total3,SC3,peaks3]=NV_070823_ScatvFLRML_6PMT(newFolder,base_name,base_name2,bead_flag);
    [count4,total4,SC4,peaks4]=NV_070823_ScatvFLRML_6PMT(newFolder,base_name2,base_name,bead_flag);
    SC_tot2 = [SC_tot2;[SC3,SC4]];
    total_tot2 = [total_tot2;[total3,total4]];
    count_tot2 = [count_tot2;[count3,count4]];
    peaks_FLR2 = [peaks_FLR2;[repmat(i,size(peaks3,1),1),peaks3]];
    peaks_Scat2 = [peaks_Scat2;[repmat(i,size(peaks4,1),1),peaks4]];
    cd(mainFolder)
end

t=total_tot2./sum(total_tot2);
perform_mat = [(1:length(Final_Fname))',SC_tot2,t];
T_final = table(Final_Fname,SC_tot2);
writetable(T_final,'SensitivityTable_Human.xlsx')
t2 = total_tot./sum(total_tot);

perform_mat2 = [SC_tot,t2];
T_final_full = table(tempnames,perform_mat2);
writetable(T_final_full,'SensitivityTable_Human_Full.xlsx')
disp('Weighted Performance')
disp('------------------------------')
final = sum(SC_tot2.*t,'omitnan')./sum(t)
disp('------------------------------')
finalSens = final(1);
finalPur = final(2);