function [finalSens, finalPur] = ScatvFLRCalculator(mainFolder,header)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  file: NV_100722_ScatvFLRML_batch.m
%  authors: Nilay Vora
%  created: 10/07/22
% modified:
%  purpose: This file will read peak data and correlate it against
%  detection schemes for peaks. Namely, Scat vs. FLR vs. Scat+FLR. This
%  data is then returned as a table to be saved per folder as a .xlsx.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
basetail='*22_*'; % This is the last number on the peak_values .mat files
% we need this to grab the base file name for iterative interrogation of
% each data type.
basetail_21='*21*';
%%
cd(mainFolder)
pv=[];

dirinfo = dir();
dirinfo(~[dirinfo.isdir]) = [];  %remove non-directories

dirinfo(ismember( {dirinfo.name}, {'.', '..'})) = [];  %remove . and ..
% dirinfo(ismember( {dirinfo.name}, {'NV_030822_Blood_Cell', ...
%                                    'NV_031522_Blood_Cell',...
%                                    'NV_032122_Blood_Cell',...
%                                    'NV_040422_Blood_Cell'})) = [];  %remove . and ..

subdirinfo = cell(length(dirinfo));

for K = 1 : length(dirinfo)
    thisdir = dirinfo(K).name;
    subdirinfo{K} = dir(fullfile(thisdir, '*_NoScat*'));
end

subdirinfo =  subdirinfo(~cellfun('isempty',subdirinfo));

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
%% Interogation
bead_flags = ones(length(subdirinfo),1);
bead_flags(6:13) = 0;
%bead_flags(12) = 1;
for i=1:length(subdirinfo)
    disp(['Day ',num2str(i), ' of ',num2str(length(subdirinfo))])
    newFolder=subdirinfo{i}.folder;
    cd(newFolder)

    dirinfo2 = dir();
    dirinfo2(~[dirinfo2.isdir]) = [];  %remove non-directories

    dirinfo2(ismember( {dirinfo2.name}, {'.', '..','Local Plots'})) = [];  %remove . and ..

    bead_flag = bead_flags(i);
    subdirinfo2 =  dirinfo2;
    for j=1:length(subdirinfo2)
        disp(['     File ',num2str(j),' of ',num2str(length(subdirinfo2))])
        newFolder2=[subdirinfo2(j).folder,'\',subdirinfo2(j).name];
        cd(newFolder2)
        base_name = [];
        files= dir(basetail);
        if length(files)>1
            for k=1:length(files)
                if contains(files(k).name,'Corrected') && contains(files(k).name,'_NoScatCell')
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
                    if contains(files(k).name,'Corrected') && contains(files(k).name,'_NoScat')
                        base_name=files(k).name;
                        break
                    end
                end
            end
        else
            files = dir(basetail_21);
            if length(files)>1
                for k=1:length(files)
                    if contains(files(k).name,'Corrected') ...
                            && contains(files(k).name,'_NoScatCell')
                        base_name=files(k).name;
                        break
                    end
                end
                for k=1:length(files)
                    if contains(files(k).name,header) ...
                            && contains(files(k).name,'_NoFLR')
                        base_name2=files(k).name;
                        break
                    end
                end
                if isempty(base_name) == 1
                    for k=1:length(files)
                        if contains(files(k).name,'Corrected') ...
                                && contains(files(k).name,'_NoScat')
                            base_name=files(k).name;
                            break
                        end
                    end
                end
            end
        end
        [count,total,SC,peaks1]=NV_020323_ScatvFLRML(newFolder2,base_name,base_name2,bead_flag);
        [count2,total2,SC2,peaks2]=NV_020323_ScatvFLRML(newFolder2,base_name2,base_name,bead_flag);
        SC_tot=[SC_tot;[SC,SC2]];
        total_tot=[total_tot;[total,total2]];
        count_tot=[count_tot;[count,count2]];
        peaks_FLR = [peaks_FLR;[repmat(i,size(peaks1,1),1),repmat(j,size(peaks1,1),1),peaks1]];
        peaks_Scat = [peaks_Scat;[repmat(i,size(peaks2,1),1),repmat(j,size(peaks2,1),1),peaks2]];
    end
    cd(newFolder)
    files= dir(basetail);
    base_name = [];
    if ~isempty(files)
        for k=1:length(files)
            if contains(files(k).name,'Corrected') && contains(files(k).name,'_NoScatCell')
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
                if contains(files(k).name,'Corrected') && contains(files(k).name,'_NoScat')
                    base_name=files(k).name;
                    break
                end
            end
        end
    else
        files = dir(basetail_21);
        for k=1:length(files)
            if contains(files(k).name,'Corrected') ...
                    && contains(files(k).name,'_NoScatCell')
                base_name=files(k).name;
                break
            end
        end
        for k=1:length(files)
            if contains(files(k).name,header) ...
                    && contains(files(k).name,'_NoFLR')
                base_name2=files(k).name;
                break
            end
        end
        if isempty(base_name)
            for k=1:length(files)
                if contains(files(k).name,'Corrected') ...
                        && contains(files(k).name,'_NoScat')
                    base_name=files(k).name;
                    break
                end
            end
        end
    end
    [count3,total3,SC3,peaks3]=NV_020323_ScatvFLRML(newFolder,base_name,base_name2,bead_flag);
    [count4,total4,SC4,peaks4]=NV_020323_ScatvFLRML(newFolder,base_name2,base_name,bead_flag);
    SC_tot2 = [SC_tot2;[SC3,SC4]];
    total_tot2 = [total_tot2;[total3,total4]];
    count_tot2 = [count_tot2;[count3,count4]];
    peaks_FLR2 = [peaks_FLR2;[repmat(i,size(peaks3,1),1),peaks3]];
    peaks_Scat2 = [peaks_Scat2;[repmat(i,size(peaks4,1),1),peaks4]];
end

t=total_tot2./sum(total_tot2);
perform_mat = [(1:30)',SC_tot2,t];
disp('Weighted Performance')
disp('------------------------------')
final = sum(SC_tot2.*t,'omitnan')./sum(t);
disp('------------------------------')
finalSens = final(1);
finalPur = final(2);
% sortrat = sortrows(perform_mat,2);
% disp('Weighted Performance of Best 26 days')
% disp('------------------------------')
% final = sum(sortrat(5:end,2:3).*sortrat(5:end,4:5), 'omitnan')./(sum(sortrat(5:end,4:5)))
% disp('------------------------------')
% % By data file
% t=total_tot./sum(total_tot);
% perform_mat = [(1:92)',SC_tot,t];
% disp('Weighted Performance')
% disp('------------------------------')
% final = sum(SC_tot.*t, 'omitnan')./sum(t)
% disp('------------------------------')
% 
% total_mat = [peaks_Scat2;peaks_FLR2];
% unique_totals = unique(total_mat,'rows');
% total_FLR = sum(total_tot2);
% disp('True Performance')
% disp('------------------------------')
% final = length(unique_totals)./total_FLR(1,1)
% disp('------------------------------')