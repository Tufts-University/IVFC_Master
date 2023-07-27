function [count_mat,tot_mat,ratio_mat,idxs]=NV_070823_ScatvFLRML_6PMT(scatpath,Base_Name,Base_Name2,bead_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  file: NV_070823_ScatvFLRML_6PMT.m
%  authors: Nilay Vora
%  created: 07/08/2023
% modified:
%  purpose: This file will read peak data and correlate it against
%  detection schemes for peaks. Namely, Scat vs. FLR vs. Scat+FLR. This
%  data is then returned as a table to be saved per folder as a .xlsx.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set folder path in which peak values are stored:
cd(scatpath)
%% Set number of chunks of data:
chunk=60;
%% Main Loop
pv1=[];
pv2=[];
idxs = [];
    
load(Base_Name,'-mat','peak_values');
if contains(Base_Name,'Scat')
    if length(peak_values)<1
        ratio_mat = 0; 
        count_mat = 0;
        tot_mat = 0;
        return
    else
        peak_values = peak_values((peak_values(:,20)>=30),:);
        peak_values((peak_values(:,5)<0.15),:) = [];
        if bead_flag == 1
            peak_values(peak_values(:,4)>0.05,:)=[];
        end
    end
else
    peak_values(peak_values(:,10)<17,:)=[];
    if bead_flag == 1
        peak_values(peak_values(:,4)>0.05,:)=[];
    end
end
pv1=[pv1;peak_values];
load(Base_Name2,'-mat','peak_values');
if contains(Base_Name2,'Scat')
    if length(peak_values)<1
        ratio_mat = 0; 
        count_mat = 0;
        tot_mat = 0;
        return
    else
        peak_values = peak_values((peak_values(:,20)>=30),:);
        peak_values((peak_values(:,5)<0.15),:) = [];
        if bead_flag == 1
            peak_values(peak_values(:,4)>0.05,:)=[];
        end
    end
else
    peak_values(peak_values(:,10)<17,:)=[];
    if bead_flag == 1
        peak_values(peak_values(:,4)>0.05,:)=[];
    end
end

pv2=[pv2;peak_values];
% Sort peaks inputted
pv1_sorted=sortrows(pv1,[24,7,8]);
pv2_sorted=sortrows(pv2,[24,7,8]);
count=0;
i_new=[];
er=0;
for j=1:max(unique(pv1_sorted(:,24))) %# of data files
    f1=pv1_sorted(:,24)==j;
    d1=pv1_sorted(f1,:);
    f2=pv2_sorted(:,24)==j;
    d2=pv2_sorted(f2,:);
    for k=1:chunk %Chunk number
        f11=d1(:,7)==k; % grab all peaks from specified chunk
        d11=d1(f11,:);
        f21=d2(:,7)==k; % grab all peaks from specified chunk
        d21=d2(f21,:);
        loc=d11(:,8); % grab all locs from specified chunk
        %%%%% Changed to account for unregisterd clusters in FLR %%%%%
        if contains(Base_Name2,'FLR')
            loc2 = [];
            loc_scat = d21(:,8);  % grab all locs from specified chunk
            width_scat = d21(:,10); % grab all clusters widths from scat channel
            for npeaks = 1:size(d21,1)
                tloc = linspace(loc_scat(npeaks) - width_scat(npeaks),...
                    loc_scat(npeaks) + width_scat(npeaks));
                loc2 = [loc2;tloc'];
            end
        else
            loc2=d21(:,8);  % grab all locs from specified chunk
        end
        
        width=d11(:,10);  % grab all width from first file only
        width(width==0)=81; % if width is 0 we make it insanely large
        for i=1:length(loc) %location in a chunk
            int=loc(i); 
            found=find(loc2>(int-width(i)*1.5) & loc2<(int+width(i)*1.5)); 
            % Search the second file for peaks that fall within 1.5*width
            % of the peak in the first file
            if ~isempty(found) % Is a corresponding peak found
                if length(found)==1
                    s=find(found+er==i_new, 1); % Checks if peak is already
                    % in the saved matrix
                    if isempty(s)
                        count=count+1;
                        i_new=[i_new;found+er]; %#ok<*AGROW> 
                        idxs = [idxs;[j,k,loc2(found)]];  %#ok<*AGROW> 
                        if ~isempty(find(diff(i_new)<0, 1))
                            disp('broke')
                        end
                    end
                else
                    % If more than one peak is found in the range try to
                    % sort out which on we care about
                    dist=[]; 
                    ii_store=[];
                    AUC=[];
                    for ii=1:length(found)
                        s=find(found(ii)+er<=i_new, 1); % Peaks must not be
                        % older than what is already in the matrix
                        if isempty(s)
                            dist=[dist,d11(i,8)-loc2(found(ii))];
                            %AUC=[AUC,d11(i,8)-d21(found(ii),8)];
                            ii_store=[ii_store,ii];
                        end
                    end
                    % We now have the distance and AUC which can be used to
                    % remove the misidentified peak
                    if ~isempty(dist)
                        restricted=[]; %#ok<*NASGU> 
                        test=[];
                        if i<length(loc) % Prevents us from breaking code
                            restricted=loc(i+1); % Grab the next peak loc
                            test=find(restricted>loc2(found(ii_store)));
                            % test checks if the next location is greater
                            % than any of the identified locations
                            if ~isempty(test)
                                ii=ii_store(test);
                                s=[];
                                sfail=[];
                                for idx=1:length(ii)
                                    if found(ii(idx))+er==i_new
                                        sfail=[sfail,find(found(ii(idx))+er==i_new)];
                                    else
                                        s=s; %#ok<*ASGSL> 
                                    end
                                end
                                found(sfail)=[];
                                if ~isempty(dist)
                                    [~,min_ind]=min(abs(dist));
                                    ii=ii_store(min_ind);
                                    s=find(found(ii)+er==i_new, 1);
                                end
                                if isempty(s)
                                    count=count+1;
                                    i_new=[i_new;found(ii)+er];
                                    if contains(Base_Name2,'FLR')
                                        idxs = [idxs;[j,k,d11(i,8)]];  %#ok<*AGROW>
                                    else
                                        idxs = [idxs;[j,k,loc2(found(ii))]];  %#ok<*AGROW>
                                    end
                                    if ~isempty(find(diff(i_new)<0, 1))
                                        disp('broke')
                                    end
                                end
                            end
                        else
                            if ~isempty(dist)
                                [~,min_ind]=min(abs(dist));
                                ii=ii_store(min_ind);
                                s=find(found(ii)+er==i_new, 1);
                                
                                if isempty(s)
                                    count=count+1;
                                    i_new=[i_new;found(ii)+er];
                                    if contains(Base_Name2,'FLR')
                                        idxs = [idxs;[j,k,d11(i,8)]];  %#ok<*AGROW>
                                    else
                                        idxs = [idxs;[j,k,loc2(found(ii))]];  %#ok<*AGROW>
                                    end
                                    if ~isempty(find(diff(i_new)<0, 1))
                                        disp('broke')
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        er=er+length(loc2);
    end
end
%count=length(unique(i_new));
ratio=count/length(pv1_sorted(:,8));
ratio_mat=ratio; 
count_mat=count;
tot_mat=length(pv1_sorted(:,8));
