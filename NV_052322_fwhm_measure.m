function[fwhm,peak_area]=NV_052322_fwhm_measure(data_fwhm,peak_height)
% Takes data and measures fwhm, and peak area for all channels
%% 052322 Update
% Code is now set to normalize peaks within a range for more accurate FWHM
% measurement as previously, low and high peaks caused undermeasurement and
% poor labeling.
%% Grab FWHMs
% step left
%[~,locs2]=findpeaks(data_fwhm,'NPeaks',3,'SortStr','descend');
[~,locs2]=findpeaks(data_fwhm,'SortStr','descend');
locs2=sort(locs2,'ascend');
peaks=data_fwhm(locs2);
locs2(peaks<0.15.*peak_height)=[];
peaks(peaks<0.15.*peak_height)=[];
data2=data_fwhm;
if length(locs2)>1
    data2(1:locs2(1))=data2(1:locs2(1))./peaks(1);
    m=(peaks(2)-peaks(1))./(locs2(2)-locs2(1));
    b=peaks(2)-locs2(2).*m;
    y=m.*((locs2(1)+1):(locs2(2)-1))+b;
    data2((locs2(1)+1:(locs2(2)-1)))=data2((locs2(1)+1:(locs2(2)-1)))./y';
    if length(locs2)>2
        data2(locs2(2))=1;
        m=(peaks(3)-peaks(2))./(locs2(3)-locs2(2));
        b=peaks(3)-locs2(3).*m;
        y=m.*((locs2(2)+1):(locs2(3)-1))+b;
        data2((locs2(2)+1:(locs2(3)-1)))=data2((locs2(2)+1:(locs2(3)-1)))./y';
        data2(locs2(3):length(data2))=data2(locs2(3):length(data2))./peaks(3);
    else
        data2(locs2(2):length(data2))=data2(locs2(2):length(data2))./peaks(2);
    end
    %data2=(data2-min(data2))./(max(data2)-min(data2));
else
    data2=data2./max(data2);
    %data2=(data2-min(data2))./(max(data2)-min(data2));
end
current_point=peak_height;
k=find(data_fwhm==current_point);
if length(k)>1
    k=min(k);
end
k1=k;
%amp=peak_height;
%half_max=0.5*amp;
half_max=0.5;
left_dist=0;
peak_area_left=0;

startpt=find(data2(1:end-1)<=half_max & data2(2:end) >= half_max,1,'first');
endpt=find(data2(1:end-1)>=half_max & data2(2:end) <= half_max,1,'last');
if isempty(startpt)
    startpt=1;
elseif isempty(endpt)
    endpt=length(data_fwhm);
end

while k1 >= startpt
    current_point=data_fwhm(k1);
    left_dist=left_dist+1;
    peak_area_left=abs(current_point)+peak_area_left;
    k1=k1-1;
end
% step right
current_point=peak_height; %#ok<NASGU> 
k1=k;
right_dist=0;
peak_area_right=0;
while k1 <= endpt
    current_point=data_fwhm(k1);
    right_dist=right_dist+1;
    peak_area_right=abs(current_point)+peak_area_right;
    k1=k1+1;
end
% Adding left and right contributions
fwhm=left_dist+right_dist;
peak_area=peak_area_left+peak_area_right-peak_height;
end