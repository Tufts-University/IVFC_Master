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
locs2(peaks<0.05.*peak_height)=[];
peaks(peaks<0.05.*peak_height)=[];
data2=data_fwhm;
init = 1;
if length(locs2)>1
    for n = 1:length(peaks)-1
        data2(init:locs2(n)) = data2(init:locs2(n))./peaks(n);
        m = (peaks(n+1)-peaks(n))./(locs2(n+1)-locs2(n));
        b =  peaks(n+1)-locs2(n+1).*m;
        y= m.*((locs2(n)+1):(locs2(n+1)-1))+b;
        data2((locs2(n)+1:(locs2(n+1)-1)))=data2((locs2(n)+1:(locs2(n+1)-1)))./y';
        init = locs2(n+1);
    end
    data2(init:end) = data2(init:end)./peaks(end);
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