function[fwhm,peak_area]=NV_101719_fwhm_measure(data_fwhm,peak_height)
% Called by Peak_Detection_07142015
% Takes data and measures fwhm, and peak area for all channels
  %% Grab FWHMs
        % step left
        %% New Code (10/16/19) N.V.
%         if peak_height==max(data_fwhm)
%             current_point=peak_height;
%         else
%             peak_height=max(data_fwhm);
%             current_point=peak_height;
%         end
       
        current_point=peak_height;
        k=find(data_fwhm==current_point);
        if length(k)>1
            k=min(k);
        end
        k1=k;
        amp=peak_height;
        half_max=0.5*amp;
        left_dist=0;
        peak_area_left=0;
        
        startpt=find(data_fwhm(1:end-1)<=half_max & data_fwhm(2:end) >= half_max,1,'first');
        endpt=find(data_fwhm(1:end-1)>=half_max & data_fwhm(2:end) <= half_max,1,'last');
        while k1 >= startpt
            current_point=data_fwhm(k1);
            left_dist=left_dist+1;
            peak_area_left=abs(current_point)+peak_area_left;
            k1=k1-1;
        end
        % step right
        current_point=peak_height;
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