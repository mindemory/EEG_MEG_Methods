function [xnew, spkindx] = despike_dilate(x, spk_thresh, exclusion_range, do_log, spk_side, permission_range)

% function xclean = despike(x, spk_thresh)
%
% Identifies spikes from time-series data;
% removes data and replaces using cubic interpolation.
%
% INPUT
%
% x = [Nsamples by Nsignals] matrix of data, containing spikes
% spk_thresh = threshold for whether a data point is a spike (expressed as multiples of inter-quartile range)
% exclusion_range  = removes N time points before and after the spike (0 = spike only)
% do_log = take logarithm of data (to normalize signals with long-tails), find spikes in log-data, but interpolate in real data 
% spk_side = flag for spike direction; 
%               0 --> positive and naegative spikes (default)
%               -1 --> negative spikes only
%               +1 --> positive spikes only
%
% OUTPUT
%
% xnew = data with spike removed, and data interpolated in place
% spkindx = index to locations at which input data was altered
%
% version: 0.1
% author: C. Honey, March 2011

% changed: s-michelmann: included permission range fixed
% dilate_binary_vector changed cubic interp to pchip 
%


if nargin < 6; permission_range = ones(size(x,1),1); end
if nargin < 5; spk_side = 0; end
if nargin < 4; do_log = false; end
if nargin < 3; exclusion_range = 0; end
if nargin < 2; spk_thresh = 2.5; end

[Nsamp, Nsig] = size(x);

xnew = x;

%now we find outliers, based on deviations from the median

if do_log
   y = log(x);  %log-transform data before finding outliers
   meds = median(y);
   iqs = iqr(y);
   if spk_side
       spk = bsxfun(@gt, spk_side*bsxfun(@minus, y, meds), spk_thresh * iqs);
   else
       spk = bsxfun(@gt, abs(bsxfun(@minus, y, meds)), spk_thresh * iqs);
   end
else
    meds = median(x);
    iqs = iqr(x);
    if spk_side
       spk = bsxfun(@gt, spk_side*bsxfun(@minus, x, meds), spk_thresh * iqs);
   else
       spk = bsxfun(@gt, abs(bsxfun(@minus, x, meds)), spk_thresh * iqs);
   end
end

spk = bsxfun(@and,spk, permission_range);

spkindx = find(spk);

if spkindx   % if we find any spikes
    for isig = 1:Nsig
        s = spk(:,isig);  %indices of spikes in this signal)
        if sum(s) > 0;  %are there any spikes?
            if exclusion_range > 0  %do we want to exclude data around the spikes?
                s = dilate_binary_vector(s, exclusion_range);  %expand the data regions to be interpolated
            end
                
            tgood = find(~s);
            
            tbad = find(s);
%             keyboard
            xnew(tbad,isig) = interp1(tgood, x(tgood,isig), tbad, 'pchip');  %interpolate over data segments containing spikes
        end
    end
end
end


function w = dilate_binary_vector(v, N)   %dilates the blocks of ones in a binary vector
    b = ones(N+1,1);
    w = (filtfilt(b,1,v)>0);
 
    %w = v | [zeros(N,1); v(1:end-N)] | [v(N+1:end); zeros(N,1)];
end



