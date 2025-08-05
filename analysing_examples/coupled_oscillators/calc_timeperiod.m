function [T,err,n] = calc_timeperiod(t,c,mindist)
% calculate the time period
% t: time
% c: the variable
% T: time period
% n: number of peaks

% find the peaks
if nargin >= 3
    [~,locs] = findpeaks(c,'MinPeakHeight',max(c)/2,'MinPeakDistance',mindist,'MinPeakProminence',mean(c));
else
    [~,locs] = findpeaks(c,'MinPeakHeight',mean(c,'omitnan'));
end

% distances between the peaks
Ti = diff(t(locs));
% calculate the average time period
T = mean(Ti);
% error
err = std(Ti);
% number of peaks
n = numel(locs);

end

