function T = calc_timeperiod_avr_cross(t,c)
% calculate the time period from the crossing of the average value
% t: time
% c: the variable
% it works only when the average is not increasing or decreasing

c_mean = mean(c);

% if the reltive change in the signal is small
% then it is considered just as noise
if (max(c)-min(c))/(c_mean) < 0.001/100
    T = nan;
    return
end

% find the points which crosses the mid value with positiove slope
midpoints = c(2:end)>c_mean & c(1:end-1)<=c_mean;
midpoints(end+1) = false;
tmidpoints = t(midpoints);

% calulate the avrage time between the crosses
if ~isempty(tmidpoints)
    T = (tmidpoints(end)-tmidpoints(1))/(length(tmidpoints)-1);
else
    T = nan;
end

end

