function [T] = calc_T2(c,t,ratio)

% we need at least two data points
% the relative amplitude should be at least ratio
if(numel(t)>2 && (max(c)-min(c))/max(c) > ratio)
    T = calc_timeperiod_avr_cross(t,c);
else
    T = nan;
end

end

