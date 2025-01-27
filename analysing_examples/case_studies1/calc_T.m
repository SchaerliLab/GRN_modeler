function [T] = calc_T(c,t,T_approx,ratio)

% we need at least two data points
% the relative amplitude should be at least ratio
if(numel(t)>2 && (max(c)-min(c))/max(c) > ratio)
    T = calc_timeperiod_acf(c,t,T_approx);
else
    T = nan;
end

end

