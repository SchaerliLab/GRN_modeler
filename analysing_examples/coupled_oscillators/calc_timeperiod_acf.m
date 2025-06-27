function [T] = calc_timeperiod_acf(c,t,T_approx)
%CALC_TIMEPERIOD_ACF calculate time period based on autocorrelation
%function
% t: time, should be equidistant, increasing
% T_approx: approximate time period to determine the necessary time lags
% T: calulated time period, nan, if we did not find a peak in the acf

dt = diff(t(1:2));

if nargin < 3
    % no T_approx, use everything
    n_lags = length(c)-1;
else
    % estimated time lags
    n_lags = min(round(3*T_approx/dt),length(c)-1);
end

acf = autocorr(c,n_lags);

[~,lks] = findpeaks(acf);

if isempty(lks)
    T = nan;
else
    T = lks(1)*dt;
end

end

