function [duty,dt_high,dt_low] = calc_duty(t,c,thres)
%CALC_DUTYcalculate duty cycle in periodic signal
% t: time
% c: the signal (concentration)
% thres: the threshold for the signal

% shift the signal to see where we cross the treshold
c = c-thres;
% we change sign where we cross
cross = c(1:end-1).*c(2:end);
% the time points for the threshold crossing
t_cross = t(cross<0);
% time spent in the different states
dt = diff(t_cross);
% high and low states
if c(find(cross<0,1,"first")+1)>0 % if we start with high states
    dt_high = mean(dt(1:2:end));
    dt_low = mean(dt(2:2:end));
else
    dt_high = mean(dt(2:2:end));
    dt_low = mean(dt(1:2:end));
end

% calc the duty cycle
duty = dt_high/(dt_high+dt_low);

end

