function [t,c,names,t_cross,pos_min] = crossing_times(Mobj,p)
%CROSSING_TIMES calculate the time points in a simulation when a
%species concentration crosses a value from positive direction to calculate
%poincare section and map
% Mobj: simbiology model object
% p.t_min: transient time
% p.c_target: target concentartion for Poincare section
% p.target_name: the name of the species which should cross a certain value
% p.n_interp: in interpolation use n_interp in both directions

% simulate
simdata = sbiosimulate(Mobj);
[t,c,names] = getdata(simdata);

% find the points where uP_N1 crosses c_target from the positive direction
y = c(:,strcmp(names,p.target_name)) - p.c_target;
dc = y(1:end-1).*y(2:end);
% find where dc is negative (the target species crosses c_target)
pos = 1:length(dc);
pos = pos(dc<0);
% get rid of the cases when uP_N1 started under c_target
pos(y(pos)<0) = [];
% we are not using the time points before t_min
pos_min = find(t>p.t_min,1,'first');
pos(pos<pos_min) = [];

% we need some extra data point for the interpolation
% we get rid of the last crossing if it happens too close to the end of the
% simulation
if length(t) < pos(end)+p.n_interp
    pos(end) = [];
end

% interpolate to find the appropriate time values
t_cross = zeros(length(pos),1);
for i = 1:numel(pos)
    range = pos(i)-p.n_interp:pos(i)+p.n_interp;
    t_cross(i) = interp1(y(range),t(range),0);
end

end

