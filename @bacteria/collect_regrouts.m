function [regrouts] = collect_regrouts(obj)
% calculate every regulator route signed e.g.:
% {{'N2','HILL','|-','P_N1'}, {'N3','HILL3','|-','P_N2'}}

regrouts = struct2cell(obj.data.regulator_routes);
regrouts = [regrouts{:}];

end

