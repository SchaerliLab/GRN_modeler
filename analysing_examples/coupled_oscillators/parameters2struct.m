function [p] = parameters2struct(Mobj)
% Create a structure from the parameters in a simbiology model object
% Mobj: simbiology model object

names = get(Mobj.Parameters,'Name');
for i = 1:length(names)
    if ~isvarname(Mobj.Parameters(i).Name)
        % rename parameters parameter_name => param1, ..., paramn
        p.(['param' int2str(i)]) = get(sbioselect(Mobj.Parameters,'Name',names{i}),'Value');
    else
        p.(names{i}) = get(sbioselect(Mobj.Parameters,'Name',names{i}),'Value');
    end
end

