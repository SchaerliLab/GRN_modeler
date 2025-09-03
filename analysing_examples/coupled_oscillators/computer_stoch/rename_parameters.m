function Mobj = rename_parameters(Mobj)
%RENAME_PARAMETERS rename parameters to get proper C++ format
% Mobj: SimBiology model object

% rename the parameters to store them in "p." struct later
% and replace not C++ compatible names
for i = 1:numel(Mobj.Parameters)
    % change '<-' to '_A_' and '|-' to '_R_'
    parname = Mobj.Parameter(i).Name;
    parname = strrep(parname, '<-', '_A_');
    parname = strrep(parname, '|-', '_R_');
    rename(Mobj.Parameter(i),['p.' parname]);
end

end

