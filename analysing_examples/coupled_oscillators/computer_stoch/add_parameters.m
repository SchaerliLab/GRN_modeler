function [Mobj] = add_parameters(Mobj,p)
% Set parameters in a simbiology model object
% Mobj: simbiology model object
% c0: struct, the parameters will be set according to the 
% names and values of p 

names = fieldnames(p);
for i = 1:length(names)
    addparameter(Mobj,names{i},p.(names{i}));
end

