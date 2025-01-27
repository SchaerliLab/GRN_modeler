function [] = rename_prop(obj,Mobj,suffix,exception)
% rename individual properties in simbiology model Mobj
% adding 'suffix' to the end of their name
% except if their name is the same as the object name or it is in exception

if nargin < 4
    exception = [];
end

if nargin < 3
    suffix = ['_',obj.name];
end

individ_pars = sbioselect(Mobj,'Notes','Individual');
for i = 1:numel(individ_pars)
    % rename if it does not have the name of the object
    if ~strcmp(individ_pars(i).Name,obj.name) && ~any(strcmp(individ_pars(i).Name,exception))
        new_name = [individ_pars(i).Name,suffix];
        rename(individ_pars(i),new_name);
    end
end

% % if we have input in the UserData, rename the species there
% if isa(Mobj.UserData,'struct') && isfield(Mobj.UserData,'input')
%     for i = 1:numel(Mobj.UserData.input)
%         % rename it if it is not in the exception
%         if ~strcmp(Mobj.UserData.input(i),obj.name) && ~any(strcmp(Mobj.UserData.input(i),exception))
%             Mobj.UserData.input{i} = [Mobj.UserData.input{i},suffix];
%         end
%     end
% end

end

