function [] = add_model(~,Mobj_in,Mobj_out,notes)
% add a simbiology model Mobj_in (specieses, parameters, reactions, rules) to
% another simbio model Mobj_out (Mobj_out += Mobj_in)
% if notes is given, the object (species parameter etc) will be added only
% with this specific notes (e.g. 'individual' or 'common')

if isempty(Mobj_in)
    % there is nothing to add
    return;
end

if nargin < 4
    % we did not prescribe the notes
    notes = [];
end

% species
copy_property(Mobj_in,Mobj_out,'Species',notes)
% parmateres
copy_property(Mobj_in,Mobj_out,'Parameters',notes)
% rules
copy_property(Mobj_in,Mobj_out,'Rules',notes)

% reactions
if nargin >=4 % specific notes are prescribed
    % reactions
    for i = 1:numel(Mobj_in.Reactions)
        % reactions should not have the same name if it is not an automatic
        % 'Reaction_XY' type name
        if strcmp(Mobj_in.Reactions(i).Notes,notes) && (contains(Mobj_in.Reactions(i).Name,'Reaction') || isempty(sbioselect(Mobj_out.Reactions,'Name',Mobj_in.Reactions(i).Name)))
            copyobj(Mobj_in.Reactions(i),Mobj_out);
        end
    end

else
    % reactions
    for i = 1:numel(Mobj_in.Reactions)
        % reactions should not have the same name if it is not an automatic
        % 'Reaction_XY' type name
        if  (contains(Mobj_in.Reactions(i).Name,'Reaction') || isempty(sbioselect(Mobj_out.Reactions,'Name',Mobj_in.Reactions(i).Name)))
            copyobj(Mobj_in.Reactions(i),Mobj_out);
        end
    end
end

% if the UserData is empty, we will copy it
if isempty(Mobj_out.UserData)
    Mobj_out.UserData = Mobj_in.UserData;
end

end

function copy_property(Mobj_in,Mobj_out,property,notes)
% copy a given property fromMobj_in to Mobj_out
if ~isempty(notes) % specific notes are prescribed
    for i = 1:numel(Mobj_in.(property))
        if strcmp(Mobj_in.(property)(i).Notes,notes) && isempty(sbioselect(Mobj_out.Species,'Name',Mobj_in.(property)(i).Name))
            copyobj(Mobj_in.(property)(i),Mobj_out);
        end
    end
else
    for i = 1:numel(Mobj_in.(property))
        if isempty(sbioselect(Mobj_out.(property),'Name',Mobj_in.(property)(i).Name))
            copyobj(Mobj_in.(property)(i),Mobj_out);
        end
    end
end
end
