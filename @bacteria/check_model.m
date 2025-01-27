function [] = check_model(~,Mobj)
% every property of the simbio model should have 'Notes' with either
% 'Individual' or 'Common' value

% Species
check_notes(Mobj,'Species');
check_units(Mobj,'Species');
% check_name(Mobj,'Species');

% Parameters
check_notes(Mobj,'Parameters');
check_units(Mobj,'Parameters');
% % no '_' if it is a 'HILL' parameter
% for i = 1:numel(Mobj.Parameters)
%     if any(contains(Mobj.Parameters(i).Name,'HILL')) && any(contains(Mobj.Parameters(i).Name,'_'))
%         warning('%s ''%s'' name contains ''_''!',prop,Mobj.(prop)(i).Name)
%     end
% end

% Reactions
check_notes(Mobj,'Reactions');

% Rules
check_notes(Mobj,'Rules');

end

function check_notes(Mobj,prop)
% 'Notes' should be 'Individual' or 'Common' for a given 'prop' property
% (e.g. 'Species', 'Parameters') for the simbiology model object (Mobj).
for i = 1:numel(Mobj.(prop))
    if ~(strcmp(Mobj.(prop)(i).Notes,'Individual') || strcmp(Mobj.(prop)(i).Notes,'Common'))
        warning('%s ''%s'' does not have proper ''Notes''!',prop,Mobj.(prop)(i).Name)
    end
end
end

function check_units(Mobj,prop)
% 'Units' should not be empty for a given 'prop' property
% (e.g. 'Species', 'Parameters') for the simbiology model object (Mobj).
for i = 1:numel(Mobj.(prop))
    if isempty(Mobj.(prop)(i).Units)
        warning('%s ''%s'' does not have ''Units''!',prop,Mobj.(prop)(i).Name)
    end
end
end

% function check_name(Mobj,prop)
% % 'Species' 'Name' should not contain '_' if its 'Tag' is 'Regulator'
% for i = 1:numel(Mobj.(prop))
%     if strcmp(Mobj.(prop)(i).Tag,'Regulator') && any(contains(Mobj.(prop)(i).Name,'_'))
%         warning('%s ''%s'' name contains ''_''!',prop,Mobj.(prop)(i).Name)
%     end
% end
% end

