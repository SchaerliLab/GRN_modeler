function [] = delete_individual(obj)
% delete individual parameters, reactions, rules, species from Mobj

% reactions
delete(sbioselect(obj.Mobj.Reactions,'Notes','Individual'));
% rules
delete(sbioselect(obj.Mobj.Rules,'Notes','Individual'));
% parameters
delete(sbioselect(obj.Mobj.Parameters,'Notes','Individual'));
% species
delete(sbioselect(obj.Mobj.Species,'Notes','Individual'));

end

