function [] = delete_individual(obj)
% delete the model

delete(obj.Reactions(strcmp({obj.Reactions.Notes},'Individual')));
delete(obj.Rules(strcmp({obj.Rules.Notes},'Individual')));
delete(obj.Parameters(strcmp({obj.Parameters.Notes},'Individual')));
delete(obj.Species(strcmp({obj.Species.Notes},'Individual')));

end
