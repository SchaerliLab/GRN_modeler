function [] = delete(obj)
% delete the model

delete(obj.Reactions);
delete(obj.Rules);
delete(obj.Parameters);
delete(obj.Species);

end
