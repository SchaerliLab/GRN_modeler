function obj = clean_up(obj)
% delete the model

obj.Reactions(~isvalid(obj.Reactions)) = [];
obj.Rules(~isvalid(obj.Rules)) = [];
obj.Parameters(~isvalid(obj.Parameters)) = [];
obj.Species(~isvalid(obj.Species)) = [];

end
