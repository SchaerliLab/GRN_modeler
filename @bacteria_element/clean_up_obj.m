function obj = clean_up_obj(obj)
% clear the deleted properties from the model

obj.Mobj = obj.Mobj.clean_up();

end

