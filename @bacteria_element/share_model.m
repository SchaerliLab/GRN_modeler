function obj = share_model(obj,Mobj)
% add model (Mobj) to the shared common model
% make a (handle class) copy for the properties in the actual
% object (obj.Mobj)

% species
obj = copy_property(obj,Mobj,obj.data.Mobj,'Species');
% parmateres
obj = copy_property(obj,Mobj,obj.data.Mobj,'Parameters');
% rules
obj = copy_property(obj,Mobj,obj.data.Mobj,'Rules');
% reactions
obj = copy_property(obj,Mobj,obj.data.Mobj,'Reactions');

% if the UserData is empty, we will copy it
if isempty(fieldnames(obj.Mobj.UserData))
    obj.Mobj.UserData = Mobj.UserData;
end

end

function obj = copy_property(obj,Mobj_in,Mobj_out,property)
% copy a given property into the shared model (if it is not there yet)
% make a handle for the property in the object (Mobj_common)

for i = 1:numel(Mobj_in.(property))
    % make a copy in the shared model if it is not there yet
    handle = sbioselect(Mobj_out.(property),'Name',Mobj_in.(property)(i).Name);
    % we copy it if it has a new name
    if isempty(handle)
        handle = copyobj(Mobj_in.(property)(i),Mobj_out);
    end
    % make a handle copy from the shared model in the actual object
    obj.Mobj.(property)(i) = handle;
end

end