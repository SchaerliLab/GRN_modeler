function obj = extend_model(obj,Mobj)
% extend the model with another Mobj_class object

obj = extend(obj,Mobj,'Species');
obj = extend(obj,Mobj,'Reactions');
obj = extend(obj,Mobj,'Rules');
obj = extend(obj,Mobj,'Parameters');

% extend the userdata
fields = fieldnames(Mobj.UserData);
for i = 1:numel(fields)
    if isempty(obj.UserData) || ~isfield(obj.UserData,fields{i})
        % add the whole field
        obj.UserData(1).(fields{i}) = Mobj.UserData.(fields{i});
    else
        % check the elements of the field
        data = Mobj.UserData.(fields{i});
        for j = 1:numel(data)
            % add the missing elements
            if ~any(contains([obj.UserData.(fields{i})],data{j}))
                obj.UserData.(fields{i}){end+1} = data{j};
            end
        end
    end
end

end

function obj = extend(obj,Mobj,property)
% extend the given property if it has a new name
for i = 1:numel(Mobj.(property))
    if isempty(sbioselect(obj.(property),'Name',Mobj.(property)(i).Name))
        obj.(property)(end+1) = Mobj.(property)(i);
    end
end
end
