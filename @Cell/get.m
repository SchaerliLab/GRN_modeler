function value = get(cell,name,PropertyName,varargin)
% get the property ('PropertyName') of a given object ('name')
% from the object listed in varargin
% if the object does not have named properties, it can be empty
% varagin{1}: obj_name
% varagin{2}: reg_name1
% varagin{3}: reg_name2
% e.g. cell.get('k5_N1','Value','N1','R1')
% varargin can be 'cell'. e.g. cell.get('k5_N1','Value',{'N1','R1'})

% use the cell element
if iscell(varargin{1}) && numel(varargin)==1
    varargin = varargin{1};
end

% find the parent object
[parent_obj,good_object] = cell.route2fieldnames(varargin);
if good_object == false
    value = [];
    return;
end

% get the object
parent_obj = eval(parent_obj);

% if it is a property of the parent object directly
if isprop(parent_obj,name)
    value = parent_obj.(name);
    return;
elseif ~isempty(sbioselect(parent_obj.Mobj.Species,'Name',name))
    % species
    value = get(sbioselect(parent_obj.Mobj.Species,'Name',name),PropertyName);
    return
elseif ~isempty(sbioselect(parent_obj.Mobj.Reactions,'Name',name))
    % reaction
    value = get(sbioselect(parent_obj.Mobj.Reactions,'Name',name),PropertyName);
    return
elseif ~isempty(sbioselect(parent_obj.Mobj.Rules,'Name',name))
    % parameter
    value = get(sbioselect(parent_obj.Mobj.Rules,'Name',name),PropertyName);
    return
elseif ~isempty(sbioselect(parent_obj.Mobj.Parameters,'Name',name))
    % rule
    value = get(sbioselect(parent_obj.Mobj.Parameters,'Name',name),PropertyName);
    return
end

% no property found
warning('We did not find the property of ''%s''.',name);
value = [];
return;

end