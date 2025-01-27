function set(cell,name,PropertyName,PropValue,varargin)
% set the property ('PropertyName') of a given object ('name')
% to a given value ('PropValue')
% from the object listed in varargin
% if the object does not have named properties, it can be empty
% varagin{1}: obj_name
% varagin{2}: reg_name1
% varagin{3}: reg_name2
% e.g. cell.set('k5_N1','Value',5,'N1','R1')
% varargin can be 'cell'. e.g. cell.set('k5_N1','Value',5,{'N1','R1'})

% use the cell element
if iscell(varargin{1}) && numel(varargin)==1
    varargin = varargin{1};
end

% add the 'input' if it is not included
varargin = cell.addinput2route(varargin);

% find the parent object
[parent_obj,good_object] = cell.route2fieldnames(varargin);
if good_object == false
    return;
end

% get the object
parent_obj = eval(parent_obj);

% if it is a property of the parent object directly
if isprop(parent_obj,name)
    parent_obj.(name) = PropValue;
    return;
elseif ~isempty(sbioselect(parent_obj.Mobj.Species,'Name',name))
    % species
    set(sbioselect(parent_obj.Mobj.Species,'Name',name),PropertyName,PropValue);
    return
elseif ~isempty(sbioselect(parent_obj.Mobj.Reactions,'Name',name))
    % reaction
    set(sbioselect(parent_obj.Mobj.Reactions,'Name',name),PropertyName,PropValue);
    return
elseif ~isempty(sbioselect(parent_obj.Mobj.Rules,'Name',name))
    % parameter
    set(sbioselect(parent_obj.Mobj.Rules,'Name',name),PropertyName,PropValue);
    return
elseif ~isempty(sbioselect(parent_obj.Mobj.Parameters,'Name',name))
    % rule
    set(sbioselect(parent_obj.Mobj.Parameters,'Name',name),PropertyName,PropValue);
    return
end

% no property found
warning('We did not find the property of ''%s''.',name);
return;

end