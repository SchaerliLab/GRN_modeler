function [cell,suffix] = add_regulator(cell,type,varargin)
% add regulator to an object or its regulators
% obj |- reg1 |- reg2 |- ...
% cell: Cell object
% for varargin we can have two types: it contains only the object names and
% in this case the first input will be used:
% varagin{1}: obj_name
% varagin{2}: reg_name1
% varagin{3}: reg_name2
% or it can contain the input names after the object name:
% varagin{1}: obj_name
% varagin{2}: obj_name_input
% varagin{3}: reg_name1
% varagin{4}: reg_name1_input
% varagin{5}: reg_name2

% if we are using cell input
if iscell(varargin{1})
    varargin = varargin{1};
end

% add a name for the regulator if it is not provided
if numel(varargin)==1
    warning('A regulator name is necessary.')
    return;
end

% original regulator name
original_name = cell.original_node_name(varargin{end});

% if the regulator is produced by a node, and just the node name is given
% we will add a prefix (the first output) except for the first one,
% where the node might be regulated directly
for i = 2:numel(varargin)
    if isfield(cell.nodes,varargin{i})
        % we are using the first prefix if nothing else is given
        varargin{i} = [cell.nodes.(varargin{i}).Mobj.UserData.output_prefix{1},varargin{i}];
    end
end

% if the 'input' names are not specified in varargin
% lets put the first input names in it
varargin = cell.addinput2route(varargin);

% find the regulated object
[regulated_obj,good_object] = cell.route2fieldnames(varargin(1:end-2));

% suffix for the parameter names, e.g. [{'N1'} {'HILL1'} {'|-'} {'R1'} {'HILL3'} {'<-'} {'R2'}]
if numel(varargin) > 3
    suffix = eval([regulated_obj '.suffix']);
    suffix{end+1} = varargin{end-1};
else
    suffix(1,1:2) = varargin(1:2);
end
% extend the suffix according to the regulation type
if contains(type,'R')
    suffix{end+1} = '|-';
else
    suffix{end+1} = '<-';
end
suffix{end+1} = varargin{end};

% get the simplified suffix of the regulated object
if numel(varargin) > 3
    % this is a regulator, so it has a suffix_simple
    suffix_simple = eval([regulated_obj '.suffix_simple']);
else
    suffix_simple = varargin(1);
end
% if we have just one input
if numel(eval([regulated_obj '.Mobj.UserData.input'])) < 2
    suffix_simple{1,end+1} = '';
else
    suffix_simple{1,end+1} = suffix{end-2}; % eg. 'HILL'
end
suffix_simple{end+1} = suffix{end-1}; % eg. '<-'
suffix_simple{end+1} = original_name;
% if we have only one input or output, we will not use them in the route
% suffix_simple = cell.simplify_suffix(suffix_simple);

% add regulator to the regulator list if it is not there yet
if ~isfield(cell.data.regulator_routes,varargin{end})
    cell.data.regulator_routes(1).(varargin{end}) = {};
end

% add the route to the regulator list (useful when we are deleting
% regulators), let it be a row vector
cell.data.regulator_routes.(varargin{end}){end+1} = suffix;% reshape(varargin,1,numel(varargin));

% % if Hill_target is empty, the first input ('HILL' parameter) will be used
% % for the regulated object
% if isempty(Hill_target)
%     Hill_target = eval([regulated_obj '.Mobj.UserData.input{1}']);
% end
if good_object==true % if the object exist
    % we will add the regulator if it is not added yet
    if ~any(strcmp(fieldnames(eval([regulated_obj '.regulators'])),varargin{end-1})) ||...
            ~any(strcmp(fieldnames(eval([regulated_obj '.regulators.(''' varargin{end-1} ''')'])),varargin{end}))
        regulator_name = varargin{end};
        eval([regulated_obj ' = ' regulated_obj '.add_regulator(''' regulator_name ''',''' type ''',suffix,suffix_simple,' regulated_obj ');']);
    end
end
end