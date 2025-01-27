function [str,success] = route2fieldnames(obj,varargin)
% create fieldnames for the given route to a certain object
% e.g. route2fieldnames('N1','HILL1','R1','HILL4','R2') 
% or route2fieldnames({'N1','HILL1','R1','HILL4','R2'})
% if the input is not given, we will use the first possible input
% => cell.nodes('N1').regulators.('HILL1').('R1').regulators.('HILL4').('R2')
% success: check whether the field exsists (true/false)

% if we are using cell input
if iscell(varargin{1})
    varargin = varargin{1};
end

% if it was called with the regulation type, get rid of that
% e.g.: {'G','HILL','|-','P_R3'} => {'G','HILL','P_R3'}
varargin(strcmp(varargin,'<-')|strcmp(varargin,'|-')) = [];

% find the first regulated object (node/protease)
if isfield(obj.nodes,varargin{1})
    str = ['cell.nodes.(''' varargin{1} ''')'];
elseif isfield(obj.proteases,varargin{1})
    str = ['cell.proteases.(''' varargin{1} ''')'];
else
    if iscell(varargin{1})
        error_name = strcat('_',varargin{1});
        error_name = [error_name{:}];
        error_name = error_name(2:end);
    else
        error_name = varargin{1};
    end
    error('We did not find the object ''%s''.',error_name)
end

% % add node prefix
% for i = 3:2:numel(varargin)
%     if isfield(obj.nodes,varargin{i})
%         % if the regulator is produced by a node, and just the node name is given
%         % we will add prefix
%         varargin{i} = [obj.nodes.(varargin{i}).Mobj.UserData.output_prefix,varargin{i}];
%     end
% end

% we have further regulators
if numel(varargin) >= 2
    names = strcat('.regulators.(''',varargin(2:2:end),''').(''',varargin(3:2:end),''')');
    names = [names{:}];
    str = [str names];
end

% check whether the field exist
success = checkstruct(obj,str);

end