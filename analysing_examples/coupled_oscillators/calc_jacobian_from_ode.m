function [f_jac] = calc_jacobian_from_ode(f,p,varargin)
% calculate Jacobian analitically from ODE
% the variables in f should be named as c(1),..,c(n)
% f: function handle for anonymous function
% p: struct for the parameters in f. (eg: p.k = 4, p.alfa = 3.4)
% f_jac: function handle for the jacobian matrix
% varargin: names of the variables which will included as input for the
% Jacobian

str = func2str(f);

% number of variables
c_var_pos = strfind(str,'c(');
nvar = 0;
for i = 1:length(c_var_pos)
    tmp = sscanf(str(c_var_pos(i)+2:end),'%d');
    if tmp > nvar
        nvar = tmp;
    end
end

% rename species eg: c(1) => c1
vars = cell(nvar,1);
for i = 1:nvar
    vars{i} = ['c',int2str(i)];
    str = strrep(str,['c(',int2str(i),')'],vars{i});
end

% add given values to the parameters
if nargin > 1
    var_names = fieldnames(p);
    for i = 1:length(var_names)
        eval([var_names{i},'= p.(var_names{i});'])
    end
end

% calulate the Jacobian analitically
J = jacobian(str2sym(str),str2sym(vars));

% make anonymous function from the symbolic function
F = matlabFunction(J);
str = func2str(F);
% delete original argument
str(1:find(str==')',1,'first')) = [];
% rename species eg: c1 => c(1)
for i = 1:nvar
    vars{i} = ['c',int2str(i)];
    str = strrep_smart(str,vars{i},['c(',int2str(i),')']);
end

% extra output names
output_vars = [];
if nargin > 2
    if iscell(varargin{1})
        varargin = varargin{1};
    end
    varargin = strcat(',',varargin);
    output_vars = [varargin{:}];
end
% new argument
str = ['@(c' output_vars ')',str,';'];
f_jac = eval(str);
end
