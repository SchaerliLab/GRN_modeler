function [f,fv] = simbio2ode(Mobj_in,argnames)
% convert simbiology model to ODE function
% Mobj: simbiology model
% argnames: cell for the parameters names what we would 
% like to see in the argument of the output ODE. eg: argnames = {'A';'IPTH'};
% f = f(c,argnames1,argnames2,...)
% fv = fv(c1,c2,...,argnames1,argnames2,...)

Mobj = copyobj(Mobj_in);

% rename species c1 => cn
names_original = get(Mobj.Species,'Name');
for i = 1:length(names_original)
    rename(sbioselect(Mobj.Species,'Name',names_original{i}),['c(' int2str(i) ')']);
end

M = getstoichmatrix(Mobj);

% isconst = cell2mat({Mobj.Species.Constant});
% M(isconst,:) = [];
% nonconstvars = (isconst==false);
% names = {Mobj.Species(nonconstvars).Name};
% names_const = {Mobj.Species(isconst).Name};
% nvar = size(names,2);
% nr = size(M,2);

% % position of the constant and non constant species
% tmp = 1:numel(Mobj.Species);
% pos_const = tmp(isconst);
% pos_nonconst = tmp(nonconstvars);

% rename parameters parameter_name => param1, ..., paramn
% if it is not an accepted variable name
for i = 1:length(Mobj.Parameters)
    if ~isvarname(Mobj.Parameters(i).Name)
        rename(Mobj.Parameters(i),['param' int2str(i)]);
    end
end

% set the parameters
for i = 1:size(Mobj.Parameters,1)
    name = get(Mobj.Parameters(i),'Name');
    value = get(Mobj.Parameters(i),'Value');
    eval([name ' = value;'])
end
% % set the constant species with a new name
% for i = 1:numel(pos_const)
%     name = ['c' int2str(pos_const(i))];
%     value = get(Mobj.Species(pos_const(i)),'Value');
%     eval([name ' = value;'])
% end

% % set constant variables
% names_const = {Mobj.Species(isconst).Name};
% for i = 1:size(names_const,2)
%     name = names_const{i};
%     value = get(sbioselect(Mobj.Species,'Name',name),'Value');
%     eval([name ' = value;'])
% end

% % get the reactions change the specieses
% reactions = cell(nr,1);
% for i = 1:nr
%     reactions{i} = get(Mobj.Reactions(i),'ReactionRate');
%     for j = 1:nvar
%         reactions{i} = strrep_smart(reactions{i},names{j},['c(',int2str(pos_nonconst(j)) ')']);
%     end
%     % % change the constans specieses as well
%     % for j = 1:numel(pos_const)
%     %     reactions{i} = strrep_smart(reactions{i},names_const{j},['c',int2str(pos_const(j))]);
%     % end
% end

% extra argument for the constant variables
arg_extra = '';
if nargin == 2
    % argnames should be a cell array
    if ischar(argnames)
        argnames = {argnames};
    end
    for i = 1:size(argnames,1)
        arg_extra = [arg_extra,',',argnames{i}];
    end
end

% creating the function
f = [];
fv_str = ['f = @(c' arg_extra ') ['];
for i = 1:numel(Mobj.Species) %nvar
    % if we have any reactions
    if any(M(i,:))
        for j = 1:size(M,2)
            if M(i,j) ~= 0
                % fv_str = [fv_str,num2str(M(i,j),'%+d'),'*(',reactions{j},')'];
                fv_str = [fv_str,num2str(M(i,j),'%+d'),'*(',get(Mobj.Reactions(j),'ReactionRate'),')'];
            end
        end
    else
        % constant variable
        fv_str = [fv_str,'0'];
    end
    fv_str = [fv_str,';'];
end
fv_str = [fv_str,'];'];
% replace rules in the function
for i = 1:size(Mobj.Rules,1)
    str = Mobj.Rules(i).Rule;
    newstr = strtrim(split(str,'='));
    % % ci => c(i)
    % for j = 1:nvar
    %     newstr{2} = strrep_smart(newstr{2},['c' int2str(j)],['c(',int2str(j) ')']);
    % end
    newstr{2} = ['(' newstr{2} ')'];
    fv_str = strrep_smart(fv_str,newstr{1},newstr{2});
end
% create the function
eval(fv_str)

if nargout == 2
    % for vectorization
    fv_str = strrep(fv_str,'*','.*');
    fv_str = strrep(fv_str,'/','./');
    fv_str = strrep(fv_str,'^','.^');
    % c(i) <-> ci
    arg = '';
    for i = 1:numel(Mobj.Species) %nvar
        arg = [arg,'c',int2str(i)];
        if i~=numel(Mobj.Species) %nvar
            arg = [arg,','];
        end
    end
    pos = find(fv_str=='c',1,'first');
    for i = 1:numel(Mobj.Species)
        fv_str = strrep(fv_str,['c(' int2str(i) ')'],['c' int2str(i)]);
    end
    fv_str = [fv_str(1:pos-1),arg,fv_str(pos+1:end)];
end
fv_str = strrep(fv_str,'f = @','fv = @');
eval(fv_str)

end

