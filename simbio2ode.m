function removed_species = simbio2ode(Mobj,ODE_name,type)
%SIMBIO2ODE create ODE function from Simbiology model
% Mobj: Simbiology model object
% ODE_name: name of the created ODE function
% type: 'normal': [dc1, dc2, dc3, dc4] = generated_model(time, c1, c2, c3, c4, p)
% type: 'cell': dc = generated_model(time, c) 
% (dc and c are cells: c1 = c{1}, dc{1} = dc1)
% type: 'range': c{1}(R), dc{1}(R) (reactions in a given range)
% removed_species: logical array, if a species is set by a rule,
% not by the ODEs, it will have a logical 1

% the default is the normal version
if nargin < 3
    type = 'normal';
end

% make a copy from the model
model = copyobj(Mobj);

% rename the species
switch type
    case 'normal'
        for i = 1:length(Mobj.Species)
            rename(model.Species(i),['c' int2str(i)]);
        end
    case 'cell'
        for i = 1:length(Mobj.Species)
            rename(model.Species(i),['c{' int2str(i) '}']);
        end
    case 'range'
        for i = 1:length(Mobj.Species)
            rename(model.Species(i),['c{' int2str(i) '}(R)']);
        end
    otherwise
        error('Unknown type')
end
% rename parameters
n_intern = 0;
n_extern = 0;
for i = 1:length(Mobj.Parameters)
    if model.Parameters(i).Constant == true
        % external parameters
        n_extern = n_extern+1;
        rename(model.Parameters(i),['p.k_ext' int2str(n_extern)]);
    else
        % internal parameters, changing in time, e.g. for Hill functions
        n_intern = n_intern+1;
        rename(model.Parameters(i),['k_int' int2str(n_intern)]);
    end
end
% rename reations
for i = 1:length(Mobj.Reactions)
    rename(model.Reactions(i),['r' int2str(i)]);
end

[odes, ~, ~] = getequations(model);
odes = strsplit(odes,'\n');

% % vectorize the system
% odes = strrep(odes,'/','./');
% odes = strrep(odes,'*','.*');
% odes = strrep(odes,'^','.^');

%% create the function

% we do not need the extension now
ODE_name = replace(ODE_name,'.m','');

% write the function
fileID = fopen([ODE_name,'.m'],'w');

% headline
switch type
    case 'normal'
        fprintf(fileID,'function [dc1');
        for i = 2:numel(model.species)
            fprintf(fileID,[', dc' int2str(i)]);
        end
        fprintf(fileID,['] = ', ODE_name,'(time']);
        for i = 1:numel(model.species)
            fprintf(fileID,[', c' int2str(i)]);
        end
        % "p" is a struct for the constant parameters
        fprintf(fileID,', p)\n');
    case 'cell'
        fprintf(fileID,['function dc = ', ODE_name,'(time, c, p)\n']);
    case 'range'
        fprintf(fileID,['function dc = ', ODE_name,'(time, c, p, R)\n']);
end

fprintf(fileID,'%% Automatically generated function from a SimBiology model\n');

% create persistent variables for internal variables (it is useful if these
% are matrices)
fprintf(fileID,'\n%% Persistent (inner) variables:\n');
for i = 1:numel(model.Reactions)
    fprintf(fileID,['persistent r' int2str(i) ';\n']);
end
for i = 1:n_intern
    fprintf(fileID,['persistent k_int' int2str(i) ';\n']);
end

% changes in the rules
if numel(model.Rules) > 0

    % remove the rules which are determining a given species and create a
    % separate function for them
    [model,removed_species] = rule_function(model,ODE_name,type);

    % reorder the rules that they can be calculated consecutively
    rule_order = reorder_rules(model);
else
    % otherwise there is no removed species from the rules and ODEs
    removed_species = false(numel(model.Species),1);
end

% print Rules
fprintf(fileID,'\n%% Rules:\n');
for i = 1:numel(model.Rules)
    fprintf(fileID,vectorize([model.Rules(rule_order(i)).Rule ';\n']));
end

% print Fluxes
fprintf(fileID,'\n%% Fluxes:\n');
for i = 1:numel(model.Reactions)
    fprintf(fileID,vectorize(['r' int2str(i) ' = ' model.Reactions(i).ReactionRate ';\n']));
end

myode = odes(2:find(strcmp(odes,'Fluxes:'),1,'first')-1);
% if a species is not changed, we write zeros in the output
if numel(myode) < numel(model.species)
    tmp = myode;
    myode = cell(numel(model.species),1);
    for i = 1:numel(model.species)
        switch type
            case 'normal'
                str = ['c' int2str(i)];
            case {'cell';'range'}
                str = ['c{' int2str(i) '}'];
        end
        position = find(contains(tmp,str)==1,1,'first');
        if ~isempty(position)
            myode{i} = tmp{position};
        else
            switch type
                case 'normal'
                    myode{i} = ['dc' int2str(i) ' = 0*c' int2str(i)];
                case 'cell'
                    myode{i} = ['[dc{' int2str(i) '}] = 0*[c{' int2str(i) '}]'];
                case 'range'
                    myode{i} = ['[dc{' int2str(i) '}] = 0*[c{' int2str(i) '}(R)]'];
            end
        end
    end
end

% print ODEs
fprintf(fileID,'\n%% ODEs:\n');
switch type
    case 'normal'
        for i = 1:numel(model.species)
            fprintf(fileID,vectorize(strrep([myode{i},';\n'],['d(c' int2str(i) ')/dt'],['dc' int2str(i)])));
        end
    case 'cell'
        for i = 1:numel(model.species)
            fprintf(fileID,vectorize(strrep([myode{i},';\n'],['d([c{' int2str(i) '}])/dt'],['dc{' int2str(i) '}'])));
        end
    case 'range'
        for i = 1:numel(model.species)
            fprintf(fileID,vectorize(strrep([myode{i},';\n'],['d([c{' int2str(i) '}(R)])/dt'],['dc{' int2str(i) '}'])));
        end
end

fprintf(fileID,'\nend\n');
fclose(fileID);

end

function str = vectorize(str)
% vectorize the system
str = strrep(str,'/','./');
str = strrep(str,'*','.*');
str = strrep(str,'^','.^');
end

