function [model,removed_species] = rule_function(model,ODE_name,type)
% remove the rules which are determining a given species and create a
% separate function for them
% ODE_name: the name of the output function: "ODE_name_rule"
% model: SimBiology model object
% removed_species: whether a species in the list is set by a rule or not
% type: 'normal': [dc1, dc2, dc3, dc4] = generated_model(time, c1, c2, c3, c4, p)
% type: 'cell': dc = generated_model(time, c) 
% (dc and c are cells: c1 = c{1}, dc{1} = dc1)
% type: 'range': c{1}(R), dc{1}(R) (reactions in a given range)

% number of the rules
nrule = numel(model.Rules);
% variables set by the rules
rule_names = cell(nrule,1);
% rules to remove
removed_rules = false(nrule,1);
% removed species
removed_species = false(numel(model.Species),1);
% name of the species
species_names = {model.Species.name};
for i = 1:nrule
    % find "=" sign
    pos = find(model.Rules(i).Rule=='=',1,'first')-2;
    rule_names{i} = model.Rules(i).Rule(1:pos);
    % get rid of the brackets
    rule_names{i} = strrep(rule_names{i},'[','');
    rule_names{i} = strrep(rule_names{i},']','');
    % if it is a species
    if any(strcmp(species_names,rule_names{i}))
        removed_rules(i) = true;
        removed_species(find(strcmp(species_names,rule_names{i})==1,1,"first")) = true;
    end
end
% create separate function for these species set by a rule
% write the function
fileID_rule = fopen([ODE_name,'_rule.m'],'w');
switch type
    case 'normal'
        fprintf(fileID_rule,'function [dc1');
        for i = 2:numel(model.species)
            fprintf(fileID_rule,[', c' int2str(i)]);
        end
        fprintf(fileID_rule,['] = ', ODE_name,'_rule(time']);
        for i = 1:numel(model.species)
            fprintf(fileID_rule,[', c' int2str(i)]);
        end
        % "p" is a struct for the constant parameters
        fprintf(fileID_rule,', p)\n');
    case 'cell'
        fprintf(fileID_rule,['function c = ', ODE_name,'_rule(time, c, p)\n']);
    case 'range'
        fprintf(fileID_rule,['function c = ', ODE_name,'_rule(time, c, p, R)\n']);
end
fprintf(fileID_rule,['%% Automatically generated function from a SimBiology model\n'...
    '%% to set the species determined by a rule\n']);

% print Rules
fprintf(fileID_rule,'\n%% Rules:\n');
for i = 1:numel(model.Rules)
    if removed_rules(i) == true
        fprintf(fileID_rule,vectorize([model.Rules(i).Rule ';\n']));
    end
end
fprintf(fileID_rule,'\nend\n');
fclose(fileID_rule);
% remove the rules
delete(model.Rules(removed_rules));
end

