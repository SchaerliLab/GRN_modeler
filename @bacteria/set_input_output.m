function [Mobj] = set_input_output(~,Mobj,obj_type)
% set the userdata of the model to see the inputs (HILL parameters) and
% outputs (species acting as regulators) of the node
% and set the node prefix

% input names 'Tag': 'INPUT'
param_tags = {Mobj.Parameters.Tag};
param_names = {Mobj.Parameters.Name};
data.input = param_names(strcmp(param_tags,'INPUT'));

switch obj_type
    case 'node'
        % Output names 'Tag': 'REGULATOR'
        species_tags = {Mobj.Species.Tag};
        species_names = {Mobj.Species.Name};
        data.output_prefix = species_names(strcmp(species_tags,'REGULATOR'));
        % add '_' to the prefix:
        data.output_prefix = strcat(data.output_prefix,'_');
%     case 'protease'
end

% save the input/output variables into the UserData
Mobj.UserData = data;
        
end

