function cell = add_node2protease(cell,node_name,prot_name)
% add node to protease
% add a node to the node_list of the protease
% put the reactions (due to the protease) into the node reactions
% extend the protease "rule"


% if there was an unmodified protease model before (without added nodes)
% we delete its individual reactions, species, rules etc, because they are not
% assigned to nodes and we start with an empty model
if isempty(cell.proteases.(prot_name).node_list)
    cell.proteases.(prot_name).Mobj.delete_individual();
    cell.proteases.(prot_name).Mobj = Mobj_class();
end

% type of the protease
type = cell.proteases.(prot_name).type;

% temporary protease object to set the properties according to
% the new node. Later on this object will be added to the whole
% protease model
tmp = protease(cell.proteases.(prot_name).name,type,node_name);


% === setting up things ===
% % add the new model with the new node to the protease model
% add_model(cell.proteases.(prot_name),tmp.Mobj_individ,cell.proteases.(prot_name).Mobj_individ);
% % add the common reactions and specieses

% add the species to the species list
for i = 1:numel(tmp.Mobj.Species)
    if ~strcmp(tmp.Mobj.Species(i).Name,prot_name)
        cell.proteases.(prot_name).species_list{end+1} = tmp.Mobj.Species(i).Name;
    end
end

% add the new model (with the new node) to the protease model
cell.proteases.(prot_name).Mobj = cell.proteases.(prot_name).Mobj.extend_model(tmp.Mobj);

% extend the substrate list in the rule
set_proteaserule(cell,prot_name);

% add the node to the node list
cell.proteases.(prot_name).node_list{end+1} = node_name;

end