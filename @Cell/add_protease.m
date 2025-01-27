function [cell,protease_name] = add_protease(cell,node_name,protease_name,type)
% add protease to the cell and connect it to a given node
% if the protease already exists we will just give the node

if isempty(node_name)
    % we are not doing anything without a node
    return;
end

% give a name to the protease if it is not provided
if nargin == 1 || isempty(protease_name)
    names = {cell.data.Mobj.Species.Name};
    node_number = 1;
    protease_name = 'PROT1';
    % we are looking for a new name
    while any(contains(names,protease_name))
        node_number = node_number+1;
        protease_name = ['PROT' int2str(node_number)];
    end
end

if ~isfield(cell.proteases,protease_name)
    % create the protease
    if nargin >= 4
        cell.proteases(1).(protease_name) = protease(protease_name,type);
    else
        cell.proteases(1).(protease_name) = protease(protease_name);
    end
end
if ~isempty(node_name)
    % if the node already contains the protease, we will not
    % add it again
    if any(strcmp(cell.nodes.(node_name).protease,protease_name))
        % if the interaction exists, we are not doing anything
        warning('There was already a node ''%s'' - protease ''%s'' interaction!',node_name,protease_name)
        return;
    end
    % add node reactions to the protease
    % cell = cell.proteases.(protease_name).add_node(node_name,cell);
    cell = cell.add_node2protease(node_name,protease_name);
    % add protease to the node
    cell.nodes.(node_name).protease{end+1} = protease_name;
end

end