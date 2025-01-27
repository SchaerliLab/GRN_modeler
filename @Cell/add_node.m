function [cell,name] = add_node(cell,name,type)
% add node to the cell

% give a name to the node if it is not provided
if nargin == 1 || isempty(name)
    names = {cell.data.Mobj.Species.Name};
    node_number = 1;
    name = 'N1';
    % we are looking for a new name
    while any(contains(names,name))
        node_number = node_number+1;
        name = ['N' int2str(node_number)];
    end
end

if any(strcmp(name,fieldnames(cell.nodes)))
    % if the node exists, we are not doing anything
    warning('There was already a node with name ''%s''!',name)
    return;
end

% create the new node
if nargin >= 3
    cell.nodes(1).(name) = node(name,type);
else
    cell.nodes(1).(name) = node(name);
end
end