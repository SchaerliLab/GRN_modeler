function [] = reg2graph(cell,obj,weight)
% add the regulator (node in the graph) and its interactions (edges) to the
% graph
% if the regulator was produced by a node, add that inteartion as well
% go through on every regulator of the regulator recursively
% weight.regulator: weight for the regulator

% inputs of the given object
input_names = fieldnames(obj.regulators);

for n_input = 1:numel(input_names)

    % structure with the regulators for the actual input
    obj_input_regs = obj.regulators.(input_names{n_input});

    % regulators for the given input of the given object
    regnames = fieldnames(obj_input_regs);

    % go through on every regulators
    for i = 1:numel(regnames)

        if obj_input_regs.(regnames{i}).type(1)=='A'
            % activation: positive weight
            weight_regulator = abs(weight.regulator);
        elseif obj_input_regs.(regnames{i}).type(1)=='R'
            % repression negative weight
            weight_regulator = -abs(weight.regulator);
        else
            % other regulation type
            weight_regulator = -2*abs(weight.regulator);
        end

        suffix_simple = obj_input_regs.(regnames{i}).suffix_simple;

        % regulated name
        regulated_name = [suffix_simple{1:end-3}];
        % regulation name
        regulation_name = [suffix_simple{1:end}];
        % regulator name
        regulator_name = [suffix_simple{end}];

        % if the regulator is part of a node we will use the node name for
        % the regulator
        name = cell.isNodeProduced(regulator_name);
        if ~isempty(name) % check if it is coming from a given node
            % use the node name
            regulator_name = name;
        end

        % add them as nodes if they are not present yet
        % add the regulator as a graph node if it is new
        if ~any(strcmp(obj.data.G.Nodes.Name,regulated_name))
            obj.data.G = addnode(obj.data.G,regulated_name);
        end
        if ~any(strcmp(obj.data.G.Nodes.Name,regulation_name))
            obj.data.G = addnode(obj.data.G,regulation_name);
        end
        if ~any(strcmp(obj.data.G.Nodes.Name,regulator_name))
            obj.data.G = addnode(obj.data.G,regulator_name);
        end

        % add the edge between regulator and regulation
        obj.data.G = addedge(obj.data.G,regulator_name,regulation_name,weight_regulator);
        % add the edge between regulation and regulated
        obj.data.G = addedge(obj.data.G,regulation_name,regulated_name,weight_regulator);

        % go through on every regulator of the regulator recursively
        cell.reg2graph(obj_input_regs.(regnames{i}),weight);

    end
end
end