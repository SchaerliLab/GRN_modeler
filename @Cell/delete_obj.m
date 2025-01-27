function cell = delete_obj(cell,obj_name)
% delete objects (iteratively for the regulators)
% with regulators, if obj_name contains the complete route of
% the regulator e.g. {'node1','HILL1','reg1','HILL','reg2'}, then 'reg2' will be
% deleted just from this regulation (and the further regulatiors for
% 'reg2' in this interaction will be deleted too),
% otherwise (e.g. obj_name='reg2') it will be deleted from every interactions
% if the first name is a node name and the second one is a protease name
% the node will be deleted only from the protease list
% e.g.: {'node_name','protease_name'}

% we need row cell arrays
if iscell(obj_name)
    obj_name = reshape(obj_name,[1,numel(obj_name)]);
end

% if the object name comtains '<-' or '|-', get rid of it
obj_name(strcmp(obj_name,'<-')|strcmp(obj_name,'|-')) = [];

% if it is a regulator name, we will have to delete it from
% every regulation chain
if isfield(cell.data.regulator_routes,obj_name)
    % delete the regulator from every route where it has a role
    routes = cell.data.regulator_routes.(obj_name);
    for i = 1:numel(routes)
        cell = delete_obj(cell,routes{i});
    end
    return;
elseif iscell(obj_name) && numel(obj_name) == 2 && cell.isProtease(obj_name{2}) && cell.isNode(obj_name{1})
    % if it is a protease - node interaction what we want to delete

    % delete the node from the protease list
    if ~isempty(strcmp(cell.nodes.(obj_name{1}).protease,obj_name{2}))

        % === delete node stuff from protease ===
        % the nodes interacting with the protease
        node_list = cell.proteases.(obj_name{2}).node_list;
        cell.proteases.(obj_name{2}).node_list = ...
            node_list(~strcmp(node_list,obj_name{1}));
        % delete the species from node
        species_list = cell.proteases.(obj_name{2}).species_list;
        cell.proteases.(obj_name{2}).species_list = ...
            species_list(~contains(species_list,obj_name{1}));
        % rebuild the substrate rule of the protease
        set_proteaserule(cell,obj_name{2});

        % === delete protease stuff from the node ===
        % we have to delete the reactions which we added to the
        % node with the given protease
        delete(cell.nodes.(obj_name{1}).Mobj_individ.Reactions(contains({cell.nodes.(obj_name{1}).Mobj_individ.Reactions.Name},...
            ['Protease_' obj_name{1} '_' obj_name{2}])));

        % delete it from the node's protease list
        cell.nodes.(obj_name{1}).protease(strcmp(cell.nodes.(obj_name{1}).protease,obj_name{2})) = [];

    end
    return;
end

% find the object
[str,success] = cell.route2fieldnames(obj_name);
if success == false
    return;
end
obj = eval(str);

% delete the regulators of the object recursively
% inputs of the given object
input_names = fieldnames(obj.regulators);
for n_input = 1:numel(input_names)
    % structure with the regulators for the actual input
    obj_input_regs = obj.regulators.(input_names{n_input});
    regulator_names = fieldnames(obj_input_regs);
    for i = 1:numel(regulator_names)
        % the route for the regulator e.g. {'node1','HILL2','reg2'}
        if iscell(obj_name)
            route = [obj_name,input_names{n_input},regulator_names{i}];
        else
            route = {obj_name,input_names{n_input},regulator_names{i}};
        end
        % delete the regulator and its regulators
        cell = delete_obj(cell,route);
    end
end

% identify the object and delete it
switch class(obj)
    case 'node'
        % delete the node from the protease list
        if ~isempty(cell.nodes.(obj_name).protease)
            prot_names = cell.nodes.(obj_name).protease;
            for i = 1:numel(prot_names)
                % delete the node from the protease list
                node_list = cell.proteases.(prot_names{i}).node_list;
                cell.proteases.(prot_names{i}).node_list = ...
                    node_list(~strcmp(node_list,obj_name));
                % delete the species from node
                species_list = cell.proteases.(prot_names{i}).species_list;
                cell.proteases.(prot_names{i}).species_list = ...
                    species_list(~contains(species_list,obj_name));
                % rebuild the substrate rule of the protease
                set_proteaserule(cell,prot_names{i});
            end
        end
        % delete every regulation where the protein from this
        % node take part
        for node_prefix = cell.nodes.(obj_name).Mobj.UserData.output_prefix
            regulator_from_the_node = [node_prefix{1}, obj_name];
            if isfield(cell.data.regulator_routes,regulator_from_the_node)
                cell = cell.delete_obj(regulator_from_the_node);
            end
        end
        % delete the individual part of the node model
        obj.delete_individual();
        % delete the node from the node list
        cell.nodes = rmfield(cell.nodes,obj_name);

    case 'protease'
        % we have to delete the reactions which we added to the
        % node with the given protease
        node_names = fieldnames(cell.nodes);
        for i = 1:numel(node_names)
            % if the node is effected by the protease
            if any(strcmp(cell.nodes.(node_names{i}).protease,obj_name))
                delete(cell.nodes.(node_names{i}).Mobj_individ.Reactions(contains({cell.nodes.(node_names{i}).Mobj_individ.Reactions.Name},...
                    ['Protease_' node_names{i} '_' obj_name])));
            end
        end
        % delete the individual part of the node model
        obj.delete_individual();
        % delete protease from the cell protease list
        cell.proteases = rmfield(cell.proteases,obj_name);
        % delete it from the node's protease list
        node_names = fieldnames(cell.nodes);
        for i = 1:numel(node_names)
            cell.nodes.(node_names{i}).protease(strcmp(cell.nodes.(node_names{i}).protease,obj_name)) = [];
        end

    case 'regulator'
        % % if the regulator was given by a node name, we should
        % % add the prefix
        % if isfield(cell.nodes,obj_name{end})
        %     obj_name{end} = [cell.data.node_prefix,obj_name{end}];
        % end

        % delete the regulator from the regulator list
        n_regulator_routes = numel(cell.data.regulator_routes.(obj_name{end}));
        position = true(n_regulator_routes,1);
        for i = 1:n_regulator_routes
            regrout = cell.data.regulator_routes.(obj_name{end}){i};
            regrout(strcmp(regrout,'<-')|strcmp(regrout,'|-')) = [];
            if numel(regrout) ~= numel(obj_name)
                position(i) = false;
            else
                position(i) = all(strcmp(regrout,obj_name));
            end
        end
        cell.data.regulator_routes.(obj_name{end})(position) = [];
        % if the regulator list is empty, we do not have more of this
        % regulator, delete the list element
        if isempty(cell.data.regulator_routes.(obj_name{end}))
            cell.data.regulator_routes = rmfield(cell.data.regulator_routes,obj_name{end});
            cell.data.regulators = rmfield(cell.data.regulators,obj_name{end});
        end

        % delete the individual part of the node model
        obj.delete_individual();

        % delete regulator from the list of the parent
        [str,success] = route2fieldnames(cell,obj_name{1:end-2});
        if success == false
            return;
        end
        eval([str,'.regulators.(''' obj_name{end-1} ''')= rmfield(' str '.regulators.(''' obj_name{end-1} '''),''' obj_name{end} ''');']);
        % if the input does not contain any more regulators, delete the
        % input field
        if isempty(fieldnames(eval([str,'.regulators.(''' obj_name{end-1} ''')'])))
            eval([str,'.regulators = rmfield(' str '.regulators,''' obj_name{end-1} ''');']);
        end

    otherwise
        warning('We do not find the object (%s) what you want to delete.',obj_name)
        return;
end

end