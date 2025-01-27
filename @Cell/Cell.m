classdef Cell < bacteria 
    %CELL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties

        % nodes in the cell
        nodes struct;

        % protease enzymes in the cell
        proteases struct;

        % we do not have cellular level regulators (yet)
        % the regulators are stored in SharedData
        regulators struct;

    end
    
    methods

        % ctor
        function obj = Cell(model_name)

            % add the GRN folder and subfolders to the search path
            % find the correct folder
            [folderPath, ~, ~] = fileparts(which('Cell'));
            % remove the end
            folderPath = folderPath(1:end-6);
            addpath(genpath(folderPath))

            % name of the reaction kinetics model for the nodes
            if nargin == 0
                % default model
                obj.data.model_name = 'Tomazou';
            else
                obj.data.model_name = model_name;
            end

            % clean up the regulator lists
            obj.data.clean_start();

            % create a new model
            obj.create_model();

            % original settings
            obj.set_model();

            % set up the models
            obj.node_model();
            obj.protease_model();
            obj.regulator_model();
        end

        % delete objects (iteratively for the regulators)
        cell = delete_obj(cell,obj_name);

        % decide whether object is a node or not
        function isnode = isNode(cell,name)
            % if it is a node object
            if isa(name,'node')
                isnode = true;
                return;
            end
            % if it is the name of a node
            if iscell(name) && numel(name)>1
                isnode = false;
                return;
            end
            isnode = isfield(cell.nodes,name);
        end

        % decide whether the name is a node output or not
        function isnode = isNode_output(cell,name)
            isnode = false;
            node_names = fieldnames(cell.nodes);
            for i = 1:numel(node_names)
                % if it might be connected to a node
                if contains(name,node_names{i})
                    output_names = cell.nodes.(node_names{i}).Mobj.UserData.output_prefix;
                    for j = 1:numel(output_names)
                        % if it has the "outputnodename" structure
                        if strcmp(name,[output_names{j},node_names{i}])
                            isnode = true;
                            return;
                        end
                    end
                end
            end
        end

        % decide whether object is a protease or not
        function isprotease = isProtease(cell,name)
            % if it is a protease object
            if isa(name,'protease')
                isprotease = true;
                return;
            end
            % if it is the name of a protease
            if iscell(name) && numel(name)>1
                isprotease = false;
                return;
            end
            isprotease = isfield(cell.proteases,name);
        end

        % decide whether object is a regulator or not
        function isregulator = isRegulator(cell,name)
            % if it is a regulator object
            if isa(name,'regulator')
                isregulator = true;
                return;
            end
            % if it is the name of a regulator
            if iscell(name)
                % we are using row cell arrays for regulator names
                name = reshape(name,[1,numel(name)]);
            end
            isregulator = isfield(cell.data.regulator_routes,name);
        end

        % decide whether a given species comes from a node output or not
        % and give back the name of the node
        function name = isNodeProduced(cell,variable_name)
            % go through on every possible node
            node_names = fieldnames(cell.nodes);
            for j = 1:numel(node_names)
                node_name = node_names{j};
                % check if it is coming from a given node
                if contains(variable_name,node_name)
                    % check whether the name has the ['prefix','nodename'] structure
                    output_prefix = cell.nodes.(node_name).Mobj.UserData.output_prefix;
                    for k = 1:numel(output_prefix)
                        if strcmp(variable_name,[output_prefix{k},node_name])
                            % return the name of the node
                            name = node_name;
                            return
                        end
                    end
                end
            end
            % if it is not coming from a node, we return an empty variable
            name = [];
            return
        end

        % add node to the cell
        [cell,node_name] = add_node(cell,name,type);

        % add regulator to an object
        [cell,suffix] = add_regulator(cell,type,unit,varargin);

        % add protease to the cell and connect it to a given node
        [cell,protease_name] = add_protease(cell,node_name,protease_name,type);

        % get the property
        value = get(cell,name,PropertyName,varargin);

        % set the property
        set(cell,name,PropertyName,PropValue,varargin);

        % === graph creation ===
        % create and plot graph
        % handle ox the axes to plot (optional)
        h = make_graph(obj,ha);
        
        % add regressior and its edges to the gragh
        reg2graph(cell,obj,weight);

        % add a node to a protease
        cell = add_node2protease(cell,node_name,prot_name);

        % set the rule for the substrates for a protease
        set_proteaserule(cell,prot_name);

        % clean up every deleted model element in every object
        cell = clean_up_cell(cell);

        % delete the extra input parameters and give the model object
        Mobj = get_model(cell);

        % set the configset in a model according to another configset
        set_configset(cell,Mobj,configsetobj_in);

    end
end

