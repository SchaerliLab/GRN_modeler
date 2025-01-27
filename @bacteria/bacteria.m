classdef bacteria

    % common variable in SharedData: the simbiology model object
    properties (Constant)
        data = SharedData();
    end

    properties (Abstract) %, SetAccess = {?SharedData,?bacteria,?Cell,?node,?regulator,?protease})
        regulators struct;
    end

    methods

        % create model
        create_model(obj);

        % set model parameters
        set_model(obj);
            
        % add regulator to the object (eg. node, protease) with a given 'name'
        % 'type': 'A': activation, 'R': repression
        % 'suffix': to create individual parameter names
        % 'suffix_simple': input, putput is present only if there are more
        % of them
        function obj = add_regulator(obj,name,type,suffix,suffix_simple,obj_regulated)

            % create the regulator
            % suffix{end-2} is the selected input of the object
            obj.regulators(1).(suffix{end-2}).(name) = regulator(name,type,suffix,suffix_simple,obj_regulated);
            
        end

        % add simbio model (Mobj_in) to another model (Mobj_out)
        % (species, reactions, parameters, rules) to the common cell model 
        % if 'notes' is given, the copied objects must have this specific
        % 'Notes'
        add_model(~,Mobj_in,Mobj_out,notes);

        % create node model
        node_model(obj);

        % create protease model
        protease_model(obj);

        % create regulator model
        regulator_model(obj);
        
        % checks wheter the fields in a nested structure or object exist or not
        success = checkstruct(obj,str,start_pos);

        % create fieldnames for the given route to a certain object
        [str,success] = route2fieldnames(obj,varargin);

        % check the model properties
        check_model(obj,Mobj);

        % set the userdata of the model
        Mobj = set_input_output(obj,Mobj,obj_type);

        % calculate the regulator routes
        regrouts = collect_regrouts(obj);

        % if the 'input' is not specified in the route, put the first
        % possible 'input' in it (except for the last one)
        route_out = addinput2route(obj,route);

        % get the original node name without a prefix
        original_name = original_node_name(obj,name);

    end
end