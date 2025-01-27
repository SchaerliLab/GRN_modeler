classdef protease < bacteria_element
    % The protein degradation can go with Michaelis-Menten kinetics

    properties (SetAccess = protected)

        % name of the object
        name;

        % protease type
        type;

    end

    properties

        % struct with model for species, reactions, rules, parameters
        % (handle classes with the properties in the common model)
        Mobj Mobj_class;

        % list of the nodes which proteins will be degraded with the given
        % protease
        node_list cell = {};

        % list of the species which will be degraded with the given protease
        species_list cell = {};

        % list of the regulators
        regulators  struct;

    end

    methods

        % ctor
        % name: the name of the protease
        % if we want to add a node, the reactions will be added to the node
        % in the add_node2protease() function!
        function obj = protease(name,type,node_name)

            % create default protease model
            if isempty(fieldnames(obj.data.protease_models))
                obj.protease_model();
            end

            % if type is not given, we will use the first type
            if nargin < 2 || isempty(type)
                protease_types = fieldnames(obj.data.protease_models);
                type = protease_types{1};
            end

            % set the properties
            obj.name = name;
            obj.type = type;

            % === rename individual properties ===
            % create new temporary model
            Mobj = sbiomodel(obj.data.name);
            % add a compartment
            addcompartment(Mobj,obj.data.compartment_name);

            % copy the protease model into the temporary model
            obj.add_model(obj.data.protease_models.(type),Mobj);

            % change the protease name
            for i = 1:numel(Mobj.Species)
                if strcmp(Mobj.Species(i).Tag,'PROTEASE')
                    rename(Mobj.Species(i),name);
                end
                % the species will be renamed according to the node name
                % except the protease
                if nargin >= 3 && ~strcmp(Mobj.Species(i).Name,name)
                    rename(Mobj.Species(i),[Mobj.Species(i).Name,'_',node_name]);
                end
            end

            % % rename every individual properties except the one with the
            % % object name and the species, they are renamed according
            % % to the node name
            % exception = {Mobj.Species.Name};
            % obj.rename_prop(Mobj,['_',obj.name],exception);

            % rename every individual properties except the one with the
            % object name,the input and the species, they will be renamed according
            % to the node name
            exception = [Mobj.UserData.input,{Mobj.Species.Name}];
            obj.rename_prop(Mobj,['_',obj.name],exception);

            % rename the input
            % we give "-" to get [] around the input name
            for i = 1:numel(Mobj.UserData.input)
                par = sbioselect(Mobj.Parameters,'Name',Mobj.UserData.input{i});
                if strcmp(par.Notes,'Individual')
                    rename(par,[par.Name,'_',obj.name,'-input']);
                end
            end

            % if the node name is given, rename the reactions
            if nargin >= 3
                for i = 1:numel(Mobj.Reactions)
                    set(Mobj.Reactions(i),'Name',['Protease_' node_name '_' name '_' int2str(i)]);
                end
            end

            % add it to the common model and make a handle copy in obj.Mobj
            obj = obj.share_model(Mobj);

        end

    end

end

