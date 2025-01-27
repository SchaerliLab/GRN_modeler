classdef node < bacteria_element
    % Building the gene regulatory network node by node

    properties (SetAccess = protected)

        % name of the object
        name;

        % node type
        type;

    end

    properties

        % struct with model for species, reactions, rules, parameters
        % (handle classes with the properties in the common model)
        Mobj Mobj_class;

        % mRNA production
        % list of the regulators
        regulators struct;

        % name of the protease objects which will degrade the proteins
        % produced by the given node
        protease cell;

    end

    methods

        % constructor
        function obj = node(name,type)

            % create default node model
            if isempty(fieldnames(obj.data.node_models))
                obj.node_model();
            end

            % if node_type is not given, we will use the first type
            if nargin < 2
                node_types = fieldnames(obj.data.node_models);
                type = node_types{1};
            end

            % set the properties
            obj.name = name;
            obj.type = type;

            % === rename individual properties ===
            % create new temporary model
            Mobj = sbiomodel(obj.data.name);
            % add a compartment
            addcompartment(Mobj,obj.data.compartment_name);

            % copy the node model into the temporary model
            obj.add_model(obj.data.node_models.(type),Mobj);

            % % rename every individual properties (except the one with the
            % % object name)
            % obj.rename_prop(Mobj);

            % rename every individual properties (except the one with the
            % object name) and the input
            exception = [Mobj.UserData.input];
            obj.rename_prop(Mobj,['_',obj.name],exception);

            % rename the input
            % we give "-" to get [] around the input name
            for i = 1:numel(Mobj.UserData.input)
                par = sbioselect(Mobj.Parameters,'Name',Mobj.UserData.input{i});
                if strcmp(par.Notes,'Individual')
                    rename(par,[par.Name,'_',obj.name,'-input']);
                end
            end

            % % change the species names according to the object name
            % for i = 1:numel(obj.Mobj_individ.Species)
            %     % due to the reactions, some 'common' species might be in
            %     % Mobj_individual, we neglect them everywhere
            %     if strcmp(Mobj.Species(i).Notes,'Individual')
            %         rename(Mobj.Species(i),[obj.Mobj.Species(i).Name,'_',obj.name]);
            %     end
            % end
            % % rename the parameters
            % for i = 1:numel(Mobj.Parameters)
            %     if strcmp(Mobj.Parameters(i).Notes,'Individual')
            %         rename(Mobj.Parameters(i),[Mobj.Parameters(i).Name,'_',obj.name]);
            %     end
            % end

            % add it to the common model and make a handle copy in obj.Mobj
            obj = obj.share_model(Mobj);

        end
       
    end
end