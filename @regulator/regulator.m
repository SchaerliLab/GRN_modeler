classdef regulator < bacteria_element
    % Properties of the regulators (activation/inhibition by TFs or small
    % molecules)

    properties (SetAccess = protected)

        % name of the object
        name;

        % regulator type R/A: repression or acitvation
        type;

        % name of the 'HILL' parameter when the regulator will act
        Hill_target (1,:) char;

    end

    properties

        % struct with model for species, reactions, rules, parameters
        % (handle classes with the properties in the common model)
        Mobj Mobj_class;

        % regulators of the regulator
        regulators struct;

        % every object in the regulation chain,
        % eg. node A|-R1|-R2: [{'A'} {'HILL1'} {'|-'} {'R1'} {'HILL4'} {'<-'} {'R2'}]
        suffix (1,:) cell;

        % simplified suffix, the input name will be not present if it has
        % just one input and the object name will be present if it has just
        % one output
        suffix_simple (1,:) cell;

    end

    methods

        % constructor
        % HILL_target is where the regulator acts on the regulated object
        % the reactions will be added to the object in add_regulator()
        % function
        % suffix e.g.: [{'N1'} {'HILL1'} {'|-'} {'R1'} {'HILL3'} {'<-'} {'R2'}]
        function obj = regulator(name,type,suffix,suffix_simple,obj_regulated)

            % % create default regulator model
            % % if it is deleted or it is not set yet
            % if ~isfield(obj.data.regulator_models,type) || isempty(obj.data.regulator_models.(type))
            %     % create the node models
            %     obj.regulator_model();
            % end

            % set the properties
            obj.name = name;
            obj.suffix = suffix;
            obj.type = type;
            obj.suffix_simple = suffix_simple;

            % check whether the INPUT of the regulated object is individual
            % if it is a common parameter, we will not deal with it
            if isa(obj_regulated,'node')
                individual_INPUT = strcmp(get(sbioselect(obj.data.node_models.(obj_regulated.type),'Name',suffix{end-2}),'Notes'),'Individual');
            elseif isa(obj_regulated,'protease')
                individual_INPUT = strcmp(get(sbioselect(obj.data.protease_models.(obj_regulated.type),'Name',suffix{end-2}),'Notes'),'Individual');
            else
                individual_INPUT = strcmp(get(sbioselect(obj.data.regulator_models.(obj_regulated.type),'Name',suffix{end-2}),'Notes'),'Individual');
            end
            
            % if the input parameter is a 'common' parameter, we will not add the new individual 
            % input parameter, it can generate regulation just through the reactions
            % rename the hill target if it is an individual parameter
            if individual_INPUT==true
                if isa(obj_regulated,'node') || isa(obj_regulated,'protease')
                    % rename it with the object name
                    obj.Hill_target = [suffix{end-2},'_',obj_regulated.name,'-input'];
                else
                    % rename it with the suffix
                    obj.Hill_target = [suffix{end-2},'_',obj_regulated.suffix_simple{:},'-input'];
                end
            end

            % === rename individual properties ===
            % create new temporary model
            Mobj = sbiomodel(obj.data.name);
            % add a compartment
            addcompartment(Mobj,obj.data.compartment_name);

            % copy the regulator model into the temporary model
            obj.add_model(obj.data.regulator_models.(type),Mobj);

            % === Renaming regulator elements ===
            % change the regulator name
            for i = 1:numel(Mobj.Species)
                if contains(Mobj.Species(i).Name,'REGULATOR')
                    new_name = strrep(Mobj.Species(i).Name,'REGULATOR',name);
                    rename(Mobj.Species(i),new_name);
                end
            end

            % rename every properties with suffix
            % except the input, OUTPUT and the species
            exception = [Mobj.UserData.input,Mobj.species.Name,{'OUTPUT'}];
            obj.rename_prop(Mobj,['_',suffix_simple{:}],exception);

            % rename the input
            % we give "-" to get [] around the input name
            for i = 1:numel(Mobj.UserData.input)
                par = sbioselect(Mobj.Parameters,'Name',Mobj.UserData.input{i});
                if strcmp(par.Notes,'Individual')
                    rename(par,[par.Name,'_',obj.suffix_simple{:},'-input']);
                end
            end

            % rename the output of the regulator according to the input of
            % the regulated object
            if individual_INPUT==true % just for individual INPUTs
                OUTPUT_name = [suffix{end-2},'_',suffix_simple{:}];
                output_parameter = sbioselect(Mobj.Parameters,'Name','OUTPUT');
                if ~isempty(output_parameter)
                    rename(output_parameter,OUTPUT_name);
                end
            end

            % rename the individual species according to the regulated object name
            % (and not with the suffix!)
            % the 'REGULATED' word will be replaced with the object name
            % e.g. dCas:REGULATOR:DNA_REGULATED => dCas:REGULATOR:DNA_N1
            for i = 1:numel(Mobj.Species)
                if strcmp(Mobj.Species(i).Notes,'Individual') && any(contains(Mobj.Species(i).Name,'REGULATED'))
                    rename(Mobj.Species(i),strrep(Mobj.Species(i).Name,'REGULATED',obj_regulated.name));
                end
            end

            % add it to the common model and make a handle copy in obj.Mobj
            obj = obj.share_model(Mobj);

            % === Renaming regulated elements ===
            % if the INPUT is individual, we are modifying it in the
            % reactions and rules
            if individual_INPUT == true
                % rename the HILL parameter according to the regulator
                % we keep a *HILL_obj.name term (after the HILL parameter of the regulator)
                % because there might be more regulators
                % pattern = [obj.Hill_target '(?!_)']; % hill_name cannot continue with '_'
                pattern = ['[' obj.Hill_target ']'];

                % change the 'input' in the reactions
                for i = 1:numel(obj_regulated.Mobj.Reactions)
                    % e.g. *HILL_A_input => *HILL_A_R1_R2*HILL_A_input
                    obj_regulated.Mobj.Reactions(i).Reactionrate =...
                        strrep(get(obj_regulated.Mobj.Reactions(i),'Reactionrate'),pattern,...
                        ['[' OUTPUT_name ']*[' obj.Hill_target ']']);
                end
                % rules
                for i = 1:numel(obj_regulated.Mobj.Rules)
                    % e.g. *HILL_A_input => *HILL_A_R1_R2*HILL_A_input
                    obj_regulated.Mobj.Rules(i).Rule =...
                        strrep(get(obj_regulated.Mobj.Rules(i),'Rule'),pattern,...
                        ['[' OUTPUT_name ']*[' obj.Hill_target ']']);
                end
            end

        end

    end
end

