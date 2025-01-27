classdef bacteria_element < bacteria
    % Common properties for the objects in the bacteria (e.g. obj,
    % regulator, protease)

    properties (SetAccess = protected, Abstract)

        name (1,:) char; % name of the object

        % type of the object
        type (1,:) char;

    end

    properties (Access = public, Abstract)

        % struct with model for species, reactions, rules, parameters
        % (handle classes with the properties in the common model)
        Mobj Mobj_class;

    end


    methods (Access = public)

        % ctor: create empty properties
        function obj = bacteria_element()
            obj.Mobj = Mobj_class();
        end

        % add model (Mobj) to the shared common model
        % make a (handle class) copy for the properties in the actual
        % object (obj.Mobj)
        obj = share_model(obj,Mobj);

        % rename individual properties
        rename_prop(obj,Mobj,suffix,exception);

    end

    methods (Access = public)

        % delete input ('HILL' parameters) from the reactions and rulesin the Mobj model
        delete_input(obj,Mobj);

        % delete individual parameters, reactions, rules, species from Mobj
        delete_individual(obj);

        % clear the deleted properties from the model
        obj = clean_up_obj(obj);

    end
end