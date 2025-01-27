classdef Mobj_class
    properties
        UserData    struct;
        Species     SimBiology.Species;
        Reactions   SimBiology.Reaction;
        Rules       SimBiology.Rule;
        Parameters  SimBiology.Parameter;
    end

    methods
        
        % extend the model with another Mobj_class object
        obj = extend_model(obj,Mobj);

        % delete model
        [] = delete(obj);

        % delete individual properties
        [] = delete_individual(obj);

        % clear the deleted properties
        obj = clean_up(obj);
        
    end

end