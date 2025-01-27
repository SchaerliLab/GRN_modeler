function regulator_model(obj)
% build the regulator models recursively

% number of regulator type
persistent n_type;
if isempty(n_type)
    n_type = 1;
else
    n_type = n_type+1;
end

% create a model
Mobj = sbiomodel('regulator_model');
% add a compartment
addcompartment(Mobj,obj.data.compartment_name);

% add the model with the actual type
[Mobj,regulator_types] = feval(obj.data.model_name,Mobj,'regulator',n_type);

% if we do not have a model (empty type) we will do nothing
if isempty(regulator_types)
    return;
end

% check the model
obj.check_model(Mobj);

% set the userdata of the model
Mobj = obj.set_input_output(Mobj,'regulator');

% add the regulator model to the common list
obj.data.regulator_models(1).(regulator_types{n_type}) = Mobj;

% call the function to put in the other models
if numel(regulator_types) > n_type
    obj.regulator_model();
else
    n_type = 0;
end

end


