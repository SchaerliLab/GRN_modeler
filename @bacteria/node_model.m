function node_model(obj)
% build the default node models recursively

% number of function calls
persistent n_call;
if isempty(n_call)
    n_call = 1;
else
    n_call = n_call+1;
end

% create a model
Mobj = sbiomodel('node_model');
% add a compartment
addcompartment(Mobj,obj.data.compartment_name);

% add the model with the actual type
[Mobj,node_types] = feval(obj.data.model_name,Mobj,'node',n_call);

% if we do not have a model (empty type) we will do nothing
if isempty(node_types)
    return;
end

% check the model
obj.check_model(Mobj);

% set the userdata of the model
Mobj = obj.set_input_output(Mobj,'node');

% add the node model to the common list
obj.data.node_models(1).(node_types{n_call}) = Mobj;
% call the function to put in the other models
if numel(node_types) > n_call
    obj.node_model();
else
    n_call = 0;
end