function protease_model(obj)
% add reactions to the protease

% number of function calls
persistent n_call;
if isempty(n_call)
    n_call = 1;
else
    n_call = n_call+1;
end

% create a model
Mobj = sbiomodel('protease_model');
% add a compartment
addcompartment(Mobj,obj.data.compartment_name);
% % PROTEASE
% addspecies(Mobj,'PROTEASE','InitialAmount',0,...
%             'Units','molecule','Notes','Individual','Constant',true);

% add the model with the actual type
[Mobj,protease_types] = feval(obj.data.model_name,Mobj,'protease',n_call);

% if we do not have a model (empty type) we will do nothing
if isempty(protease_types)
    return;
end

% check the model
obj.check_model(Mobj);

% set the userdata of the model
Mobj = obj.set_input_output(Mobj,'protease');

% add the protease model to the common list
obj.data.protease_models(1).(protease_types{n_call}) = Mobj;
% call the function to put in the other models
if numel(protease_types) > n_call
    obj.protease_model();
else
    n_call = 0;
end

end

