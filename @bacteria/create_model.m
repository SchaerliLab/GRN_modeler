function create_model(obj)
% create model

% ========== build the model ==========

% create new model
obj.data.Mobj = sbiomodel(obj.data.name);
% add a compartment
addcompartment(obj.data.Mobj,obj.data.compartment_name);

% ======== add common parameters =======
% after sbioreset simbio models might be deleted, set up them
obj.data.set_Mobj();

end