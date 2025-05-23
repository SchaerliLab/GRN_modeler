function Mobj = get_model(cell)
% delete the extra input parameters and give the model object

% create new model
Mobj = sbiomodel(cell.data.name);
% copy the compartment
copyobj(cell.data.Mobj.Compartment,Mobj);

% copy the model
cell.add_model(cell.data.Mobj,Mobj);

% delete the extra input parameters
% in the regulators
regulator_names = fieldnames(cell.data.regulator_routes);
for i = 1:numel(regulator_names)
    regulation = cell.data.regulator_routes.(regulator_names{i});
    for j = 1:numel(regulation)
        % find the object
        regulated_obj = eval(cell.route2fieldnames(regulation{j}));
        regulated_obj.delete_input(Mobj);
    end
end
% nodes
node_names = fieldnames(cell.nodes);
for i = 1:numel(node_names)
    % delete the inputs of the object
    cell.nodes.(node_names{i}).delete_input(Mobj);
end
% proteases
protease_names = fieldnames(cell.proteases);
for i = 1:numel(protease_names)
    % delete the inputs of the object
    cell.proteases.(protease_names{i}).delete_input(Mobj);
end

% delete every unused property
delete(findUnusedComponents(Mobj));

% if we have an external model, we will add that too
if isvalid(cell.data.Mobj_extern)

    % if we have only one compartment, we will add it to the same
    % compartment as in our model
    if numel(cell.data.Mobj_extern.Compartments)==1 && numel(Mobj.Compartments)==1
        % rename the compartment
        cell.data.Mobj_extern.Compartments.Name = Mobj.Compartments(1).Name;
    end

    % copy the external model into our model
    cell.add_model(cell.data.Mobj_extern,Mobj);

end

% copy the configset
configsetobj_in = getconfigset(cell.data.Mobj);
cell.set_configset(Mobj,configsetobj_in);

% accelerate
if cell.data.Accelerate==true
    sbioaccelerate(Mobj);
end

end
