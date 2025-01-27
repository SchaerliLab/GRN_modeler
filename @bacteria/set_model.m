function set_model(obj)

% basic parameters
obj.data.Mobj.Name = obj.data.name;
obj.data.Mobj.Compartment(1).Name = obj.data.compartment_name;
obj.data.Mobj.Compartment(1).Value = obj.data.volume.Value;
obj.data.Mobj.Compartment(1).Units = obj.data.volume.Units;

% % simulation settings
% configsetobj = getconfigset(obj.data.Mobj);
% set(get(getconfigset(obj.data.Mobj),'CompileOptions'),...
%     'Dimensionalanalysis',obj.data.Dimensionalanalysis,'UnitConversion',false)
% set(configsetobj,'SolverType',obj.data.SolverType,'Stoptime',obj.data.Stoptime)
% set(get(configsetobj,'SolverOptions'),'AbsoluteTolerance',obj.data.AbsoluteTolerance,...
%     'RelativeTolerance',obj.data.RelativeTolerance)
% 
% if obj.data.Accelerate == true
%     sbioaccelerate(obj.data.Mobj)
% end

% % followed species
% % if a species is not in the model, we will delet it from StatesToLog
% nonmodel_name = false(numel(obj.data.StatesToLog),1);
% for i = 1:numel(obj.data.StatesToLog)
%     nonmodel_name(i) = ~any(strcmp(obj.data.StatesToLog(i),{obj.data.Mobj.Species.Name}));
% end
% obj.data.StatesToLog(nonmodel_name) = [];
% set(get(configsetobj,'RuntimeOptions'),'StatesToLog',obj.data.StatesToLog)

end