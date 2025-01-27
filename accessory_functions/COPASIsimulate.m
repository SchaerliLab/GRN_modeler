function [t,c,names] = COPASIsimulate(Mobj,solver,export_name)
%COPASISIMULATE run simulation in copasi
% Mobj: simbiology model object, which contains the simulation settings
% solver: type of the solver used in COPASI
% export_name: optional, the name of the exported SBML file in the tmp
% folder
% t: output time
% c: output concentrations
% name: smecies names according to the columns of c

% Mobj = correct_modifiers(convert2irrev(Mobj));

% % to avoid version problems
% pyenv("ExecutionMode","OutOfProcess");

%% Set up the necessary variables
% name of the exported model with full path
if nargin >=3
    model_name = [tempdir export_name];
else
    model_name = [tempdir 'model.sbml'];
end

% export SBML file into the temporary directory of the system
sbmlexport(Mobj, model_name);

% get simulation settings
configsetobj = getconfigset(Mobj);

% tolerance values
abstol = configsetobj.SolverOptions.AbsoluteTolerance;
reltol = configsetobj.SolverOptions.RelativeTolerance;

% followed species and time
StatesTolog = {'Time',configsetobj.RuntimeOptions.StatesToLog.Name};

% simulation times
StopTime = configsetobj.StopTime;
if ~isempty(configsetobj.SolverOptions.OutputTimes)
    StartTime = configsetobj.SolverOptions.OutputTimes(1);
else
    StartTime = 0;
end
times = configsetobj.SolverOptions.OutputTimes;

% maximum simulation time for the solver
MaximumWallClock = configsetobj.MaximumWallClock;

%% run the simulation

simdata = pyrunfile("run_COPASI.py","simdata",model_name=model_name,abstol=abstol,reltol=reltol,StartTime=StartTime,StopTime=StopTime,times=times,MaximumWallClock=MaximumWallClock,StatesTolog=StatesTolog,solver=solver);

% Extract column names
names = cellfun(@char, cell(simdata.columns.tolist()), 'UniformOutput', false);

% first name is Time
names(1) = [];

% Extract the data as a cell array
c = double(simdata.values);

% % get time data
% t = double(simdata.index.tolist()).';

% if we do not start from 0 time, we do not need the initial condition
% again
if c(1,1)<StartTime || (~isempty(times) && c(1,1)<times(1,1))
    c(1,:) = [];
end

% the first column is the time
t = c(:,1);
c(:,1) = [];

end

