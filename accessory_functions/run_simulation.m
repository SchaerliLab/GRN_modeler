function [t,c,names] = run_simulation(Mobj,solver,name)
%RUN_SIMULATION
% Mobj: simbiology model object
% solver: type of the solver, if it is not given, the solver given in Mobj
% will be used
% name: in the temporary epxorted model for the COPASI simulation we will 
% use this name

if nargin == 1
    % if no solver is selected, we will use the one selected in Mobj
    matlabsolver = true;
elseif any(contains({'ode15s','ode45','sundials'},solver)) 
    % MATLAB simulation with the given solver
    matlabsolver = true;
    configset = getconfigset(Mobj);
    if ~strcmp(configset.SolverType,solver)
        set(configset,'SolverType',solver);
    end
else
    % run COPASI simulation
    matlabsolver = false;
end


% run simulation
if matlabsolver == true

    % simbiology simulation
    [t,c,names] = sbiosimulate(Mobj);

else 

    % COPASI simulation
    if nargin >= 3
        [t,c,names] = COPASIsimulate(Mobj,solver,name);
    else
        [t,c,names] = COPASIsimulate(Mobj,solver);
    end

end

end

