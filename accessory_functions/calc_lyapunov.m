function exponents = calc_lyapunov(Mobj,t_ortho,n_exp)
%CALC_LYAPUNOV calculate Lyapunov exponents with COPASI
% Mobj: simbiology model object, we will use the settings from this
% t_ortho: Orthonormalization interval
% n_exp: number of lyapunov exponents
% (https://copasi.org/Support/User_Manual/Methods/Lyapunov_Exponents_Calculation/)

problem = py.dict();
method = py.dict();

% name to export the model for COPASI
model_name = [tempdir 'model.sbml'];

% export SBML file into the temporary directory of the system
sbmlexport(Mobj, model_name);

% get simulation settings
configsetobj = getconfigset(Mobj);

% tolerance values
method{"Absolute Tolerance"} = configsetobj.SolverOptions.AbsoluteTolerance;
method{"Relative Tolerance"} = configsetobj.SolverOptions.RelativeTolerance;

% maximum number of internal steps
method{"Max Internal Steps"} = 1e4;

% number of independent variables = number of lyapunov exponents
% (the number of species that are not constant minus the number of mass conservation relations)
% the number of calculatd ly.exps
if nargin >= 3
    problem{"ExponentNumber"} = n_exp;
else
    problem{"ExponentNumber"} = 1;
end

% we doo not nees divergence calculation
problem{"DivergenceRequested"} = "False";

% orthogonalization interval
if nargin >= 2 && exist('t_ortho','var')
    method{"Orthonormalization Interval"} = t_ortho;
else
    method{"Orthonormalization Interval"} = 1;
end
% transient time
if ~isempty(configsetobj.SolverOptions.OutputTimes)
    problem{"TransientTime"} = configsetobj.SolverOptions.OutputTimes(1);
else
    problem{"TransientTime"} = 0;
end
% simulation times
method{"Overall time"} = configsetobj.StopTime;

exponents = pyrunfile("lyapunov_COPASI.py","exponents",model_name=model_name,problem=problem,method=method);

end

