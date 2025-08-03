function [hf,exponents,pks] = calc_bifurcation_Lyapunov(Mobj,p)
%CALC_BIFURCATION_LYAPUNOV calculate the bifurcation diagram
% and the Lyapunov exponents and plot them
% Mobj: simbiology model object
% p.n_workers: number of nodes in the parallelization
% p.target_spec: name of the followed species in the bifurcation diagram
% p.t_min,p.t_end: start and end of the simulation
% p.c: concentrations for the changed variable
% p.c_name: the name of the changed variable 
% p.iszero: there is a zero Lyap exponent (true/false)

% simulation settings:
configsetobj = getconfigset(Mobj);
set(configsetobj,'SolverType','sundials','Stoptime',p.t_end)

disp('accelerate')
tic
sbioaccelerate(Mobj);
toc

%% simulate: the effect of 'PROT2', bifurcation diagram

% % Get current parallel pool
% p_pool = gcp('nocreate');
% % parallel computing
% if isempty(p_pool) || p_pool.NumWorkers ~= p.n_workers
%     delete(p_pool);
%     % start a new parallel pool
%     parpool("Processes",p.n_workers);
% end

pks = cell(length(p.c),1);

% number of the Lyapunov exponents
% % number of the non constant species
% n_exp = sum(cellfun(@(a)~a,get(Mobj.Species,'Constant')));
n_exp = 2;
exponents = zeros(length(p.c),n_exp);

tic
% parfor i = 1:numel(p.c)
for i = 1:numel(p.c)

    % for parallel computations
    Mobj_par = copyobj(Mobj);
    p_par = p;

    configset = getconfigset(Mobj_par);
    % we follow every species for Lyap exp calculaton later on
    configset.RuntimeOptions.StatesToLog = Mobj_par.Species;

    sbioaccelerate(Mobj_par)

    %  change the 'PROT2' concentration
    set(sbioselect(Mobj_par.Species,'Name',p_par.c_name),'Value',p_par.c(i))
    
    % simulate
    simdata = sbiosimulate(Mobj_par);
    [t,c,names] = getdata(simdata);
    % get rid of the beginning
    start_pos = find(t>=p_par.t_min,1,'first');

    pks{i} = findpeaks(c(start_pos:end,strcmp(names,p_par.target_spec)),'MinPeakProminence',1,'MinPeakHeight',100);

    % set the output concentrations as initial concentrations
    % to get rid of the transient next time in the ly.exp calculation
    for spec = 1:numel(Mobj_par.Species)
        set(sbioselect(Mobj_par.Species,'Name',names{spec}),'Value',c(end,spec))
    end

    % calculate largest Ly. exp.
    % start and stop time
    start_time = 1e3;
    stop_time = start_time+1e4;

    % set the start and stop time in the model
    set(get(configset,'SolverOptions'),'OutputTimes',linspace(start_time,stop_time,2));
    set(configset,'StopTime',stop_time);

    configset.SolverOptions.AbsoluteTolerance = 1e-6;
    configset.SolverOptions.RelativeTolerance = 1e-4;

    % convert model that COPASI would like the SBML
    Mobj_par = convert2irrev(Mobj_par);
    Mobj_par = correct_modifiers(Mobj_par);
    set_configset(Mobj_par,configset)

    % calulate the largest lyapunov exponent
    t_ortho = 1e2;

    exponents(i,:) = calc_lyapunov(Mobj_par,t_ortho,n_exp);

    disp([int2str(i) '/' int2str(numel(p.c))])

end
toc

%% plot the bifurcation diagram

 
hf = figure;
hold on
for i = 1:length(pks)
    scatter(p.c(i)*ones(length(pks{i}),1),pks{i},'b.')
end

if p.iszero % if there should be a zero Lyapunov exponent
    if length(p.c) > 1
        % calculate the accuracy of the Ly. exponent calculation based on
        % the standard deviation of the zero Ly exp
        zero_exp = zeros(length(p.c),1);
        for i = 1:length(p.c)
            zero_exp(i) = exponents(i,min(abs(exponents(i,:)),[],2)==abs(exponents(i,:)));
        end
        % stanard deviation
        err = std(zero_exp);
    else
        % general error value, if we do not have more data
        err = 1e-3;
    end
end


% decide whether it is chaotic
chaotic = false(length(p.c),1);
for i = 1:length(p.c)
    
    nonzero_exp = exponents(i,:);
    % set the exponent (closes to zero) to zero
    if p.iszero == true
        nonzero_exp(min(abs(exponents(i,:)),[],2)==abs(exponents(i,:))) = 0;
    end

    % % find the exponent further away from zero
    % % because here we suppose that we should have a zero Ly.exp
    % nonzero_exp = exponents(i,max(abs(exponents(i,:)),[],2)==abs(exponents(i,:)));
    chaotic(i) = max(nonzero_exp) > 3*err;
end

if numel(p.c) > 1
    % color chaotic region background
    dc = p.c(2)-p.c(1);
    % chaotic = exponents(:,1) > 1e-2;
    ax = axis;
    % color the region around the positive lyapunov exponents
    % first chaotic point in the region
    start = find(chaotic(1:end)>0,1,'first');
    while ~isempty(start)
        % the last chaotic point in the region
        stop = start-1+find(chaotic(start:end)==0,1,'first')-1;
        if isempty(stop)
            stop = length(chaotic);
        end
        start_c = p.c(1)+(start-1-0.5)*dc;
        stop_c = p.c(1)+(stop-1+0.5)*dc;
        hp = patch([start_c,stop_c,stop_c,start_c],[ax(3),ax(3),ax(4),ax(4)],'red');
        hp.FaceAlpha = 0.2; % make face a little transparent
        hp.EdgeAlpha = 0; % make edge totally transparent
        if stop == length(chaotic)
            break
        end
        start = stop+find(chaotic(stop+1:end)>0,1,'first');
    end
end

end

function set_configset(Mobj,configsetobj_in)
%SET_CONFIGSET set the configset in Mobj (simbiology model object)
%according to configsetobj_in

configsetobj_out = getconfigset(Mobj);
set(configsetobj_out,'SolverType',configsetobj_in.SolverType, ...
    'Stoptime',configsetobj_in.Stoptime,'TimeUnit',configsetobj_in.TimeUnit,...
    'MaximumWallClock',configsetobj_in.MaximumWallClock)
set(get(configsetobj_out,'CompileOptions'),...
    'Dimensionalanalysis',configsetobj_in.CompileOptions.Dimensionalanalysis,...
    'UnitConversion',configsetobj_in.CompileOptions.UnitConversion)
set(get(configsetobj_out,'SolverOptions'), ...
    'AbsoluteTolerance',configsetobj_in.SolverOptions.AbsoluteTolerance,...
    'RelativeTolerance',configsetobj_in.SolverOptions.RelativeTolerance, ...
    'OutputTimes',configsetobj_in.SolverOptions.OutputTimes)
% set(get(configsetobj_out,'RuntimeOptions'),'StatesToLog',configsetobj_in.RuntimeOptions.StatesToLog)

end