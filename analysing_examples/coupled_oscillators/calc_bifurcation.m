function [hf,pks] = calc_bifurcation(Mobj,p)
%CALC_BIFURCATION calculate the bifurcation diagram and plot it
% Mobj: simbiology model object
% p.n_workers: number of nodes in the parallelization
% p.target_spec: name of the followed species in the bifurcation diagram
% p.t_min,p.t_end: start and end of the simulation
% p.c: concentrations for the changed variable
% p.c_name: the name of the changed variable 

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

tic
% parfor i = 1:numel(p.c)
for i = 1:numel(p.c)

    % % for parallel computations
    % Mobj_par = copyobj(Mobj);
    % p_par = p;

    % sbioaccelerate(Mobj_par)

    %  change the given parameter concentration
    set(sbioselect(Mobj,'Name',p.c_name),'Value',p.c(i))
    
    % simulate
    simdata = sbiosimulate(Mobj);
    [t,c,names] = getdata(simdata);
    % get rid of the beginning
    start_pos = find(t>=p.t_min,1,'first');

    pks{i} = findpeaks(c(start_pos:end,strcmp(names,p.target_spec)),'MinPeakProminence',1,'MinPeakHeight',100);

end
toc

%% plot the bifurcation diagram

hf = figure;
hold on
for i = 1:length(pks)
    scatter(p.c(i)*ones(length(pks{i}),1),pks{i},'b.')
end

end
