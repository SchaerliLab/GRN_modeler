%% calculate the time period of Tomazou repressilator in the function of nodes
clc
clear
close all

% we evaluate the simulations only after this time / min
start_time = 1e3;
stop_time = 1e4;

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

% number of nodes
nodes_repressilator = 3:2:20;
% Time period / min
T = zeros(numel(nodes_repressilator,2));

%% repressilator

tic
for i = 1:numel(nodes_repressilator)
    % build the model
    Ecoli = build_repressilator('Elowitz',nodes_repressilator(i));
    
    T(i,1) = set_run(Ecoli,start_time,stop_time);

    disp([int2str(i) '/' int2str(numel(nodes_repressilator))])
end
toc

%% repressilator inner 3 node version

tic
for i = 1:numel(nodes_repressilator)
    if i > 7
        % build the model
        Ecoli = build_repressilator('Elowitz',nodes_repressilator(i));

        % add crossing to the model
        Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL','P_N7');
        
        T(i,2) = set_run(Ecoli,start_time,stop_time);
    else
        T(i,2) = nan;
    end

    disp([int2str(i) '/' int2str(numel(nodes_repressilator))])
end
toc

%% plot the data

!mkdir -p output

% fit line
fitobject = fit(nodes_repressilator(2:end).',T(2:end,1),'poly1');

save('output/repr_line','fitobject','nodes_repressilator','T')

% time period
figure
hold on
plot(nodes_repressilator,T(:,1),'LineWidth',2,'LineStyle','--','Marker','o')
plot(nodes_repressilator,T(:,2),'LineWidth',2,'LineStyle','--','Marker','o')
% plot(fitobject)
xlabel('\bf number of nodes','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$T$ / min','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(gcf,'output/T_repressilator.pdf')

%% set and run simulations
function T = set_run(Ecoli,start_time,stop_time)

% simulation settings
Ecoli.set('P_N1','InitialAmount',1.000000e+02,'N1');
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-8)
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-6)
set(getconfigset(Ecoli.data.Mobj),'Stoptime',stop_time);
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];

% run simulation
Mobj = Ecoli.get_model();
[t,c,~] = run_simulation(Mobj,'ode15s');

% get rid of the initial transient part
t_min = find(t>start_time,1,'first');
t(1:t_min) = [];
c(1:t_min,:) = [];

% calculate time period
T = calc_timeperiod(t,c);
end