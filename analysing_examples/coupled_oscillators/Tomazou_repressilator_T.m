%% calculate the time period of Tomazou repressilator in the function of nodes
clc
clear
close all

% we evaluate the simulations only after this time / min
start_time = 1e3;
stop_time = 1e4;

%% repressilator

% number of nodes
nodes_repressilator = 3:2:20;
% Time period / min
T = zeros(numel(nodes_repressilator,1));

tic
for i = 1:numel(nodes_repressilator)
    % build the model
    Ecoli = build_repressilator('Tomazou',nodes_repressilator(i));
    
    % simulation settings
    Ecoli.set('P_N1','InitialAmount',1.000000e+03,'N1');
    set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-8)
    set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-6)
    set(getconfigset(Ecoli.data.Mobj),'Stoptime',stop_time);
    Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];
    % change the protease conc. so we have the same amount/node
    Ecoli.set('PROT1','InitialAmount',100,'PROT1');

    % run simulation
    Mobj = Ecoli.get_model();
    [t,c,~] = run_simulation(Mobj,'ode15s');

    % get rid of the initial transient part
    t_min = find(t>start_time,1,'first');
    t(1:t_min) = [];
    c(1:t_min,:) = [];

    % store the time period
    thres = (max(c)-min(c))/2;
    [~,dt_high,dt_low] = calc_duty(t,c,thres);
    T(i,1) = dt_high+dt_low;
    % T(i,1) = calc_timeperiod(t,c);

    disp([int2str(i) '/' int2str(numel(nodes_repressilator))])
end
toc

%% plot the data

!mkdir -p output

% time period
figure
hold on
plot(nodes_repressilator,T,'LineWidth',2,'LineStyle','--','Marker','o')
xlabel('\bf number of nodes','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$T$ / min','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(gcf,'output/T_compare.pdf')