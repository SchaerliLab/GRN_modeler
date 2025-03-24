%% Toggle switch
clc
close all
clear

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

%% build the model
clean_up_GRN
Ecoli = Cell('Elowitz');
solver = 'ode15s';
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-8)
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-6)
set(getconfigset(Ecoli.data.Mobj),'MaximumWallClock',10)
Ecoli = Ecoli.add_node('N1','type1');
Ecoli = Ecoli.add_node('N2','type1');
Ecoli = Ecoli.add_regulator('Repression_in','N2','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL','P_N2');
set(getconfigset(Ecoli.data.Mobj),'StopTime',1e3);
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N2'}];
Ecoli.set('P_N2','InitialAmount',1.000000e+02,'N2');
Ecoli.set('P_N1','InitialAmount',1.000000e+02,'N1');

Mobj = Ecoli.get_model();

% speed up simulations
sbioaccelerate(Mobj)

%% Run simulation

% change the initial condition for the P_N1 protein
set(sbioselect(Mobj.Species,'InitialAmount','P_N1'),'Value',90)

% simulation
[t,c,names] = run_simulation(Mobj);

hf = figure;
hold on
for i = 1:numel(Ecoli.data.StatesToLog)
    plot(t,c(:,strcmp(names,Ecoli.data.StatesToLog{i})),'LineWidth',2)
end
xlabel('\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$c$ / molecule','interpreter','latex','Fontsize',18)
legend(Ecoli.data.StatesToLog,'Interpreter','none')
set(gca,'LineWidth',2,'Fontsize',16)
!mkdir -p output
exportgraphics(hf,'output/toggleswitch.pdf','Resolution',300)

%% parameter scan to show bistability

% number of steps
n = 100;
% minimal and maximal protein concentration
c_min = 10;
c_max = 2*100-c_min;

% tolerance values for finding the steady state
AbsTol = 1e-8; % for algebraic method
RelTol = 1e-6; % for simulation method
Method = 'auto'; % 'auto', 'simulation' or 'algebraic'

% initial condition
c_in = linspace(c_min,c_max,n);

% store the output concentration
c_out = nan(n,2);

% calculate steady states
tic
for i = 1:n
    % change the initial condition for the P_N1 protein
    set(sbioselect(Mobj.Species,'Name','P_N1'),'InitialAmount',c_in(i))
    % calculate steady state
    [success, ~, Mobj_out] = sbiosteadystate(Mobj, ...
        'Method',Method,'AbsTol',AbsTol,'RelTol',RelTol);
    % if the simulation was successful
    if success == 1
        c_out(i,1) = get(sbioselect(Mobj_out.Species,'Name','P_N1'),'InitialAmount');
        c_out(i,2) = get(sbioselect(Mobj_out.Species,'Name','P_N2'),'InitialAmount');
    end
end
toc

%% plot the steady states
hf = figure;
plot(c_in,c_out,'LineWidth',2)
xlabel('\bf\boldmath$\left[\textrm{P}_{N1}\right]_0$ / molecule','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$c$ / molecule','interpreter','latex','Fontsize',18)
legend(Ecoli.data.StatesToLog,'Interpreter','none','Location','best')
set(gca,'LineWidth',2,'Fontsize',16)
!mkdir -p output
exportgraphics(hf,'output/toggleswitch_steadystate.pdf','Resolution',300)

%% make graph
figure
Ecoli.data.inner_regulator_size = 10;
Ecoli.data.node_fontsize = 16;
Ecoli.make_graph();
axis off
exportgraphics(gcf,'output/toggle_graph.pdf','Resolution',300)