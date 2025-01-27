%% Build mixed repressilators when certain edges are replaced with activations
clc
clear
close all

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

%% build the model

% number of nodes
N = 6;
% edges for activation
A = [1,2,3];%[1,3,5];
% Hill exponent
n = 6;

% generated name
name = ['N', int2str(N), '_A', strrep(int2str(A),' ','')];

% build the circuit
Ecoli = build_mixed_repressilator('Elowitz',N,A);

% set the Hil exponent
set(sbioselect(Ecoli.data.Mobj,'Name','n'),'Value',n);

% plot the circuit
Ecoli.make_graph();
box off
axis off
exportgraphics(gcf,['output/mixed_circuit_', name, '.pdf'],'Resolution',300)

%% run the simulation

% simulation settings
set(getconfigset(Ecoli.data.Mobj),'StopTime',1e3);
Ecoli.set('P_N1','InitialAmount',1.000000e+02,'N1');

solver = 'ode15s';
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-8)
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-6)

% run and plot
Mobj = Ecoli.get_model();
set(getconfigset(Mobj),'SolverType',solver);
[t,c,names] = run_simulation(Mobj,solver);
hf = figure;
hold on
for i = 1:numel(Ecoli.data.StatesToLog)
    plot(t,c(:,strcmp(names,Ecoli.data.StatesToLog{i})),'LineWidth',2)
end
xlabel('\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$c$ / molecule','interpreter','latex','Fontsize',18)
legend(Ecoli.data.StatesToLog,'Interpreter','none')
set(gca,'LineWidth',2,'Fontsize',16)

% save plot
exportgraphics(hf,['output/mixed_', name, '_n', int2str(n), '.pdf'],'Resolution',300)