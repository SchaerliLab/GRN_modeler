%% Stricker oscillator
clc
close all
clear

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))
clean_up_GRN
Ecoli = Cell('CRISPR');
solver = 'ode15s';
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-8)
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-6)
set(getconfigset(Ecoli.data.Mobj),'MaximumWallClock',10)
Ecoli = Ecoli.add_node('N1','Elowitz');
Ecoli = Ecoli.add_node('N2','Elowitz');
Ecoli = Ecoli.add_regulator('Activation_in_TF','N2','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Activation_in_TF','N1','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Repression_in_TF','N2','HILL','P_N2');
Ecoli = Ecoli.add_regulator('Repression_in_TF','N1','HILL','P_N2');
set(getconfigset(Ecoli.data.Mobj),'StopTime',1.000000e+03);
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N2'}];
set(sbioselect(Ecoli.data.Mobj,'Name','n_molecule_N2<-N1'),'Value',6.000000e+00);
set(sbioselect(Ecoli.data.Mobj,'Name','n_molecule_N1|-N2'),'Value',6.000000e+00);

%% Run simulation
Mobj = Ecoli.get_model();
[t,c,names] = run_simulation(Mobj);

% plot
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
exportgraphics(hf,'output/stricker.pdf','Resolution',300)

%% make graph
figure
Ecoli.data.inner_regulator_size = 10;
Ecoli.data.node_fontsize = 16;
Ecoli.make_graph();
axis off
exportgraphics(gcf,'output/stricker_graph.pdf','Resolution',300)