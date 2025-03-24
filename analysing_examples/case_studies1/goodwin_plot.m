%% Goodwin oscillator
clc
close all
clear

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

%% build the model
clean_up_GRN
Ecoli = Cell('Tomazou');
Ecoli = Ecoli.add_node('N1','type1');
Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL','P_N1');
Ecoli = Ecoli.add_protease('N1','PROT1','type1');
Ecoli.set('PROT1','InitialAmount',5.000000e+01,'PROT1');
set(getconfigset(Ecoli.data.Mobj),'Stoptime',1e2);
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-10);
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-8);
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'uP_N1'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'mRNA_N1'}];
Ecoli.data.StatesToLog(strcmp(Ecoli.data.StatesToLog,'mRNA_N1')) = [];
Ecoli.data.StatesToLog(strcmp(Ecoli.data.StatesToLog,'uP_N1')) = [];
Ecoli.data.StatesToLog(strcmp(Ecoli.data.StatesToLog,'P_N1')) = [];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'mRNA_N1'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'uP_N1'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];

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
exportgraphics(hf,'output/goodwin.pdf','Resolution',300)

%% make graph
figure
Ecoli.data.inner_regulator_size = 10;
Ecoli.data.node_fontsize = 16;
Ecoli.make_graph();
axis off
exportgraphics(gcf,'output/goodwin_graph.pdf','Resolution',300)