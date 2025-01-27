%% actolator2
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
Ecoli = Ecoli.add_regulator('Mixed_in','N2','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Mixed_in','N1','HILL','P_N2');
set(getconfigset(Ecoli.data.Mobj),'StopTime',1e3);
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N2'}];
Ecoli.set('P_N1','InitialAmount',1.000000e+02,'N1');
set(sbioselect(Ecoli.data.Mobj,'Name','n_N2<-N1'),'Value',6.000000e+00);
set(sbioselect(Ecoli.data.Mobj,'Name','n_N1<-N2'),'Value',6.000000e+00);
%% Run simulation
Mobj = Ecoli.get_model();
set(getconfigset(Mobj),'SolverType',solver);
[t,c,names] = run_simulation(Mobj,solver);
figure
hold on
for i = 1:numel(Ecoli.data.StatesToLog)
    plot(t,c(:,strcmp(names,Ecoli.data.StatesToLog{i})),'LineWidth',2)
end
xlabel('\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$c$ / molecule','interpreter','latex','Fontsize',18)
legend(Ecoli.data.StatesToLog,'Interpreter','none')
set(gca,'LineWidth',2,'Fontsize',16)
!mkdir -p output
exportgraphics(gcf,'output/actolator2.pdf','Resolution',300)
