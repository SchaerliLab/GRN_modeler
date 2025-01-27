%% reptolator6
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
set(getconfigset(Ecoli.data.Mobj),'StopTime',2.000000e+03);
Ecoli = Ecoli.add_node('N1','type1');
Ecoli = Ecoli.add_node('N2','type1');
Ecoli = Ecoli.add_node('N3','type1');
Ecoli = Ecoli.add_node('N4','type1');
Ecoli = Ecoli.add_node('N5','type1');
Ecoli = Ecoli.add_node('N6','type1');
Ecoli = Ecoli.add_regulator('Repression_in','N2','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Repression_in','N3','HILL','P_N2');
Ecoli = Ecoli.add_regulator('Repression_in','N4','HILL','P_N3');
Ecoli = Ecoli.add_regulator('Repression_in','N5','HILL','P_N4');
Ecoli = Ecoli.add_regulator('Repression_in','N6','HILL','P_N5');
Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL','P_N6');
Ecoli = Ecoli.add_regulator('Repression_in','N4','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL','P_N4');
Ecoli = Ecoli.add_regulator('Repression_in','N5','HILL','P_N2');
Ecoli = Ecoli.add_regulator('Repression_in','N2','HILL','P_N5');
Ecoli = Ecoli.add_regulator('Repression_in','N6','HILL','P_N3');
Ecoli = Ecoli.add_regulator('Repression_in','N3','HILL','P_N6');
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N2'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N3'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N4'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N5'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N6'}];
Ecoli.set('P_N1','InitialAmount',1.000000e+03,'N1');

%% Run simulation
Mobj = Ecoli.get_model();
configset = getconfigset(Mobj);
Mobj = convert2irrev(Mobj);
Mobj = correct_modifiers(Mobj);
Ecoli.set_configset(Mobj,configset)
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
exportgraphics(gcf,'output/reptolator6.pdf','Resolution',300)
