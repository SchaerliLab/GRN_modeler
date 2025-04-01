%% create graphs for different regulation systems
clc
clear
close all

!mkdir -p output

%% act-act
clean_up_GRN
Ecoli = Cell('Elowitz');
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-8)
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-6)
set(getconfigset(Ecoli.data.Mobj),'MaximumWallClock',10)
Ecoli = Ecoli.add_node('N1','type1');
Ecoli = Ecoli.add_node('N2','type1');
Ecoli = Ecoli.add_node('N3','type1');
Ecoli = Ecoli.add_regulator('Activation_in','N1','HILL','P_N2');
Ecoli = Ecoli.add_regulator('Activation_in','N1','HILL','P_N3');
Ecoli.make_graph();
exportgraphics(gcf,'output/act_act.pdf','Resolution',300)

%% act-rep
clean_up_GRN
Ecoli = Cell('Elowitz');
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-8)
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-6)
set(getconfigset(Ecoli.data.Mobj),'MaximumWallClock',10)
Ecoli = Ecoli.add_node('N1','type1');
Ecoli = Ecoli.add_node('N2','type1');
Ecoli = Ecoli.add_node('N3','type1');
Ecoli = Ecoli.add_regulator('Activation_in','N1','HILL','P_N2');
Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL','P_N3');
Ecoli.make_graph();
exportgraphics(gcf,'output/act_rep.pdf','Resolution',300)
%% rep-rep
clean_up_GRN
Ecoli = Cell('Elowitz');
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-8)
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-6)
set(getconfigset(Ecoli.data.Mobj),'MaximumWallClock',10)
Ecoli = Ecoli.add_node('N1','type1');
Ecoli = Ecoli.add_node('N2','type1');
Ecoli = Ecoli.add_node('N3','type1');
Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL','P_N3');
Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL','P_N2');
Ecoli.make_graph();
exportgraphics(gcf,'output/rep_rep.pdf','Resolution',300)
%% rep-rep-act
clean_up_GRN
Ecoli = Cell('Elowitz');
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-8)
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-6)
set(getconfigset(Ecoli.data.Mobj),'MaximumWallClock',10)
Ecoli = Ecoli.add_node('N1','type1');
Ecoli = Ecoli.add_node('N2','type1');
Ecoli = Ecoli.add_node('N3','type1');
Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL','P_N3');
Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL','P_N2');
Ecoli = Ecoli.add_node('N4','type1');
Ecoli = Ecoli.add_regulator('Activation_in','N1','HILL','P_N4');
Ecoli.make_graph();
exportgraphics(gcf,'output/rep_rep_act.pdf','Resolution',300)