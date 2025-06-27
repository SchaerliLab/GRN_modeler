%% frequency modulation

clc
close all
clear

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

% modulation type: AM/FM
modulation = 'AM';

%%
clean_up_GRN
Ecoli = Cell('Tomazou');
solver = 'ode15s';
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-10)
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-8)
set(getconfigset(Ecoli.data.Mobj),'MaximumWallClock',60)
Ecoli = Ecoli.add_node('N1','type1');
Ecoli = Ecoli.add_node('N2','type1');
Ecoli = Ecoli.add_node('N3','type1');
Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Repression_in','N2','HILL','P_N2');
Ecoli = Ecoli.add_regulator('Repression_in','N3','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Repression_in','N3','HILL','P_N2');
Ecoli = Ecoli.add_protease('N1','PROT1','type1');
Ecoli = Ecoli.add_protease('N2','PROT2','type1');
Ecoli.set('PROT2','InitialAmount',1.000000e+01,'PROT2');
Ecoli.set('PROT1','InitialAmount',1.000000e+01,'PROT1');
set(getconfigset(Ecoli.data.Mobj),'StopTime',3e4);
set(sbioselect(Ecoli.data.Mobj,'Name','a1_N1'),'Value',1e2);%6.223776e+01);
set(sbioselect(Ecoli.data.Mobj,'Name','a1_N2'),'Value',9.5e1);%6.223776e+01);

% Tomazou circuit
% Add nodes
Ecoli = Ecoli.add_node('N4');
Ecoli = Ecoli.add_node('N5');
Ecoli = Ecoli.add_node('N6');
Ecoli = Ecoli.add_node('G');
% Add protease
Ecoli = Ecoli.add_protease('N4','C');
Ecoli = Ecoli.add_protease('N5','C');
Ecoli = Ecoli.add_protease('N6','C');
Ecoli = Ecoli.add_protease('G','L');
% Add regulators
Ecoli = Ecoli.add_regulator('Repression_in','N4','N6');
Ecoli = Ecoli.add_regulator('Repression_in','N5','N4');
Ecoli = Ecoli.add_regulator('Repression_in','N6','N5');
Ecoli = Ecoli.add_regulator('Repression_in','G','N6');
Ecoli = Ecoli.add_regulator('Activation_in','N5','Y');
Ecoli.set('Y','Constant',true,'N5','Y');
Ecoli = Ecoli.add_regulator('Activation_in','G','U');
Ecoli.set('U','Constant',true,'G','U')
switch modulation
    case 'AM'
        Ecoli = Ecoli.add_regulator('Activation_out','N5','Y','I1');
        Ecoli = Ecoli.add_regulator('Activation_in','G','U','P_N3');
    case 'FM'
        Ecoli = Ecoli.add_regulator('Activation_in','N5','Y','P_N3');
        Ecoli = Ecoli.add_regulator('Activation_out','G','U','I1');
end
% Initial conentrations
Ecoli.set('P_N5','InitialAmount',100,'N5')
Ecoli.set('C','InitialAmount',100,'C')
Ecoli.set('U','Value',100,'G','U')
Ecoli.set('Y','Value',100,'N5','Y')
Ecoli.set('L','InitialAmount',300,'L')
% Parameters
% G copy number 
Ecoli.set('n_copy_G','Value',50,'G')
Ecoli.set('a1_G','Value',1000,'G')
Ecoli.set('K_molecule_G<-U','Value',50,'G','U')
Ecoli.set('n_copy_N5','Value',50,'N5')
Ecoli.set('a1_N5','Value',10000,'N5')
Ecoli.set('K_molecule_N5<-Y','Value',50,'N5','Y')

Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_G','P_N4','P_N3','P_N1'}];

Mobj = Ecoli.get_model();

set(sbioselect(Mobj.Species,'Name','I1'),'InitialAmount',30); % uM
% set(sbioselect(Mobj.Species,'Name','I2'),'InitialAmount',100); % uM

%% change the rule for the N3 -| N4 repression
% we need just a mild/partial repression here

% % the repression will be effective just in 'ratio' ratio:
% addparameter(Mobj,'ratio',1,'Units','dimensionless');
% 
% regulation = sbioselect(Mobj.Rule,'Name','Rule_HILL_N5<-Y<-N3');
% regulation.Rule = [regulation.Rule '*ratio + (1-ratio)'];

% higher K
switch modulation
    case 'AM'
        set(sbioselect(Mobj.Parameters,'Name','K_molecule_G<-U<-N3'),'Value',5e2)
    case 'FM'
        set(sbioselect(Mobj.Parameters,'Name','K_molecule_N5<-Y<-N3'),'Value',1e3)
end

% regulation = sbioselect(Mobj.Rule,'Name','Rule_HILL_N4|-N1');
% regulation.Rule = [regulation.Rule '*ratio + (1-ratio)'];
% regulation = sbioselect(Mobj.Rule,'Name','Rule_HILL_N4|-N2');
% regulation.Rule = [regulation.Rule '*ratio + (1-ratio)'];
%% run simulation

% simulate
tic
[t,c,names] = sbiosimulate(Mobj);
toc
%% plot

hf = figure;
plot(t,c(:,strcmp(names,'P_G')),'LineWidth',2)
set(gca,'XLim',[1.9e4,2.7e4])
ylim = get(gca,'YLim');
set(gca,'YLim',[0,ylim(2)])

xlabel('\bf\boldmath$t$ / minute','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$P_{G}$ / molecule','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
box off
exportgraphics(hf,['output/' modulation '_modul_indep.pdf'],'Resolution',300)