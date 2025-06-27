%% frequency modulation
% repressilator (-|) goodwin
% unidirectionally coupled

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
Ecoli = Ecoli.add_node('N2','type1');
Ecoli = Ecoli.add_node('N3','type1');
Ecoli = Ecoli.add_node('N4','type1');
Ecoli = Ecoli.add_regulator('Repression_in','N2','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Repression_in','N3','HILL','P_N2');
Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL','P_N3');
Ecoli = Ecoli.add_regulator('Repression_in','N4','HILL','P_N4');
Ecoli = Ecoli.add_regulator('Repression_in','N4','HILL','P_N3');
Ecoli = Ecoli.add_protease('N1','PROT1','type1');
Ecoli = Ecoli.add_protease('N2','PROT1','type1');
Ecoli = Ecoli.add_protease('N3','PROT1','type1');
Ecoli = Ecoli.add_protease('N4','PROT2','type1');
Ecoli.data.Accelerate = true;
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1.000000e-8);
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1.000000e-10);
set(getconfigset(Ecoli.data.Mobj),'Stoptime',1.000000e+03);
set(getconfigset(Ecoli.data.Mobj),'SolverType','sundials');
Ecoli.set('P_N1','InitialAmount',1.000000e+03,'N1');
Ecoli.set('P_N4','InitialAmount',1.000000e+02,'N4');
Ecoli.set('PROT1','InitialAmount',1.000000e+02,'PROT1');
Ecoli.set('PROT2','InitialAmount',1.000000e+01,'PROT2');
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N4','P_N3'}];

% Ecoli = Ecoli.add_node('N5','type1');
% Ecoli = Ecoli.add_regulator('Repression_in','N5','HILL','P_N4');
% Ecoli = Ecoli.add_protease('N5','PROT3','type1');
% Ecoli.set('PROT3','InitialAmount',3.000000e+01,'PROT3');

Mobj = Ecoli.get_model();

%% change the rule for the N3 -| N4 repression
% we need just a mild/partial repression here

% the repression will be effective just in 'ratio' ratio:
addparameter(Mobj,'ratio',0.9,'Units','dimensionless');

regulation = sbioselect(Mobj.Rule,'Name','Rule_HILL_N4|-N3');
regulation.Rule = [regulation.Rule '*ratio + (1-ratio)'];

% higher K
set(sbioselect(Mobj.Parameters,'Name','K_molecule_N4|-N3'),'Value',50)
% set(sbioselect(Mobj.Parameters,'Name','n_molecule_N4|-N3'),'Value',4)

% %higher n in readout node
% set(sbioselect(Mobj.Parameters,'Name','K_molecule_N5|-N4'),'Value',2)
% set(sbioselect(Mobj.Parameters,'Name','n_molecule_N5|-N4'),'Value',4)
%% run simulation

% simulate
[t,c,names] = sbiosimulate(Mobj);

%% plot

% minimum and maximum time for plot
tmin = 425;
tmax = 750;

% time position in the data
pos1 = find(t>=tmin,1,'first');
pos2 = find(t<=tmax,1,'last');

% get rid of the excess
t = t(pos1:pos2);
c = c(pos1:pos2,:);

hf = figure;
% hold on
% yyaxis left
plot(t,c(:,strcmp(names,'P_N4')),'LineWidth',2)
% yyaxis right
% plot(t,c(:,strcmp(names,'P_N3')))
set(gca,'Xlim',[tmin,tmax])
set(gca,'LineWidth',2,'Fontsize',16)
xlabel('\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$\left[P_{N4}\right]$ / molecule','interpreter','latex','Fontsize',18)
exportgraphics(hf,'output/FM_time.pdf','Resolution',300)

% find the peaks
[pks_max,pos] = findpeaks(c(:,strcmp(names,'P_N4')));
[pks_min,~] = findpeaks(-c(:,strcmp(names,'P_N4')));
pks_min = -pks_min;

% Time period
hf = figure;
plot(diff(t(pos)),'LineWidth',2)
set(gca,'LineWidth',2,'Fontsize',16)
xlabel('\bf\boldmath$$\#$peaks','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$T$ / min','interpreter','latex','Fontsize',18)
exportgraphics(hf,'output/FM_T.pdf','Resolution',300)

% we should have the same number of maximums and minimums
if numel(pks_max) < numel(pks_min)
    pks_min(end) = [];
elseif numel(pks_max) > numel(pks_min)
    pks_max(end) = [];
end
% Amplitude
hf = figure;
plot(pks_max-pks_min,'LineWidth',2)
set(gca,'LineWidth',2,'Fontsize',16)
xlabel('\bf\boldmath$$\#$peaks','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$A$ / molecule','interpreter','latex','Fontsize',18)
exportgraphics(hf,'output/FM_A.pdf','Resolution',300)