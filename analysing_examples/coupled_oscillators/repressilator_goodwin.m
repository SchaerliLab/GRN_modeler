%% repressilator |- goodwin
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
Ecoli = Ecoli.add_regulator('Repression_in','N3','HILL','P_N4');
Ecoli = Ecoli.add_protease('N1','PROT1','type1');
Ecoli = Ecoli.add_protease('N2','PROT1','type1');
Ecoli = Ecoli.add_protease('N3','PROT1','type1');
Ecoli = Ecoli.add_protease('N4','PROT2','type1');
Ecoli.data.Accelerate = true;
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1.000000e-8);
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1.000000e-10);
set(getconfigset(Ecoli.data.Mobj),'Stoptime',1.000000e+04);
set(getconfigset(Ecoli.data.Mobj),'SolverType','sundials');
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];
Ecoli.set('P_N1','InitialAmount',1.000000e+03,'N1');
Ecoli.set('P_N4','InitialAmount',1.000000e+03,'N4');
Ecoli.set('PROT1','InitialAmount',1.000000e+02,'PROT1');
Ecoli.set('PROT2','InitialAmount',1.000000e+02,'PROT2');

%% calculate bifurcation diagram and Lyapunov exponents

% get the model
Mobj = Ecoli.get_model();

% settings
p.n_workers = 1; % number of nodes in the parallelization
p.target_spec = 'P_N1'; % name of the followed species in the bifurcation diagram
% start and end of the simulation
p.t_min = 1e5;
p.t_end = 1e5+p.t_min;
p.c = linspace(1,20,200); % concentrations for the changed variable
p.c_name = 'PROT2'; % the name of the changed variable 
p.iszero = true;

% calculate bif diagr and Ly exps
[hf,exponents,pks] = calc_bifurcation_Lyapunov(Mobj,p);

% % plot should be the last one (at the top of other things)
% set(gca,'children',flipud(get(gca,'children')))
% set(gca, "Layer", "top")
% xlim([p.c(1),p.c(end)])
xlabel('\bf\boldmath$\left[\textrm{PROT}_2\right]$ / molecule','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$\max\left[(P_{N1})\right]$ / molecule','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(hf,'output/repressilator_goodwin_bifurcation_Lyap.pdf','Resolution',300)
save('output/repressilator_goodwin_bifurcation_Lyap.mat','exponents','pks')

%% Poincare map at a given point

for prot2 = [3,10]

% get the model
Mobj = Ecoli.get_model();
sbioaccelerate(Mobj)

% change the 'PROT2' concentration
% prot2 = 10; % 3 10
set(sbioselect(Mobj.Species,'Name','PROT2'),'Value',prot2)

% target concentartion for Poincare section
p.c_target = 1e4;
% name of the target species
p.target_name = 'uP_N1';
% transient time
p.t_min = 1e5;
% in inverse interpolation use n_interp in both directions
p.n_interp = 3;

% end of the simulation
t_end = p.t_min+1e5;

% simulation settings:
configsetobj = getconfigset(Mobj);
set(get(configsetobj,'RuntimeOptions'),'StatesToLog',{'P_N1','uP_N1','mRNA_N1'})
set(configsetobj,'SolverType','sundials','Stoptime',t_end)

% simulate
tic
[t,c,names,t_cross,pos_min] = crossing_times(Mobj,p);
toc

% get the concentration values from the data series in the given time
% points
c_cross_P_N1 = interp1(t,c(:,strcmp(names,'P_N1')),t_cross);
c_cross_mRNA_N1 = interp1(t,c(:,strcmp(names,'mRNA_N1')),t_cross);

% Poincare section
hf = figure;
scatter(c_cross_mRNA_N1,c_cross_P_N1,'blue','filled')
xlabel('\bf\boldmath$\left[mRNA_{N1}\right]$ / molecule','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$\left[P_{N1}\right]$ / molecule','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(hf,['output/repressilator_goodwin_poincare_section_' num2str(prot2) '.pdf'],'Resolution',300)

% Poincare map
hf = figure;
scatter(c_cross_P_N1(1:end-1),c_cross_P_N1(2:end),'blue','filled')
xlabel('\bf\boldmath$\left[P_{N1}\right]$(n) / molecule','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$\left[P_{N1}\right]$(n+1) / molecule','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(hf,['output/repressilator_goodwin_poincare_map_' num2str(prot2) '.pdf'],'Resolution',300)

% phase space of a node
hf = figure;
plot3(c(pos_min:end,strcmp(names,'P_N1')),c(pos_min:end,strcmp(names,'mRNA_N1')),c(pos_min:end,strcmp(names,'uP_N1')))
xlabel('\bf\boldmath$\left[P_{N1}\right]$','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$\left[mRNA_{N1}\right]$','interpreter','latex','Fontsize',18)
zlabel('\bf\boldmath$\left[uP_{N1}\right]$','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
box on
exportgraphics(hf,['output/repressilator_goodwin_phasespace_' num2str(prot2) '.pdf'],'Resolution',300)

% time series
hf = figure;
plot(t,c(:,strcmp(names,'P_N1')))
set(gca,'XLim',[0,1e4])
xlabel('\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$\left[P_{N1}\right]$ / molecule','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(hf,['output/repressilator_goodwin_timeseries_' num2str(prot2) '.pdf'],'Resolution',300)

end