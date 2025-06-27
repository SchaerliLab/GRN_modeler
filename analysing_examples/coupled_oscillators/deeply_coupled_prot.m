%% deeply coupled oscillators
% 4 node circuit (containing two repressilator) with one protease

clc
clear
close all

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

%% create the model (examples/node4_deeply_coupled.mat)

clean_up_GRN
Ecoli = Cell('Tomazou');
Ecoli = Ecoli.add_node('N1','type1');
Ecoli = Ecoli.add_node('N2','type1');
Ecoli = Ecoli.add_node('N3','type1');
Ecoli = Ecoli.add_node('N4','type1');
Ecoli = Ecoli.add_regulator('Repression_in','N2','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Repression_in','N3','HILL','P_N2');
Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL','P_N3');
Ecoli = Ecoli.add_regulator('Repression_in','N4','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Repression_in','N3','HILL','P_N4');
Ecoli = Ecoli.add_protease('N1','PROT1','type1');
Ecoli = Ecoli.add_protease('N2','PROT1','type1');
Ecoli = Ecoli.add_protease('N3','PROT1','type1');
Ecoli = Ecoli.add_protease('N4','PROT1','type1');
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-10);
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-8);
set(getconfigset(Ecoli.data.Mobj),'Stoptime',1e4);
Ecoli.data.Accelerate = true;
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N2'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N3'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N4'}];
Ecoli.set('P_N1','InitialAmount',1.000000e+03,'N1');
Ecoli.set('uP_N1','InitialAmount',1.000000e+03,'N1');
% Ecoli.set('PROT1','InitialAmount',3.120000e+02,'PROT1');
% Ecoli = Ecoli.add_protease('N4','PROT2','type1');
% set(sbioselect(Ecoli.data.Mobj,'Name','PROT2'),'Value',1.207692e+02);

%%

% get the model
Mobj = Ecoli.get_model();

% settings
p.n_workers = 1; % number of nodes in the parallelization
p.target_spec = 'P_N4'; % name of the followed species in the bifurcation diagram
% start and end of the simulation
p.t_min = 1e5;
p.t_end = 1e5+p.t_min;
p.c = linspace(10,500,491); % concentrations for the changed variable
p.c_name = 'PROT1'; % the name of the changed variable 
p.iszero = true;

% calculate bif diagr and Ly exps
[hf,exponents,pks] = calc_bifurcation_Lyapunov(Mobj,p);

%% plot the bifurcation diagram

% hf = figure;
% hold on
% for i = 1:length(pks)
%     scatter(c_prot(i)*ones(length(pks{i}),1),pks{i},'b.')
% end
% % plot should be the last one (at the top of other things)
% set(gca,'children',flipud(get(gca,'children')))
% set(gca, "Layer", "top")
xlabel('\bf\boldmath$\left[\textrm{PROT}_1\right]$ / molecule','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$\max\left[(P_{N4})\right]$ / molecule','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(hf,'output/deeply_coupled_bifurcation_Lyap.pdf','Resolution',300)
save('output/deeply_coupled_bifurcation_Lyap.mat','exponents','pks')

%% Poincare map at a given point

! mkdir -p output

prot1 = 310; % period-3: 290
% change the 'PROT1' concentration
set(sbioselect(Mobj.Species,'Name','PROT1'),'Value',prot1)

% start and end of the simulation
t_min = 1e5;
t_end = 1e5+t_min;


% target concentartion for Poincare section
c_target = 1e2;

% in inverse interpolation use n_interp in both directions
n_interp = 3;

% simulation settings:
configsetobj = getconfigset(Mobj);
set(get(configsetobj,'RuntimeOptions'),'StatesToLog',{'P_N4','uP_N4','mRNA_N4'})
set(configsetobj,'SolverType','sundials','Stoptime',t_end)

% simulate
simdata = sbiosimulate(Mobj);
[t,c,names] = getdata(simdata);


% find the points where uP_N4 crosses c_target from the positive direction
y = c(:,strcmp(names,'uP_N4')) - c_target;
dc = y(1:end-1).*y(2:end);
% find where dc is negative (uP_N4 crosses c_target)
pos = 1:length(dc);
pos = pos(dc<0);
% get rid of the cases when uP_N4 started under c_target
pos(y(pos)<0) = [];
% we are not using the time points before t_min
pos_min = find(t>t_min,1,'first');
pos(pos<pos_min) = [];

% we need some extra data point for the interpolation
if length(t) < pos(end)+n_interp
    pos(end) = [];
end

% interpolate to find the appropriate time values
t_cross = zeros(length(pos),1);
for i = 1:numel(pos)
    range = pos(i)-n_interp:pos(i)+n_interp;
    t_cross(i) = interp1(y(range),t(range),0);
end

% get the concentration values from the data series in the given time
% points
c_cross_P_N4 = interp1(t,c(:,strcmp(names,'P_N4')),t_cross);
c_cross_mRNA_N4 = interp1(t,c(:,strcmp(names,'mRNA_N4')),t_cross);

% Poincare section
hf = figure;
scatter(c_cross_mRNA_N4,c_cross_P_N4,'blue','filled')
xlabel('\bf\boldmath$\left[mRNA_{N4}\right]$ / molecule','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$\left[P_{N4}\right]$ / molecule','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(hf,['output/deeply_coupled_poincare_section_' num2str(prot1) '.pdf'],'Resolution',300)

% Poincare map
hf = figure;
scatter(c_cross_P_N4(1:end-1),c_cross_P_N4(2:end),'blue','filled')
xlabel('\bf\boldmath$\left[P_{N4}\right]$(n) / molecule','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$\left[P_{N4}\right]$(n+1) / molecule','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(hf,['output/deeply_coupled_poincare_map_' num2str(prot1) '.pdf'],'Resolution',300)

% phase space of a node
hf = figure;
plot3(c(pos_min:end,strcmp(names,'P_N4')),c(pos_min:end,strcmp(names,'mRNA_N4')),c(pos_min:end,strcmp(names,'uP_N4')))
xlabel('\bf\boldmath$\left[P_{N4}\right]$','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$\left[mRNA_{N4}\right]$','interpreter','latex','Fontsize',18)
zlabel('\bf\boldmath$\left[uP_{N4}\right]$','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
box on
exportgraphics(hf,['output/deeply_coupled_phasespace_' num2str(prot1) '.pdf'],'Resolution',300)

% time series
hf = figure;
plot(t,c(:,strcmp(names,'P_N4')))
set(gca,'XLim',[0,1e4])
xlabel('\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$\left[P_{N4}\right]$ / molecule','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(hf,['output/deeply_coupled_timeseries_' num2str(prot1) '.pdf'],'Resolution',300)
