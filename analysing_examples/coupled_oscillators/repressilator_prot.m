%% repressilator
% the effect of the protease concentration

clc
clear
close all

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

%% create the model

clean_up_GRN
Ecoli = Cell('Tomazou');
Ecoli = Ecoli.add_node('N1','type1');
Ecoli = Ecoli.add_node('N2','type1');
Ecoli = Ecoli.add_node('N3','type1');
Ecoli = Ecoli.add_regulator('Repression_in','N2','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Repression_in','N3','HILL','P_N2');
Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL','P_N3');
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N2'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N3'}];
set(getconfigset(Ecoli.data.Mobj),'Stoptime',1e3);
Ecoli = Ecoli.add_protease('N1','PROT1','type1');
Ecoli = Ecoli.add_protease('N2','PROT1','type1');
Ecoli = Ecoli.add_protease('N3','PROT1','type1');
Ecoli.set('P_N1','InitialAmount',1.000000e+03,'N1'); % Tomazou
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1.000000e-8);
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1.000000e-10);

%% calculate bifurcation diagram and Lyapunov exponents

for i = 1:2
    if i == 1
        Ecoli.set('P_N1','InitialAmount',1.000000e+03,'N1');
    else
        Ecoli.set('P_N1','InitialAmount',0,'N1');
    end

    % get the model
    Mobj = Ecoli.get_model();

    % settings
    p.n_workers = 1; % number of nodes in the parallelization
    p.target_spec = 'P_N1'; % name of the followed species in the bifurcation diagram
    % start and end of the simulation
    p.t_min = 1e5;
    p.t_end = 1e5+p.t_min;
    p.c = linspace(10,400,391); % concentrations for the changed variable
    p.c_name = 'PROT1'; % the name of the changed variable
    p.iszero = true;
    % calculate bif diagr
    [hf,exponents,pks] = calc_bifurcation_Lyapunov(Mobj,p);

    % xlim([p.c(1),p.c(end)])
    xlabel('\bf\boldmath$\left[\textrm{PROT}_1\right]$ / molecule','interpreter','latex','Fontsize',18)
    ylabel('\bf\boldmath$\max\left[(P_{N1})\right]$ / molecule','interpreter','latex','Fontsize',18)
    set(gca,'LineWidth',2,'Fontsize',16)
    if i == 1
        exportgraphics(hf,'output/repressilator_bifurcation_Lyap.pdf','Resolution',300)
        save('output/repressilator_bifurcation_Lyap.mat','exponents','pks')
    else
        exportgraphics(hf,'output/repressilator_bifurcation_zerostart_Lyap.pdf','Resolution',300)
        save('output/repressilator_bifurcation_zerostart_Lyap.mat','exponents','pks')
    end
end

% change back the initial condition
Ecoli.set('P_N1','InitialAmount',1.000000e+03,'N1');

%% Poincare map at a given point

% change the 'PROT1' concentration
for prot1 = [100,310]

% get the model
Mobj = Ecoli.get_model();

set(sbioselect(Mobj.Species,'Name','PROT1'),'Value',prot1)

% target concentartion for Poincare section
p.c_target = 1e2;
% name of the target species
p.target_name = 'uP_N1';
% transient time
p.t_min = 1e5;
% in inverse interpolation use n_interp in both directions
p.n_interp = 3;

% start and end of the simulation
t_end = p.t_min+1e5;

% simulation settings:
configsetobj = getconfigset(Mobj);
set(get(configsetobj,'RuntimeOptions'),'StatesToLog',{'P_N1','uP_N1','mRNA_N1'})
set(configsetobj,'SolverType','sundials','Stoptime',t_end)

sbioaccelerate(Mobj)

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
if prot1 == 100
    axis([0 2 2e4 4e4])
end
exportgraphics(hf,['output/repressilator_poincare_section_' num2str(prot1) '.pdf'],'Resolution',300)

% Poincare map
hf = figure;
scatter(c_cross_P_N1(1:end-1),c_cross_P_N1(2:end),'blue','filled')
xlabel('\bf\boldmath$\left[P_{N1}\right]$(n) / molecule','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$\left[P_{N1}\right]$(n+1) / molecule','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
if prot1 == 100
    axis([2e4 4e4 2e4 4e4])
end
exportgraphics(hf,['output/repressilator_poincare_map_' num2str(prot1) '.pdf'],'Resolution',300)

% phase space of a node
hf = figure;
plot3(c(pos_min:end,strcmp(names,'P_N1')),c(pos_min:end,strcmp(names,'mRNA_N1')),c(pos_min:end,strcmp(names,'uP_N1')))
xlabel('\bf\boldmath$\left[P_{N1}\right]$','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$\left[mRNA_{N1}\right]$','interpreter','latex','Fontsize',18)
zlabel('\bf\boldmath$\left[uP_{N1}\right]$','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
box on
exportgraphics(hf,['output/repressilator_phasespace_' num2str(prot1) '.pdf'],'Resolution',300)

% time series
hf = figure;
plot(t,c(:,strcmp(names,'P_N1')))
set(gca,'XLim',[0,1e4])
xlabel('\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$\left[P_{N1}\right]$ / molecule','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(hf,['output/repressilator_timeseries_' num2str(prot1) '.pdf'],'Resolution',300)

end

%% Example siulations

% 'PROT1' concentrations to show
prot1_values = [100,310,400];

hf = figure;
ht = tiledlayout(3,1,'TileSpacing','tight');

set(get(configsetobj,'RuntimeOptions'),'StatesToLog',{'P_N1','P_N2','P_N3'})

for i = 1:numel(prot1_values)

    prot1 = prot1_values(i);

    % start and end of the simulation
    if i == 3
        t_end = 1.1e3;
    else
        t_end = 2e3;
    end
    t_min = 1e3;

    % simulation time
    set(getconfigset(Mobj),'Stoptime',t_end)

    % change the 'PROT1' concentration
    set(sbioselect(Mobj.Species,'Name','PROT1'),'Value',prot1)

    % simulate
    simdata = sbiosimulate(Mobj);
    [t,c,names] = getdata(simdata);
    pos_min = find(t>t_min,1,'first');
    t(1:pos_min) = [];
    c(1:pos_min,:) = [];

    nexttile
    hold on
    plot(t,c(:,strcmp(names,'P_N1')),'LineWidth',2)
    plot(t,c(:,strcmp(names,'P_N2')),'LineWidth',2)
    plot(t,c(:,strcmp(names,'P_N3')),'LineWidth',2)
    set(gca,'LineWidth',2,'Fontsize',16)

    if i == 1
        % to create a legend at the end so it will be on the top
        hax1 = gca;
    end

    % % no labels on the x axis
    % if i == 1
    %     set(gca,'XTickLabel',[])
    % end

end
xlabel(ht,'\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel(ht,'\bf\boldmath$c$ / molecule','interpreter','latex','Fontsize',18)
legend(hax1,'[P_1]','[P_2]','[P_3]')
exportgraphics(hf,'output/repressilator_timeseries.pdf','Resolution',300)