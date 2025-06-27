%% coupled repressilators

clc
clear
close all

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

%% create the model (examples/coupled_repressilators2.mat)
clean_up_GRN
app.Ecoli = Cell();
app.Ecoli = app.Ecoli.add_node('N1','type1');
app.Ecoli = app.Ecoli.add_node('N2','type1');
app.Ecoli = app.Ecoli.add_node('N3','type1');
app.Ecoli = app.Ecoli.add_node('N4','type1');
app.Ecoli = app.Ecoli.add_node('N5','type1');
app.Ecoli = app.Ecoli.add_node('N6','type1');
app.Ecoli = app.Ecoli.add_regulator('Repression_in','N2','HILL','P_N1');
app.Ecoli = app.Ecoli.add_regulator('Repression_in','N3','HILL','P_N2');
app.Ecoli = app.Ecoli.add_regulator('Repression_in','N1','HILL','P_N3');
app.Ecoli = app.Ecoli.add_regulator('Repression_in','N5','HILL','P_N4');
app.Ecoli = app.Ecoli.add_regulator('Repression_in','N6','HILL','P_N5');
app.Ecoli = app.Ecoli.add_regulator('Repression_in','N4','HILL','P_N6');
app.Ecoli = app.Ecoli.add_protease('N1','PROT1','type1');
app.Ecoli = app.Ecoli.add_protease('N2','PROT1','type1');
app.Ecoli = app.Ecoli.add_protease('N3','PROT1','type1');
app.Ecoli = app.Ecoli.add_protease('N4','PROT2','type1');
app.Ecoli = app.Ecoli.add_protease('N5','PROT2','type1');
app.Ecoli = app.Ecoli.add_protease('N6','PROT2','type1');
app.Ecoli = app.Ecoli.add_protease('N2','PROT3','type1');
app.Ecoli = app.Ecoli.add_protease('N3','PROT3','type1');
app.Ecoli = app.Ecoli.add_protease('N5','PROT3','type1');
app.Ecoli = app.Ecoli.add_protease('N6','PROT3','type1');
app.Ecoli = app.Ecoli.add_protease('N1','PROT3','type1');
app.Ecoli = app.Ecoli.add_protease('N4','PROT3','type1');
app.Ecoli.data.Accelerate = true;
set(get(getconfigset(app.Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1.000000e-08);
set(get(getconfigset(app.Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1.000000e-10);
set(getconfigset(app.Ecoli.data.Mobj),'Stoptime',1.000000e+04);
app.Ecoli.data.StatesToLog = [app.Ecoli.data.StatesToLog, {'P_N1'}];
app.Ecoli.data.StatesToLog = [app.Ecoli.data.StatesToLog, {'P_N4'}];
set(sbioselect(app.Ecoli.data.Mobj,'Name','PROT3'),'Value',1.000000e+02);
set(sbioselect(app.Ecoli.data.Mobj,'Name','PROT1'),'Value',1.678322e+01);
set(sbioselect(app.Ecoli.data.Mobj,'Name','PROT2'),'Value',0.000000e+00);
app.Ecoli.set('P_N1','InitialAmount',1.000000e+03,'N1');
app.Ecoli.set('P_N4','InitialAmount',1.000000e+03,'N4');

%% simulation settings

% get the model
Mobj = app.Ecoli.get_model();

% start and end of the simulation
t_end = 2e4;
t_min = 1e4;

% % number of workers
% n_workers = 4;

% simulation settings:
configsetobj = getconfigset(Mobj);
set(configsetobj,'SolverType','sundials','Stoptime',t_end)
set(get(configsetobj,'SolverOptions'),'AbsoluteTolerance',1e-10,...
    'RelativeTolerance',1e-8)

set(get(configsetobj,'RuntimeOptions'),'StatesToLog',{'P_N1','P_N4'})

disp('accelerate')
tic
sbioaccelerate(Mobj);
toc

% prescribed 'PROT3' concentrations (common protease)
c_prot = 0:10:100;

% store data: T1, err, T4, err
data = zeros(length(c_prot),4);

% we set the individual protease concentrations to a given (different)
% value, to get different time periods
set(sbioselect(Mobj.Species,'Name','PROT1'),'Value',50)
set(sbioselect(Mobj.Species,'Name','PROT2'),'Value',200)

%% simulate

% % Get current parallel pool
% p = gcp('nocreate');
% % parallel computing
% if isempty(p) || p.NumWorkers ~= n_workers
%     delete(p);
%     % start a new parallel pool
%     parpool("Processes",n_workers);
% end

tic
for i = 1:length(c_prot)

    t_end_par = t_end;
    t_min_par = t_min;

    % get the model
    Mobj_par = Mobj;

    %  change the 'PROT1' concentration
    set(sbioselect(Mobj_par.Species,'Name','PROT3'),'Value',c_prot(i))


    % accelerate if it was not accelerated already
    sbioaccelerate(Mobj_par);

    % simulate
    simdata = sbiosimulate(Mobj_par);
    [t,c,names] = getdata(simdata);
    % get rid of the beginning
    start_pos = find(t>=t_min,1,'first');
    t(1:start_pos) = [];
    c(1:start_pos,:) = [];

    % calculate the time period of the coupled oscillators
    [T1,err1] = calc_timeperiod(t,c(:,strcmp(names,'P_N1')));
    [T4,err4] = calc_timeperiod(t,c(:,strcmp(names,'P_N4')));

    % save data
    data(i,:) = [T1,err1,T4,err4];

end
toc

%% plot the time period, synchronization

!mkdir -p output

hf = figure;
hold on
errorbar(c_prot,data(:,1),data(:,2),'LineWidth',2)
errorbar(c_prot,data(:,3),data(:,4),'LineWidth',2)
set(gca,'LineWidth',2,'Fontsize',16)
xlabel('\bf\boldmath$\left[\textrm{PROT}_3\right]$ / molecule','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$T$ / min','interpreter','latex','Fontsize',18)
legend('N_1','N_4')
exportgraphics(hf,'output/synchronization.pdf','Resolution',300)

%% example simulation

c_prot = [0,40,100];

% start and end of the simulation
t_end = 2e4;
t_min = 1e4;

configsetobj = getconfigset(Mobj);
set(configsetobj,'Stoptime',t_end)

hf = figure;
ht = tiledlayout(3,1,'TileSpacing','tight');

for i = 1:numel(c_prot)

    % change the 'PROT2' concentration
    set(sbioselect(Mobj.Species,'Name','PROT3'),'Value',c_prot(i))

    % simulate
    simdata = sbiosimulate(Mobj);
    [t,c,names] = getdata(simdata);
    % get rid of the beginning
    start_pos = find(t>=t_min,1,'first');
    t(1:start_pos) = [];
    c(1:start_pos,:) = [];

    % plot the stuff
    nexttile
    hold on
    plot(t,c(:,strcmp(names,'P_N1'))/1e5,'LineWidth',2)

    % title(['PROT3 = ' num2str(c_prot(i)) ' molecule'],'Fontsize',16)
    % no labels on the x axis
    if i ~= 3
        set(gca,'XTickLabel',[])
    end
    set(gca,'LineWidth',2,'Fontsize',16)

    if i == 1
        % to create a legend at the end so it will be on the top
        hax1 = gca;
    end

    % we show just the top, because this is the interesting part
    set(gca,'Ylim',[.6,1])

end
xlabel(ht,'\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel(ht,'\bf\boldmath$\left[P_1\right]$ / $10^5$ molecule','interpreter','latex','Fontsize',18)
% legend(hax1,'[P_1]','[P_2]','[P_3]','[P_4]')
exportgraphics(hf,'output/synchronise_time.pdf','Resolution',300)

%% phase synchronization

c_prot = [0,1,2,4,8,16]; % molecule

% change the initial conditions
set(sbioselect(Mobj,'Name','PROT1'),'Value',50) % same PROT1, PROT2
set(sbioselect(Mobj,'Name','PROT2'),'Value',50)
set(sbioselect(Mobj,'Name','PROT3'),'Value',0)
set(sbioselect(Mobj,'Name','P_N1'),'Value',1000) % slightly different N1,N4
set(sbioselect(Mobj,'Name','P_N4'),'Value',600)

% start and end of the simulation
t_end = 3e4;
t_min = 0;

% simulation settings:
configsetobj = getconfigset(Mobj);
set(configsetobj,'SolverType','sundials','Stoptime',t_end)
set(get(configsetobj,'RuntimeOptions'),'StatesToLog',{'P_N1','P_N4'})

disp('accelerate')
tic
sbioaccelerate(Mobj);
toc

hf = figure;
tiledlayout(2,1,"TileSpacing","tight")
nexttile
hold on
for i = 1:numel(c_prot)

    % PROT3
    set(sbioselect(Mobj,'Name','PROT3'),'Value',c_prot(i))

    % simulate
    [t,c,names] = sbiosimulate(Mobj);

    % find the peaks
    [~,locs1] = findpeaks(c(:,strcmp(names,'P_N1')),t);
    [~,locs2] = findpeaks(c(:,strcmp(names,'P_N4')),t);
    % maximum number of peaks: 
    nmax = min(numel(locs1),numel(locs2));
    % difference in the peak locations
    dlocs = abs(locs1(1:nmax)-locs2(1:nmax));
    % we do not plot the first peak
    plot(locs1(1:nmax),[nan;dlocs(2:end)],'LineWidth',2)
end
xlabel('\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$\Delta t$ / min','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
hl = legend(num2str(c_prot(:)),'Fontsize',8);
title(hl,'\bf\boldmath$\left[\textrm{PROT}_3\right]$','interpreter','latex')

% one extra time series for demonstration
t_end = 1e3;
nexttile
hold on
% % PROT3
% set(sbioselect(Mobj,'Name','PROT3'),'Value',4)
% % simulate
% [t,c,names] = sbiosimulate(Mobj);
t_endpos = find(t<t_end,1,'last');
t(t_endpos:end) = [];
c(t_endpos:end,:) = [];
plot(t,c(:,strcmp(names,'P_N1')),'LineWidth',2)
plot(t,c(:,strcmp(names,'P_N4')),'LineWidth',2)
xlabel('\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$\left[P_1\right]$ / mol.','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)

exportgraphics(hf,'output/synchronise_phase.pdf','Resolution',300)