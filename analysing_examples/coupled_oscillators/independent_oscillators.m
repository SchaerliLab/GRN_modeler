%% independent oscillators
% interaction between two independent repressilator and a not gate
% beat phenomenon

clc
clear
close all

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

%% create the model (examples/coupled_repressilators_not.mat)

clean_up_GRN
Ecoli = Cell('Tomazou');
Ecoli = Ecoli.add_node('N1','type1');
Ecoli = Ecoli.add_node('N2','type1');
Ecoli = Ecoli.add_node('N3','type1');
Ecoli = Ecoli.add_node('N4','type1');
Ecoli = Ecoli.add_node('N5','type1');
Ecoli = Ecoli.add_node('N6','type1');
Ecoli = Ecoli.add_regulator('Repression_in','N2','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Repression_in','N3','HILL','P_N2');
Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL','P_N3');
Ecoli = Ecoli.add_regulator('Repression_in','N5','HILL','P_N4');
Ecoli = Ecoli.add_regulator('Repression_in','N6','HILL','P_N5');
Ecoli = Ecoli.add_regulator('Repression_in','N4','HILL','P_N6');
Ecoli = Ecoli.add_protease('N1','PROT1','type1');
Ecoli = Ecoli.add_protease('N2','PROT1','type1');
Ecoli = Ecoli.add_protease('N3','PROT1','type1');
Ecoli = Ecoli.add_protease('N4','PROT2','type1');
Ecoli = Ecoli.add_protease('N5','PROT2','type1');
Ecoli = Ecoli.add_protease('N6','PROT2','type1');
Ecoli.data.Accelerate = true;
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1.000000e-08);
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1.000000e-10);
Ecoli = Ecoli.add_node('N7','type1');
Ecoli = Ecoli.add_regulator('Repression_in','N7','HILL','P_N1');
set(getconfigset(Ecoli.data.Mobj),'Stoptime',1.000000e+04);
Ecoli.set('P_N1','InitialAmount',1.000000e+03,'PROT1');
Ecoli.set('P_N4','InitialAmount',1.000000e+03,'PROT2');
Ecoli = Ecoli.add_regulator('Repression_in','N7','HILL','P_N4');
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N7'}];
Ecoli.set('PROT1','InitialAmount',1.640000e+02,'PROT1');

%% simulation settings

% get the model
Mobj = Ecoli.get_model();

% start and end of the simulation
t_end = 1e4;
t_min = 0.1*t_end;
% number of data points
n = 1e3;

% minimum number of peaks
n_minpeak = 3;

% number of workers
n_workers = 6;

% simulation settings:
configsetobj = getconfigset(Mobj);
set(configsetobj,'SolverType','sundials','Stoptime',t_end)
set(get(configsetobj,'SolverOptions'),'AbsoluteTolerance',1e-10,...
    'RelativeTolerance',1e-8)
% % output time points 
% configsetobj.SolverOptions.OutputTimes = linspace(t_min,t_end,n);

set(get(configsetobj,'RuntimeOptions'),'StatesToLog',{'P_N1','P_N4','P_N7'})

disp('accelerate')
tic
sbioaccelerate(Mobj);
toc

% prescribed 'PROT1' concentrations
c_prot = 10:5:200;
% c_prot should not be 'PROT2'
c_prot(c_prot==100) = [];

% store data: T1, err, T4, err, T7, err
data = zeros(length(c_prot),6);

%% simulate

% Get current parallel pool
p = gcp('nocreate');
% parallel computing
if isempty(p) || p.NumWorkers ~= n_workers
    delete(p);
    % start a new parallel pool
    parpool("Processes",n_workers);
end

tic
parfor i = 1:length(c_prot)

    t_end_par = t_end;
    t_min_par = t_min;
    n_par = n;

    % get the model
    Mobj_par = Mobj;

    %  change the 'PROT1' concentration
    set(sbioselect(Mobj_par.Species,'Name','PROT1'),'Value',c_prot(i))

    % number of peaks
    n_peaks = 0;

    % we will increase the length of the simulation if we have too few
    % peaks
    while n_peaks < n_minpeak || n_peaks == 0

        sbioaccelerate(Mobj_par);

        % [t,c,names] = sbiosimulate(Mobj);
        simdata = sbiosimulate(Mobj_par);
        [t,c,names] = getdata(simdata);
        % get rid of the beginning
        start_pos = find(t>=t_min,1,'first');
        t(1:start_pos) = [];
        c(1:start_pos,:) = [];

        % calculate the time period of the independent oscillators
        [T1,err1] = calc_timeperiod(t,c(:,strcmp(names,'P_N1')));
        [T4,err4] = calc_timeperiod(t,c(:,strcmp(names,'P_N4')));

        % resample th simulation
        newSimData = resample(simdata,linspace(t_min_par,t_end_par,n_par));
        dt = (t_end_par-t_min_par)/(n_par-1);
        [t_resampled,c_resampled,~] = getdata(simdata);
        % smooth over the time period of the fast oscillators
        a = smooth(c_resampled(:,strcmp(names,'P_N7')),100*ceil(max(T1,T4)/dt));
        % plot(t_resampled,a)
        [T7,err7,n_peaks] = calc_timeperiod(t_resampled,a,1/(abs(1/T1-1/T4))/4);

        % disp(1/abs(1/T1-1/T4))
        % disp(T7)
        % disp(n_peaks)

        % save data
        data(i,:) = [T1,err1,T4,err4,T7,err7];

        % set up a 10 times longer simulation, we will run it again if we
        % do not have enough peaks
        ratio = ceil(n_minpeak/n_peaks);
        if n_peaks == 0
            ratio = 10;
        end
        t_end_par = ratio*t_end_par;
        t_min_par = ratio*t_min_par;
        n_par = ratio*n_par;

        set(getconfigset(Mobj_par),'Stoptime',t_end_par)
    end
end
toc

%% plot the beats

!mkdir -p output

% base oscillator time period
T4 = mean(data(:,3));
% half position
pos = find(data(:,1)>T4,1,'last');

hf = figure;
fplot(@(x) 1./abs(1./T4-1./x),[min(data(:,1)),max(data(:,1))],'LineWidth',2)
hold on
errorbar(data(1:pos,1),data(1:pos,5),data(1:pos,6),'color',"#D95319",'LineWidth',2)
errorbar(data(pos+1:end,1),data(pos+1:end,5),data(pos+1:end,6),'color',"#D95319",'LineWidth',2)
set(gca,'LineWidth',2,'Fontsize',16)
xlabel('\bf\boldmath$T$ / min','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$T_{beat}$ / min','interpreter','latex','Fontsize',18)
legend('Theoretical','Simulated')
exportgraphics(hf,'output/beat.pdf','Resolution',300)

%% example simulation

% change the 'PROT1' concentration
set(sbioselect(Mobj.Species,'Name','PROT1'),'Value',164)

% simulate
simdata = sbiosimulate(Mobj);
[t,c,names] = getdata(simdata);

%% plot

% example time series
tmin1 = 4500;
tmax1 = 5000;

tmin2 = 6750;
tmax2 = 7250;

hf = figure;
tlayout = tiledlayout(1,2,'TileSpacing','tight');
bgAx = axes(tlayout,'XTick',[],'YTick',[],'Box','off');
bgAx.Layout.TileSpan = [1 2];
ax1 = axes(tlayout);
hold on
plot(ax1,t,c(:,strcmp(names,'P_N4')),'LineWidth',2)
plot(ax1,t,c(:,strcmp(names,'P_N1')),'LineWidth',2)
xline(ax1,tmax1,'--','LineWidth',2);
ax1.Box = 'off';
xlim(ax1,[tmin1 tmax1])
set(ax1,'LineWidth',2,'Fontsize',16)

% Create second plot
ax2 = axes(tlayout);
ax2.Layout.Tile = 2;
hold on
plot(ax2,t,c(:,strcmp(names,'P_N4')),'LineWidth',2)
plot(ax2,t,c(:,strcmp(names,'P_N1')),'LineWidth',2)
xline(ax2,tmin2,'--','LineWidth',2);
ax2.YAxis.Visible = 'off';
ax2.Box = 'off';
xlim(ax2,[tmin2 tmax2])
legend('[P_4]','[P_1]')
set(ax2,'LineWidth',2,'Fontsize',16)

% Link the axes
linkaxes([ax1 ax2], 'y')

tlayout_children = get(tlayout,'children');
set(tlayout_children(end),'LineWidth',2,'Fontsize',16)

xlabel(tlayout,'\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel(tlayout,'\bf\boldmath$c$ / molecule','interpreter','latex','Fontsize',18)
exportgraphics(hf,'output/beat_oscillators.pdf','Resolution',300)

% hf = figure;
% hold on
% plot(t,c(:,strcmp(names,'P_N4')),'LineWidth',2)
% plot(t,c(:,strcmp(names,'P_N1')),'LineWidth',2)
% xlabel('\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
% ylabel('\bf\boldmath$c$ / molecule','interpreter','latex','Fontsize',18)
% legend('[P_4]','[P_1]')
% set(gca,'LineWidth',2,'Fontsize',16)
% exportgraphics(hf,'output/beat_oscillators.pdf','Resolution',300)

% beat
hf = figure;
hold on
plot(t,c(:,strcmp(names,'P_N7')),'LineWidth',2,'Color','#EDB120')
xlabel('\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$c$ / molecule','interpreter','latex','Fontsize',18)
legend('[P_7]')
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(hf,'output/beat_output.pdf','Resolution',300)