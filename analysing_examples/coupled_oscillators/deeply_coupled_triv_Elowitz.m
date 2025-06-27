%% deeply coupled oscillators
% 4 node circuit containing two repressilator

clc
clear
close all

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

%% create the model (examples/node4_deeply_coupled.mat) with an extra protease

clean_up_GRN
Ecoli = Cell('Elowitz');
Ecoli = Ecoli.add_node('N1','type1');
Ecoli = Ecoli.add_node('N2','type1');
Ecoli = Ecoli.add_node('N3','type1');
Ecoli = Ecoli.add_node('N4','type1');
Ecoli = Ecoli.add_regulator('Repression_in','N2','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Repression_in','N3','HILL','P_N2');
Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL','P_N3');
Ecoli = Ecoli.add_regulator('Repression_in','N4','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Repression_in','N3','HILL','P_N4');
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-10);
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-8);
set(getconfigset(Ecoli.data.Mobj),'Stoptime',1e4);
Ecoli.data.Accelerate = true;
Ecoli.set('P_N1','InitialAmount',1.000000e+03,'N1');
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N2'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N3'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N4'}];

%% simulation settings

% get the model
Mobj = Ecoli.get_model();

% start and end of the simulation
t_end = 1e4;
t_min = 0.1*t_end;
% number of data points
n = 1e3;

% simulation settings:
configsetobj = getconfigset(Mobj);
set(configsetobj,'SolverType','sundials','Stoptime',t_end)
set(get(configsetobj,'SolverOptions'),'AbsoluteTolerance',1e-10,...
    'RelativeTolerance',1e-8)

disp('accelerate')
tic
sbioaccelerate(Mobj);
toc

% prescribed promoter strength for N4 node
k1_N4 = linspace(0,50,100);

% store data: T1, err, T4, err, Amplitude in P_N1
data = zeros(length(k1_N4),5);

%% simulate: the effect of 'k1_N4'

tic
for i = 1:numel(k1_N4)
    
    % change the promoter strength
    set(sbioselect(Mobj,'Name','k1_N4'),'Value',k1_N4(i));
    
    % simulate
    simdata = sbiosimulate(Mobj);
    [t,c,names] = getdata(simdata);
    % get rid of the beginning
    start_pos = find(t>=t_min,1,'first');
    t(1:start_pos) = [];
    c(1:start_pos,:) = [];

    % calculate the time period of the independent oscillators
    [T1,err1] = calc_timeperiod(t,c(:,strcmp(names,'P_N1')));
    [T4,err4] = calc_timeperiod(t,c(:,strcmp(names,'P_N4')));

    % amplitude in P_N2
    A = max(c(:,strcmp(names,'P_N2')))-min(c(:,strcmp(names,'P_N2')));

    % save data
    data(i,:) = [T1,err1,T4,err4,A];
end
toc

%% plot the effect of 'k1_N4' on the time period

% load('output/repr_line','fitobject','nodes_repressilator','T');

!mkdir -p output

hf = figure;
hold on
plot(k1_N4,data(:,1),'LineWidth',2)
set(gca,'LineWidth',2,'Fontsize',16)
ylabel('\bf\boldmath$T$ / min','interpreter','latex','Fontsize',18)
xlabel('\bf\boldmath$k_{1,N4}$ / molecule/min','interpreter','latex','Fontsize',18)
scatter([0,15,30],interp1(k1_N4,data(:,1),[0,15,30]),'filled','o','SizeData',100)
exportgraphics(hf,'output/deeply_coupled_T_Elowitz.pdf','Resolution',300)

hf = figure;
plot(k1_N4,data(:,5),'LineWidth',2)
set(gca,'LineWidth',2,'Fontsize',16)
ylabel('\bf\boldmath$A\left(\left[P_{N2}\right]\right)$ / molecule','interpreter','latex','Fontsize',18)
xlabel('\bf\boldmath$k_{1,N4}$ / molecule/min','interpreter','latex','Fontsize',18)
exportgraphics(hf,'output/deeply_coupled_A_Elowitz.pdf','Resolution',300)

%% example simulation

k1_N4 = [30,15,0];

% start and end of the simulation
t_end = 2e3;
t_min = 1e3;

configsetobj = getconfigset(Mobj);
set(configsetobj,'Stoptime',t_end)

hf = figure;
ht = tiledlayout(3,1,'TileSpacing','tight');

for i = 1:numel(k1_N4)

    % change the promoter strength
    set(sbioselect(Mobj,'Name','k1_N4'),'Value',k1_N4(i));

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

    plot(t,c(:,strcmp(names,'P_N1'))/1e3,'LineWidth',2)
    plot(t,c(:,strcmp(names,'P_N2'))/1e3,'LineWidth',2)
    plot(t,c(:,strcmp(names,'P_N3'))/1e3,'LineWidth',2)
    plot(t,c(:,strcmp(names,'P_N4'))/1e3,'LineWidth',2)

    % title(['PROT2 = ' num2str(c_prot(i)) ' molecule'],'Fontsize',16)
    % no labels on the x axis
    if i ~= 3
        set(gca,'XTickLabel',[])
    end
    set(gca,'LineWidth',2,'Fontsize',16)

    if i == 1
        % to create a legend at the end so it will be on the top
        hax1 = gca;
    end

end
xlabel(ht,'\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel(ht,'\bf\boldmath$c$ / $10^3$ molecule','interpreter','latex','Fontsize',18)
legend(hax1,'[P_1]','[P_2]','[P_3]','[P_4]')
exportgraphics(hf,'output/deeply_coupled_Elowitz.pdf','Resolution',300)
