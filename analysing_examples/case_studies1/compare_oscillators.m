%% compare the repressilator, reptolator and the actolator
clc
clear
close all

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

% we evaluate the simulations only after this time / min
start_time = 1e3;
stop_time = 1e4;

%% repressilator

% number of nodes
nodes_repressilator = 3:2:20;
% output data: duty cycle, time period / min
data_repressilator = zeros(numel(nodes_repressilator,2));

for i = 1:numel(nodes_repressilator)
    % build the model
    Ecoli = build_repressilator('Elowitz',nodes_repressilator(i));
    
    % set and run simulations
    [duty,dt_high,dt_low] = set_run(Ecoli,start_time,stop_time);

    % store the data
    data_repressilator(i,1) = duty*100;
    data_repressilator(i,2) = dt_high+dt_low;
end

%% reptolator

% number of nodes
nodes_reptolator = 4:4:20;
% output data: duty cycle, time period / min
data_reptolator = zeros(numel(nodes_reptolator,2));

for i = 1:numel(nodes_reptolator)
    % build the model
    Ecoli = build_reptolator('Elowitz',nodes_reptolator(i));
    
    % set and run simulations
    [duty,dt_high,dt_low] = set_run(Ecoli,start_time,stop_time);

    % store the data
    data_reptolator(i,1) = duty*100;
    data_reptolator(i,2) = dt_high+dt_low;
end

%% actolator

% number of nodes
nodes_actolator = 4:2:20;
% output data: duty cycle, time period / min
data_actolator = zeros(numel(nodes_actolator,2));

for i = 1:numel(nodes_actolator)
    % build the model
    Ecoli = build_actolator('Elowitz',nodes_actolator(i));
    
    % set and run simulations
    [duty,dt_high,dt_low] = set_run(Ecoli,start_time,stop_time);

    % store the data
    data_actolator(i,1) = duty*100;
    data_actolator(i,2) = dt_high+dt_low;
end

%% plot the data

!mkdir -p output

% duty cycle
figure
hold on
plot(nodes_repressilator,data_repressilator(:,1),'LineWidth',2,'LineStyle','--','Marker','o')
plot(nodes_reptolator,data_reptolator(:,1),'LineWidth',2,'LineStyle','--','Marker','o')
plot(nodes_actolator,data_actolator(:,1),'LineWidth',2,'LineStyle','--','Marker','o')
xlabel('\bf number of nodes','interpreter','latex','Fontsize',18)
ylabel('\bf duty cycle / \%','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(gcf,'output/duty_compare.pdf')

% time period
figure
hold on
plot(nodes_repressilator,data_repressilator(:,2),'LineWidth',2,'LineStyle','--','Marker','o')
plot(nodes_reptolator,data_reptolator(:,2),'LineWidth',2,'LineStyle','--','Marker','o')
plot(nodes_actolator,data_actolator(:,2),'LineWidth',2,'LineStyle','--','Marker','o')
xlabel('\bf number of nodes','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$T$ / min','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(gcf,'output/T_compare.pdf')

%% set and run simulations
function [duty,dt_high,dt_low] = set_run(Ecoli,start_time,stop_time)
% simulation settings
Ecoli.set('P_N1','InitialAmount',1.000000e+02,'N1');
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-8)
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-6)
set(getconfigset(Ecoli.data.Mobj),'Stoptime',stop_time);
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];

% run simulation
Mobj = Ecoli.get_model();
[t,c,~] = run_simulation(Mobj,'ode15s');

% get rid of the initial transient part
t_min = find(t>start_time,1,'first');
t(1:t_min) = [];
c(1:t_min,:) = [];

% calculate time period and duty cycle
% thres = get(sbioselect(Ecoli.data.regulator_models.Repression_in.Parameters,'Name','K'),'Value');
thres = (max(c)-min(c))/2;
[duty,dt_high,dt_low] = calc_duty(t,c,thres);
end