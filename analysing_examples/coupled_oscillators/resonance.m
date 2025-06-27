%% repressilator |- repressilator
% unidirectionally coupled

clc
clear
close all

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

%% build the circuit
clean_up_GRN
Ecoli = Cell('Elowitz');
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
Ecoli = Ecoli.add_regulator('Repression_in','N3','HILL','P_N4');
% Ecoli = Ecoli.add_protease('N1','PROT1','type1');
% Ecoli = Ecoli.add_protease('N2','PROT1','type1');
% Ecoli = Ecoli.add_protease('N3','PROT1','type1');
% Ecoli = Ecoli.add_protease('N4','PROT2','type1');
% Ecoli = Ecoli.add_protease('N5','PROT2','type1');
% Ecoli = Ecoli.add_protease('N6','PROT2','type1');
Ecoli.data.Accelerate = true;
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1.000000e-8);
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1.000000e-10);
set(getconfigset(Ecoli.data.Mobj),'Stoptime',1.000000e+04);
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1','P_N6'}];
Ecoli.set('P_N1','InitialAmount',1.000000e+03,'N1');
Ecoli.set('P_N4','InitialAmount',1.000000e+03,'N4');
% Ecoli.set('PROT1','InitialAmount',1.000000e+02,'PROT1');
% Ecoli.set('PROT2','InitialAmount',1.000000e+02,'PROT2');

% get the model
Mobj = Ecoli.get_model();

%% change the rule for the N3 -| N4 repression
% we need just a mild/partial repression here

% % the repression will be effective just in 'ratio' ratio:
% addparameter(Mobj,'ratio',1,'Units','dimensionless');
% 
% regulation = sbioselect(Mobj.Rule,'Name','Rule_HILL_N3|-N4');
% regulation.Rule = [regulation.Rule '*ratio + (1-ratio)'];

% higher K
set(sbioselect(Mobj.Parameters,'Name','K_N3|-N4'),'Value',1e2)

%% calculate time period and amplitude

target_spec = 'P_N1'; % name of the followed species in the bifurcation diagram
% start and end of the simulation
t_min = 1e4;
t_end = 1e4+t_min;
k1 = linspace(5,100,96*2); % values for the changed variable
c_name = 'k1_N6'; % the name of the changed variable 

% simulation settings:
configsetobj = getconfigset(Mobj);
set(configsetobj,'SolverType','sundials','Stoptime',t_end)

disp('accelerate')
sbioaccelerate(Mobj)

% T, A data in the columns
data = zeros(length(k1),2);

tic
disp('simulation')
for i = 1:numel(k1)

    %  change the given parameter concentration
    set(sbioselect(Mobj,'Name',c_name),'Value',k1(i))
    
    % simulate
    simdata = sbiosimulate(Mobj);
    [t,c,names] = getdata(simdata);
    % get rid of the beginning
    start_pos = find(t>=t_min,1,'first');
    c(1:start_pos,:) =  [];
    t(1:start_pos,:) =  [];

    % time period
    data(i,1) = calc_timeperiod(t,c(:,strcmp(names,'P_N6')));
    % amplitude
    data(i,2) = max(c(:,strcmp(names,target_spec)))-min(c(:,strcmp(names,target_spec)));

end
toc

%% plot
figure
plot(k1,data(:,2),'LineWidth',2)
xlabel('\bf\boldmath$k_{1,N6}$ / molecule/minute','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$A\left(\left[P_{N1}\right]\right)$ / molecule','interpreter','latex','Fontsize',18)
% Create textarrow
annotation(gcf,'textarrow',[0.419298245614035 0.373684210526316],...
    [0.632971291866029 0.796650717703349],'String',{'Resonance'},'Color','red','LineWidth',2,'FontSize',18);
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(gcf,'output/resonance.pdf','Resolution',300)

figure
plot(k1,data(:,1),'LineWidth',2)
xlabel('\bf\boldmath$k_{1,N6}$ / molecule/minute','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$T\left(\left[P_{N6}\right]\right)$ / minute','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(gcf,'output/resonance_T.pdf','Resolution',300)
