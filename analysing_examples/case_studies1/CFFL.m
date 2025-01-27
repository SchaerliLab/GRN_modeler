%% Coherent Feed Forward Loop
clc
close all
clear

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

clean_up_GRN

%% build the model

Ecoli = Cell('Elowitz');
solver = 'ode15s';
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-8)
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-6)
set(getconfigset(Ecoli.data.Mobj),'MaximumWallClock',10)
Ecoli = Ecoli.add_node('N1','type1');
Ecoli = Ecoli.add_node('N2','type1');
Ecoli = Ecoli.add_node('N3','type1');
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N3'}];
set(getconfigset(Ecoli.data.Mobj),'StopTime',1e3);
Ecoli = Ecoli.add_regulator('Activation_in','N2','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Activation_in','N3','HILL','P_N2');
Ecoli = Ecoli.add_regulator('Activation_in','N3','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Activation_out','N1','HILL','R1');

%% Parameter scan

Mobj = Ecoli.get_model();
set(getconfigset(Mobj),'SolverType',solver);

n_div = 100; % number of divisions
R1_values = linspace(0,30,n_div);
P_N3 = zeros(n_div,1); % output concentration

% simulation
for i = 1:n_div
    set(sbioselect(Mobj,'Name','R1'),'Value',R1_values(i));
    [t,c,names] = run_simulation(Mobj,solver);
    % save the data
    P_N3(i) = c(end);
end

% plot the data
figure
hold on
plot(R1_values,P_N3,'LineWidth',2)
xlabel('\bf\boldmath$\left[\textrm{R}_1\right]$ / \boldmath$\mu$M','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$\left[\textrm{P}_{N3}\right]$ / molecule','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)

!mkdir -p output
exportgraphics(gcf,'output/FFL.pdf','Resolution',300)
