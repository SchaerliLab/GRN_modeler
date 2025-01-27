%% Repressilator
clc
close all
clear

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

% select the type of the model
modeltype = 'deterministic'; % 'stochastic' 'deterministic'

%% build the model
clean_up_GRN
Ecoli = Cell('Elowitz');
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
Ecoli.set('P_N1','InitialAmount',1.000000e+02,'N1'); % Elowitz

switch modeltype
    case 'deterministic'
        solver = 'deterministic';
    case 'stochastic'
        solver = 'adaptivesa';
end

%% Run simulation

Mobj = Ecoli.get_model();

switch modeltype
    case 'stochastic'
        configset = getconfigset(Mobj);
        Mobj = convert2irrev(Mobj);
        Mobj = correct_modifiers(Mobj);
        Ecoli.set_configset(Mobj,configset)
end

[t,c,names] = run_simulation(Mobj,solver);

%% plot

hf = figure;
hold on
for i = 1:numel(Ecoli.data.StatesToLog)
    plot(t,c(:,strcmp(names,Ecoli.data.StatesToLog{i})),'LineWidth',2)
end
xlabel('\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$c$ / molecule','interpreter','latex','Fontsize',18)
legend(Ecoli.data.StatesToLog,'Interpreter','none')
set(gca,'LineWidth',2,'Fontsize',16)
!mkdir -p output
exportgraphics(hf,['output/repressilator_' modeltype '.pdf'],'Resolution',300)


%% make graph
figure
Ecoli.data.node_fontsize = 16;
Ecoli.data.inner_regulator_size = 10;
Ecoli.make_graph();
axis off
exportgraphics(gcf,'output/repressilator_graph.pdf','Resolution',300)