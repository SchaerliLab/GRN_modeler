%% compare solvers on the example of the repressilator with different node numbers
clc
clear
close all

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

%% repressilator

% number of nodes
nodes_repressilator = 3:2:20;

% solvers types what we will examine
solvers = {'ode15s', 'ode45', 'sundials', 'ode15s_acc', 'ode45_acc', 'sundials_acc','lsoda','radau5','adaptivesa'};%,'tauleap'};
% solvers = {'lsoda','radau5','adaptivesa','stochastic','directMethod','tauleap','hybrid','hybridlsoda','hybridode45','sde'};


% simulation time / second
simtime = zeros(length(nodes_repressilator),size(solvers,2));

% simulation length
stop_time = 1e5;

for i = 1:numel(nodes_repressilator)
    % build the model
    Ecoli = build_repressilator('Elowitz',nodes_repressilator(i));

    % simulation settings
    Ecoli.set('P_N1','InitialAmount',1.000000e+02,'N1');
    set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-8)
    set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-6)
    set(getconfigset(Ecoli.data.Mobj),'Stoptime',stop_time);
    Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];
    % output times
    set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'OutputTimes',[0,stop_time]);

    for j = 1:size(solvers,2)

        % get the model
        Mobj = Ecoli.get_model();

        % accelerate the "acc" versions
        if contains(solvers{j},'acc')
            sbioaccelerate(Mobj)
        end

        % for COPASI solvers
        if ~any(contains(solvers{j},'ode') | contains(solvers{j},'sundials'))
            configset = getconfigset(Mobj);
            Mobj = convert2irrev(Mobj);
            Mobj = correct_modifiers(Mobj);
            Ecoli.set_configset(Mobj,configset)
        end

        % measure the simulation time
        simtime(i,j) = timeit(@() run_simulation(Mobj,solvers{j}));
        % tic
        % run_simulation(Mobj,solvers{j});
        % simtime(i,j) = toc;

        % against memory leakage with python
        if ~any(contains(solvers{j},'ode') | contains(solvers{j},'sundials'))
            terminate(pyenv)
            pyenv('ExecutionMode', 'OutOfProcess');
        end

    end
end


%% plot the data

!mkdir -p output


% cmap = lines(size(solvers,2));  % Generates a colormap with distinct colors

cmap = [...
    [0.1216, 0.4667, 0.7059]; % Blue
    [1.0000, 0.6824, 0.1647]; % Orange
    [0.1725, 0.6275, 0.1725]; % Green
    [0.8392, 0.1529, 0.1569]; % Red
    [0.5804, 0.4039, 0.7412]; % Purple
    [0.5490, 0.3373, 0.2941]; % Brown
    [0.8902, 0.4667, 0.7608]; % Pink
    [0.4980, 0.4980, 0.4980]; % Gray
    [0.9333, 0.8667, 0.1412]  % Yellow
];

% duty cycle
figure
hold on
for i = 1:size(solvers,2)
    plot(nodes_repressilator,simtime(:,i),'LineWidth',2,'Color',cmap(i,:))
end
xlabel('\bf number of nodes','interpreter','latex','Fontsize',18)
ylabel('\bf simulation time / seconds','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
% set(gca,'Xlim',[3,20])
legend(solvers,'Interpreter','none','Location','best','FontSize',12)
exportgraphics(gcf,'output/solver_time.pdf')

save('output/simtime','simtime')