%% Creating computers using coupled CRISPRlators
clc
clear
close all

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))


%% build the model

clean_up_GRN
Ecoli = Cell('CRISPR');  
Ecoli = Ecoli.add_node('N1','type1');
Ecoli = Ecoli.add_node('N2','type1');
Ecoli = Ecoli.add_node('N3','type1');
Ecoli = Ecoli.add_regulator('Repression_in','N2','NOHILL','sgRNA_N1');
Ecoli = Ecoli.add_regulator('Repression_in','N3','NOHILL','sgRNA_N2');
Ecoli = Ecoli.add_regulator('Repression_in','N1','NOHILL','sgRNA_N3');
Ecoli = Ecoli.add_node('N4','type1');
Ecoli = Ecoli.add_node('N5','type1');
Ecoli = Ecoli.add_node('N6','type1');
Ecoli = Ecoli.add_regulator('Repression_in','N5','NOHILL','sgRNA_N4');
Ecoli = Ecoli.add_regulator('Repression_in','N6','NOHILL','sgRNA_N5');
Ecoli = Ecoli.add_regulator('Repression_in','N4','NOHILL','sgRNA_N6');

Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N2'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N3'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N4'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N5'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N6'}];

Ecoli.set('dCas','InitialAmount',2.0000e+02,'N6'); % change dCas here
% Ecoli.set('a1_N1','Value',1.1,'N1');
Ecoli.set('a1_N2','Value',1.1,'N2');
Ecoli.set('a1_N6','Value',1.2,'N6');

set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-10);
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-8);
set(getconfigset(Ecoli.data.Mobj),'Stoptime',1e5);
Ecoli.data.Accelerate = true;

set(get(getconfigset(Ecoli.data.Mobj),'CompileOptions'), 'UnitConversion', false);
set(get(getconfigset(Ecoli.data.Mobj),'CompileOptions'), 'DimensionalAnalysis', false);


solver = 'sundials';%'adaptivesa';

% get the simbiology model
Mobj = Ecoli.get_model();

addparameter(Mobj,'K_activity',1e4);

%% coupling through burden

% scale the activity
addparameter(Mobj,'activity','Constant',false);
addrule(Mobj,'activity = 1/(1+(([P_N1]+[P_N2]+[P_N3]+[P_N4]+[P_N5]+[P_N6])/K_activity)^2)','repeatedAssignment');

% % production rates a1, k_P
% set(sbioselect(Mobj.Parameters,'Name','k_P'),'Constant',false)
% val = get(sbioselect(Mobj.Parameters,'Name','k_P'),'Value');
% addrule(Mobj,['k_P = ' num2str(val) '*activity'],'repeatedAssignment');
for i = 1:6
    name = ['a1_N',int2str(i)];
    set(sbioselect(Mobj.Parameters,'Name',name),'Constant',false)
    val = get(sbioselect(Mobj.Parameters,'Name',name),'Value');
    addrule(Mobj,[name ' = ' num2str(val) '*activity'],'repeatedAssignment');
end

% set(get(getconfigset(Mobj),'RuntimeOptions'),'StatesToLog',{'P_N1','P_N2','P_N3','P_N4','P_N5','P_N6','activity'})

% % if we just decrease the rate constants
% for i = 1:6
%     name = ['a1_N',int2str(i)];
%     val = get(sbioselect(Mobj.Parameters,'Name',name),'Value');
%     set(sbioselect(Mobj.Parameters,'Name',name),'Value',val/4)
% end

%% The effeft of dCas concentration and the burden on coupling

% transient part of the simulation
transient = 0.15;

% output time
configsetobj = getconfigset(Mobj);
% set(configsetobj,'SolverType',solver);
configsetobj.SolverOptions.OutputTimes = linspace(transient*1e5,1e5,1e5);

% if strcmp(solver,'adaptivesa')
%     % convert to irreversible
%     configset = getconfigset(Mobj);
%     Mobj = convert2irrev(Mobj);
%     Mobj = correct_modifiers(Mobj);
%     Ecoli.set_configset(Mobj,configset)
% else
% 
%     sbioaccelerate(Mobj);
% end

sbioaccelerate(Mobj);

% changed parameters with their range
dcas = 200:100:900;
K_activity = logspace(3,5,3);

[DCAS,K] = meshgrid(dcas,K_activity);
% maximum abs corr coefficient
C = zeros(size(DCAS));
% phase locking value
PLV = zeros(size(DCAS));

tic
for i = 1:size(DCAS,1)
    % change the parameters
    set(sbioselect(Mobj.Parameters,'Name','K_activity'),'Value',K(i,1));
    for j = 1:size(DCAS,2)
        % change dCas
        set(sbioselect(Mobj.Species,'Name','dCas'),'Value',DCAS(i,j));

        % simulation
        [t,c,names] = run_simulation(Mobj,solver);

        % calculate the correlation matrix between the time series
        M = corr(c);
        C(i,j) = max(abs(M(1,4:end)));

        % time period
        T = calc_timeperiod_avr_cross(t,c(:,1));

        if isnan(T)
            PLV(i,j) = nan;
        else

            % Get the main frequency components
            fs = length(t)/t(end);  % Correct sampling rate calculation
            f0 = 1/T;               % Fundamental frequency
            bw = 0.5*f0;            % Bandwidth

            % Normalized frequency range for Butterworth filter
            low_cutoff = max(0, (f0 - bw/2)/(fs/2));  % Ensure cutoff > 0
            high_cutoff = min(1, (f0 + bw/2)/(fs/2)); % Ensure cutoff < 1

            % Design 3rd-order Butterworth bandpass filter
            [b, a] = butter(3, [low_cutoff, high_cutoff], 'bandpass');

            % filtering
            c1 = filtfilt(b,a,c(:,1));  % zero-phase filtering
            c4 = filtfilt(b,a,c(:,4));  % zero-phase filtering

            % phase calculation
            phase1 = angle(hilbert(c1));
            phase4 = angle(hilbert(c4));

            %
            PLV(i,j) = abs(mean(exp(1i * (phase1- phase4))));
        end

    end
end
toc

% when it is not oscillating, it should be nan
C(isnan(PLV)) = nan;

% change back dCAs
set(sbioselect(Mobj.Species,'Name','dCas'),'Value',2e2);

%% plot

!mkdir -p output

figure
hold on
for i = 1:numel(K_activity)
    plot(dcas,C(i,:),'-o','LineWidth',2)
end
set(gca,'LineWidth',1.5,'Fontsize',16)
xlabel('\bf [dCas] / molecule','Interpreter','latex','FontSize',18)
ylabel('\bf correlation strength','Interpreter','latex','FontSize',18)
hl = legend(num2str(log10(K_activity).'),'Location','southwest');
title(hl,'\boldmath$\log_{10}(K)$','Interpreter','latex')
ytickformat('%.1f')
set(gca,'YLim',[0,1])
exportgraphics(gcf,'output/coupling_burden.pdf','resolution',300)

figure
hold on
for i = 1:numel(K_activity)
    plot(dcas,PLV(i,:),'-o','LineWidth',2)
end
set(gca,'LineWidth',1.5,'Fontsize',16)
xlabel('\bf [dCas] / molecule','Interpreter','latex','FontSize',18)
ylabel('\bf Phase Locking Value','Interpreter','latex','FontSize',18)
hl = legend(num2str(log10(K_activity).'),'Location','southwest');
title(hl,'\boldmath$\log_{10}(K)$','Interpreter','latex')
ytickformat('%.1f')
set(gca,'YLim',[0,1])
exportgraphics(gcf,'output/coupling_burden_PLV.pdf','resolution',300)

%% Deterministic simulation for comparision

% configsetobj = getconfigset(Mobj);
% 
% configsetobj.SolverOptions.OutputTimes = linspace(0,1e4,1e4);
% 
% Mobj.Parameters(17).Value = 1e4;
% 
% [t,c,names] = run_simulation(Mobj,solver);

%% create C++ files for stochastic simulations

% copy the model when we modify the names
Mobj2 = copyobj(Mobj);

% we will write the rules explicitely
delete(Mobj2.Rules);    

% Get the stoichiometric matrix
% when we convert reversible reactions to irreversible
Mobj2 = convert2irrev(Mobj2);
Mobj2 = correct_modifiers(Mobj2);

% === parameters ===
% create parameters.h header file
write_parameters_header(Mobj2)
% create parameters.cpp
write_parameters_cpp(Mobj2)

% === reactions.h ===
% scale with activity
text = {'// 1.0/(1.0+pow((c[P_N1]+c[P_N2]+c[P_N3]+c[P_N4]+c[P_N5]+c[P_N6])/K_activity,2.0));';...
    'double activity = 1.0/(1.0+pow((c[5]+c[10]+c[15]+c[23]+c[28]+c[33])/p.K_activity,2.0));'};
% create the reaction header file
write_reaction_header(Mobj2,text)

% change 'a1' with activity
% read the whole file into a string
txt = fileread('reactions.h');
% replace every occurrence of "p.a1" with "activity*p.a1"
txt = regexprep(txt, '\<p\.a1_', 'activity*p.a1_');
% write back to the same file
fid = fopen('reactions.h','w');
fwrite(fid, txt);
fclose(fid);


%% create the mex file

mex CXXFLAGS="$CXXFLAGS -std=c++20 -Wall -pedantic -O3 -fopenmp" ...
    LDFLAGS="$LDFLAGS -fopenmp" ...
    -output reaction_sim reaction_sim.cpp parameters.cpp

%% parameter settings

% parameters into the p structure
for i = 1:numel(Mobj.Parameters)
    p.(Mobj.Parameters(i).Name) = Mobj.Parameters(i).Value;
end
% initial concentrations
c0 =( zeros(numel(Mobj.Species),1));
for i = 1:numel(Mobj.Species)
    c0(i) = Mobj.Species(i).Value;
end

% simulation settings
p.sigma_inf = 1e-2;
p.tau = 10;
p.t_sim = 1e5; % end of simulation
n_data = p.t_sim ;% output data number
p.n_step = n_data*1e2; % number of steps during the solution
p.n_data = p.n_step/round(p.n_step/n_data); % output data number

% ratio of the transient part
transient = 0.15;
start_pos = round(transient*n_data+1);

p.K_activity = 1e4;

% position of the proteins in the concentration vector
% we are following only these species
protein_ID = find(contains({Mobj.Species.Name},'P_N')).'-1; % -1: C matrix indexing

%%
% 
% clear reaction_sim % clear static variables in functions
% 
% % run stochastic simulation
% tic
% [c_stoch,t_stoch] = reaction_sim(p,c0,protein_ID);
% toc
% 
% figure
% hold on
% plot(t,c(:,1))
% plot(t_stoch,c_stoch(:,1))
% % plot(t_stoch,c_stoch(:,4))

%% Effect of the coupling through burden

% changing activity half saturation and noise stregth
K_activity = logspace(3,5,3);
sigma_inf = logspace(-2,-4,6);

[S,K] = meshgrid(sigma_inf,K_activity);
C = zeros(size(K));
PLV = zeros(size(K));

tic
for i = 1:size(K,1)
    for j = 1:size(K,2)
        p.K_activity = K(i,j);
        p.sigma_inf = S(i,j);
        
        % simulation
        clear reaction_sim
        [c_stoch,t_stoch] = reaction_sim(p,c0,protein_ID([1,4:6]));

        % get rid of transient
        t_stoch(1:start_pos-1) = [];
        c_stoch(1:start_pos-1,:) = [];

        % correlation matrix
        M = corr(c_stoch);

        C(i,j) = max(abs(M(1,2:end)));

        T = calc_timeperiod_avr_cross(t_stoch,c_stoch(:,1));

        % Get the main frequency components
        fs = length(t_stoch)/t_stoch(end);  % Correct sampling rate calculation
        f0 = 1/T;               % Fundamental frequency
        bw = 0.5*f0;            % Bandwidth

        % Normalized frequency range for Butterworth filter
        low_cutoff = max(0, (f0 - bw/2)/(fs/2));  % Ensure cutoff > 0
        high_cutoff = min(1, (f0 + bw/2)/(fs/2)); % Ensure cutoff < 1

        % Design 3rd-order Butterworth bandpass filter
        [b, a] = butter(3, [low_cutoff, high_cutoff], 'bandpass');

        % filtering
        c1 = filtfilt(b,a,c_stoch(:,1));  % zero-phase filtering
        c4 = filtfilt(b,a,c_stoch(:,2));  % zero-phase filtering

        % phase calculation2
        phase1 = angle(hilbert(c1));
        phase4 = angle(hilbert(c4));

        % phase locking value
        PLV(i,j) = abs(mean(exp(1i * (phase1- phase4))));

    end
end
toc

%% plot

!mkdir -p output

% figure
% hold on
% for i = 1:numel(K_activity)
%     plot(sigma_inf,C(i,:),'-o','LineWidth',2)
% end
% set(gca,'LineWidth',1.5,'Fontsize',16)
% xlabel('\boldmath$\sigma_{\infty}$','Interpreter','latex','FontSize',18)
% ylabel('\bf correlation strength','Interpreter','latex','FontSize',18)
% hl = legend(num2str(log10(K_activity).'),'Location','southwest');
% title(hl,'\boldmath$\log_{10}(K)$','Interpreter','latex')
% ytickformat('%.1f')
% set(gca,'YLim',[0,1],'XScale','log')
% exportgraphics(gcf,'output/coupling_noise.pdf','resolution',300)

figure
hold on
for i = 1:numel(K_activity)
    plot(sigma_inf,PLV(i,:),'-o','LineWidth',2)
end
set(gca,'LineWidth',1.5,'Fontsize',16)
xlabel('\boldmath$\sigma_{\infty}$','Interpreter','latex','FontSize',18)
ylabel('\bf Phase Locking Value','Interpreter','latex','FontSize',18)
hl = legend(num2str(log10(K_activity).'),'Location','southwest');
title(hl,'\boldmath$\log_{10}(K)$','Interpreter','latex')
ytickformat('%.1f')
set(gca,'YLim',[0,1],'XScale','log')
exportgraphics(gcf,'output/coupling_noise_PLV.pdf','resolution',300)

%% Example stochastic simulations


clear reaction_sim % clear static variables in functions

p.t_sim = 1e5;%t(end); % end of simulation
n_data = p.t_sim ;% output data number
p.n_step = n_data*1e2;%*10; % number of steps during the solution
p.n_data = p.n_step/round(p.n_step/n_data); % output data number

p.K_activity = 1e4;

% ratio of the transient part
transient = 0.15;
start_pos = round(transient*n_data+1);

for sigma_inf = [1e-2,1e-3,1e-4]

    % change noise
    p.sigma_inf = sigma_inf;

    % run stochastic simulation
    tic
    [c_stoch,t_stoch] = reaction_sim(p,c0,protein_ID([1,4:6]));
    toc

    % % get rid of transient
    % t_stoch(1:start_pos-1) = [];
    % c_stoch(1:start_pos-1,:) = [];

    T = calc_timeperiod_avr_cross(t_stoch,c_stoch(:,1));

    % Get the main frequency components
    fs = length(t_stoch)/t_stoch(end);  % Correct sampling rate calculation
    f0 = 1/T;               % Fundamental frequency
    bw = 0.5*f0;            % Bandwidth

    % Normalized frequency range for Butterworth filter
    low_cutoff = max(0, (f0 - bw/2)/(fs/2));  % Ensure cutoff > 0
    high_cutoff = min(1, (f0 + bw/2)/(fs/2)); % Ensure cutoff < 1

    % Design 3rd-order Butterworth bandpass filter
    [b, a] = butter(3, [low_cutoff, high_cutoff], 'bandpass');

    % filtering
    c1 = filtfilt(b,a,c_stoch(:,1));  % zero-phase filtering
    c4 = filtfilt(b,a,c_stoch(:,2));  % zero-phase filtering

    % phase calculation
    phase1 = angle(hilbert(c1));
    phase4 = angle(hilbert(c4));

    %% plot

    figure
    hold on
    plot(t_stoch,c_stoch(:,1))
    plot(t_stoch,c_stoch(:,2))
    set(gca,'LineWidth',1.5,'FontSize',16)
    xlabel('\bf\boldmath$t$ / min','Interpreter','latex','FontSize',18)
    % ylabel('\bf\boldmath$\left[P_1\right]$ / molecule','interpreter','latex','Fontsize',18)
    ylabel('\bf\boldmath$c$ / molecule','interpreter','latex','Fontsize',18)
    legend('P_N1','P_N4','interpreter','none','Location','southwest')
    set(gca,'Xlim',[0,10000])
    % set(gca,'Ylim',[2000, 4000])
    box on
    exportgraphics(gcf,['output/noise_traj_' num2str(log10(p.sigma_inf)) '.pdf'],'resolution',300)

    % phase difference
    figure
    dphi = angle(exp(1i*(unwrap(phase1) - unwrap(phase4))));
    plot(t_stoch,dphi,'LineWidth',2)
    xlabel('\bf\boldmath$t$ / min','Interpreter','latex','FontSize',18)
    ylabel('\bf\boldmath$\Delta \phi$','interpreter','latex','Fontsize',18)
    set(gca,'LineWidth',1.5,'FontSize',16)
    set(gca,'Ylim',[-3.5,3.5])
    exportgraphics(gcf,['output/noise_PLV_' num2str(log10(p.sigma_inf)) '.pdf'],'resolution',300)


    % Lissajoules curves
    % get rid of transient
    t_stoch(1:start_pos-1) = [];
    c_stoch(1:start_pos-1,:) = [];

    figure
    hold on
    plot(c_stoch(:,1),c_stoch(:,2))
    set(gca,'LineWidth',1.5,'FontSize',16)
    xlabel('\bf\boldmath$[P_1]$ / molecule','Interpreter','latex','FontSize',18)
    % ylabel('\bf\boldmath$\left[P_1\right]$ / molecule','interpreter','latex','Fontsize',18)
    ylabel('\bf\boldmath$[P_4]$ / molecule','interpreter','latex','Fontsize',18)
    % set(gca,'Xlim',[0,4000])
    % set(gca,'Ylim',[0, 4000])
    box on
    exportgraphics(gcf,['output/noise_Lisa_' num2str(log10(p.sigma_inf)) '.pdf'],'resolution',300)
end