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
% set the dcas:sgRNA + DNA rate constant to individual
set(sbioselect(Ecoli.data.regulator_models.Repression_in.Parameters,'Name','kfdsd'),'Notes','Individual')

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

% extra inter oscillator repression
Ecoli = Ecoli.add_regulator('Repression_in','N4','NOHILL','sgRNA_N1');

Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N2'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N3'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N4'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N5'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N6'}];

Ecoli.set('dCas','InitialAmount',2.0000e+02,'N6'); % change dCas here
Ecoli.set('a1_N6','Value',1.2,'N6');
Ecoli.set('a1_N2','Value',1.1,'N2');

set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-10);
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-8);
set(getconfigset(Ecoli.data.Mobj),'Stoptime',1e5);
Ecoli.data.Accelerate = true;

set(get(getconfigset(Ecoli.data.Mobj),'CompileOptions'), 'UnitConversion', false);
set(get(getconfigset(Ecoli.data.Mobj),'CompileOptions'), 'DimensionalAnalysis', false);


solver = 'sundials';%'adaptivesa';

% addparameter(Ecoli.data.Mobj,'K_activity',1e4);

kfdsd_original = get(sbioselect(Ecoli.data.regulator_models.Repression_in.Parameters,'Name','kfdsd'),'Value');

% get the model
Mobj = Ecoli.get_model();

%% create C++ files for stochastic simulations

% copy the model when we modify the names
Mobj2 = copyobj(Mobj);

% convert reversible reactions to irreversible
Mobj2 = convert2irrev(Mobj2);
Mobj2 = correct_modifiers(Mobj2);

% === reactions.h ===
% create the reaction header file
write_reaction_header(Mobj2)

% === parameters ===
% create parameters.h header file
write_parameters_header(Mobj2)
% create parameters.cpp
write_parameters_cpp(Mobj2)


%% create the mex file

mex CXXFLAGS="\$CXXFLAGS -std=c++20 -Wall -pedantic -O3 -fopenmp" ...
    LDFLAGS="\$LDFLAGS -fopenmp" ...
    -output reaction_sim reaction_sim.cpp parameters.cpp

%% parameter settings

% parameters into the p structure
for i = 1:numel(Mobj.Parameters)
    % change '<-' to '_A_' and '|-' to '_R_'
    parname = Mobj.Parameter(i).Name;
    parname = strrep(parname, '<-', '_A_');
    parname = strrep(parname, '|-', '_R_');
    p.(parname) = Mobj.Parameters(i).Value;
end
% initial concentrations
c0 =( zeros(numel(Mobj.Species),1));
for i = 1:numel(Mobj.Species)
    c0(i) = Mobj.Species(i).Value;
end

% position of the proteins in the concentration vector
% we are following only these species
protein_ID = find(contains({Mobj.Species.Name},'P_N')).'-1; % -1: C matrix indexing

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

%% the effect of the repression strength between the oscillators and the noise

% interoscillator repression and noise stregth
repression_strength = [0,0.01,0.1];%logspace(-4,-1,4);
sigma_inf = logspace(-2,-4,6); % 6

[S,R] = meshgrid(sigma_inf,repression_strength);
C = zeros(size(R));
PLV = zeros(size(R));

tic
for i = 1:size(R,1)
    for j = 1:size(R,2)

        % weak interoscillator repression
        p.kfdsd_N4NOHILL_R_N1 = kfdsd_original*R(i,j);
        % noise strength
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

% % max abs correlation
% figure
% hold on
% for i = 1:numel(repression_strength)
%     plot(sigma_inf,C(i,:),'-o','LineWidth',2)
% end
% set(gca,'LineWidth',1.5,'Fontsize',16)
% xlabel('\boldmath$\sigma_{\infty}$','Interpreter','latex','FontSize',18)
% ylabel('\bf correlation strength','Interpreter','latex','FontSize',18)
% hl = legend(num2str(repression_strength.'),'Location','southwest');
% title(hl,'\boldmath$\beta$','Interpreter','latex')
% ytickformat('%.1f')
% set(gca,'YLim',[0,1],'XScale','log')
% exportgraphics(gcf,'output/coupling_noise_repr.pdf','resolution',300)

% Phase Locking Value
figure
hold on
for i = 1:numel(repression_strength)
    plot(sigma_inf,PLV(i,:),'-o','LineWidth',2)
end
set(gca,'LineWidth',1.5,'Fontsize',16)
xlabel('\boldmath$\sigma_{\infty}$','Interpreter','latex','FontSize',18)
ylabel('\bf Phase Locking Value','Interpreter','latex','FontSize',18)
hl = legend(num2str(repression_strength.'),'Location','southwest');
title(hl,'\boldmath$\beta$','Interpreter','latex')
ytickformat('%.1f')
set(gca,'YLim',[0,1],'XScale','log')
exportgraphics(gcf,'output/coupling_noise_PLV_repr.pdf','resolution',300)

%% Example stochastic simulations

p.t_sim = 1e5;%t(end); % end of simulation
n_data = p.t_sim ;% output data number
p.n_step = n_data*1e2;%*10; % number of steps during the solution
p.n_data = p.n_step/round(p.n_step/n_data); % output data number

% ratio of the transient part
transient = 0.15;
start_pos = round(transient*n_data+1);

% change the strength of the repression between the oscillators
for beta = [0,0.01,0.1]

    clear reaction_sim % clear static variables in functions

    p.kfdsd_N4NOHILL_R_N1 = kfdsd_original*beta;
    p.sigma_inf = 1e-2;

    % run stochastic simulation
    tic
    [c_stoch,t_stoch] = reaction_sim(p,c0,protein_ID([1,4:6]));
    toc

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
    set(gca,'Xlim',[0,1e4])
    % set(gca,'Ylim',[2000, 4000])
    box on
    exportgraphics(gcf,['output/noise_traj_' num2str(log10(p.sigma_inf)) '_' num2str(beta) '.pdf'],'resolution',300)

    % phase difference
    figure
    dphi = angle(exp(1i*(unwrap(phase1) - unwrap(phase4))));
    plot(t_stoch,dphi,'LineWidth',2)
    xlabel('\bf\boldmath$t$ / min','Interpreter','latex','FontSize',18)
    ylabel('\bf\boldmath$\Delta \phi$','interpreter','latex','Fontsize',18)
    set(gca,'LineWidth',1.5,'FontSize',16)
    set(gca,'Ylim',[-3.5,3.5])
    exportgraphics(gcf,['output/noise_PLV_' num2str(log10(p.sigma_inf)) '_' num2str(beta) '.pdf'],'resolution',300)


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
    exportgraphics(gcf,['output/noise_Lisa_' num2str(log10(p.sigma_inf)) '_' num2str(beta) '.pdf'],'resolution',300)
end