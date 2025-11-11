%% independent oscillators
% interaction between an independent Stricker oscillator and an independent
% CRISPRlator through a not gate, AM

clc
clear
close all

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

% weaker repression for N6-|N3
repression_strength = 0.0001; % 1

%% create the model (examples/coupled_mixed_not_Stricker.mat)

clean_up_GRN
Ecoli = Cell('CRISPR');
% we need individuel kfdsd (connection to the DNA)
if repression_strength ~= 1
    target_parameter = sbioselect(Ecoli.data.regulator_models.Repression_in.Parameters,'Name','kfdsd');
    target_parameter.Notes = 'Individual';
end
solver = 'ode15s';
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-8)
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-6)
set(getconfigset(Ecoli.data.Mobj),'MaximumWallClock',20)
Ecoli = Ecoli.add_node('N1','Elowitz');
Ecoli = Ecoli.add_node('N2','Elowitz');
Ecoli = Ecoli.add_regulator('Activation_in_TF','N2','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Activation_in_TF','N1','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Repression_in_TF','N2','HILL','P_N2');
Ecoli = Ecoli.add_regulator('Repression_in_TF','N1','HILL','P_N2');
Ecoli = Ecoli.add_node('N3','type1');
Ecoli = Ecoli.add_node('N4','type1');
Ecoli = Ecoli.add_node('N5','type1');
Ecoli = Ecoli.add_node('N6','mixed');
Ecoli = Ecoli.add_regulator('Repression_in','N4','NOHILL','sgRNA_N3');
Ecoli = Ecoli.add_regulator('Repression_in','N5','NOHILL','sgRNA_N4');
Ecoli = Ecoli.add_regulator('Repression_in','N3','NOHILL','sgRNA_N5');
Ecoli = Ecoli.add_regulator('Repression_in','N6','NOHILL','sgRNA_N3');
Ecoli = Ecoli.add_regulator('Repression_in_TF','N6','HILL','P_N2');
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N2'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N3'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N6'}];
set(getconfigset(Ecoli.data.Mobj),'StopTime',2.000000e+03);



%% simulation settings

% Higher nonlinearity for oscillation
set(sbioselect(Ecoli.data.Mobj,'Name','n_molecule_N2<-N1'),'Value',6.000000e+00);
set(sbioselect(Ecoli.data.Mobj,'Name','n_molecule_N1|-N2'),'Value',6.000000e+00);

% get the model
Mobj = Ecoli.get_model();

% speed up the repressilator by faster protein degradation
set(sbioselect(Mobj.Parameters,'Name','d_P'),'Value',0.2)

% changing the repression strength
if repression_strength ~= 1
    target_parameter = sbioselect(Mobj,'Name','kfdsd_N6NOHILL|-N3');
    target_parameter.Value = target_parameter.Value*repression_strength;
end

% start and end of the simulation
t_end = 1e5;
t_min = 1e3;
% number of data points
n = 1e6;

configsetobj = getconfigset(Mobj);
set(get(configsetobj,'RuntimeOptions'),'StatesToLog',{'P_N2','P_N3','P_N6'})
set(configsetobj,'SolverType','sundials','Stoptime',t_end)
set(get(configsetobj,'SolverOptions'),'AbsoluteTolerance',1e-12,...
    'RelativeTolerance',1e-10)

% set(get(getconfigset(Mobj),'SolverOptions'),'OutputTimes',linspace(t_min,t_end,n));

disp('accelerate')
tic
sbioaccelerate(Mobj);
toc

% simulate
tic
simdata = sbiosimulate(Mobj);
toc
% resample the simulation
newSimData = resample(simdata,linspace(t_min,t_end,n));
[t,c,names] = getdata(newSimData);


%% Fourier transform

% the basic oscillator frequencies
T2 = calc_timeperiod(t,c(:,strcmp(names,'P_N2'))); % CRISPR (slow)
T3 = calc_timeperiod(t,c(:,strcmp(names,'P_N3'))); % repressilator (fast)

f2 = 1/T2;
f3 = 1/T3;

% FFT
x = c(:,strcmp(names,'P_N6'));

Y = fft(x);
L = length(t);
Fs = 1/(t(2)-t(1));

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;

%% plot FFT
hf = figure;
hp = plot(f(2:end),P1(2:end),'LineWidth',2);
set(gca,'XLim',[f2-f3*4.5,f2+f3*4.5])
xlabel('\bf\boldmath$\nu$ / min$^{-1}$','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$|\mathcal{F}([P_6])|$','interpreter','latex','Fontsize',18) % / (molecule)$^2$min
set(gca,'LineWidth',2,'Fontsize',16)
% legend('FFT','[P_4]')

% add extra lines where we would expect the peaks
ylim = get(gca,'YLim');
hold on
plot([f2,f2],ylim,'r--','LineWidth',2)
for i = 1:4
    plot([f2-i*f3,f2-i*f3],ylim,'k--','LineWidth',2)
    plot([f2+i*f3,f2+i*f3],ylim,'k--','LineWidth',2)
end

hl = legend('FFT','f_2','f_2\pm n\cdotf_3');
hl.Direction = 'reverse';

% reverse order, so the FFT fill be above the others
hc = get(gca,'Children');
set(gca,'Children',hc(end:-1:1))
exportgraphics(hf,'output/fft_Stricker.pdf','Resolution',300)

%% example simulation

% simulation time
set(getconfigset(Mobj),'Stoptime',2e3)

% simulate
simdata = sbiosimulate(Mobj);
[t,c,names] = getdata(simdata);

hf = figure;
hold on
% yyaxis left
plot(t,c(:,strcmp(names,'P_N3')),'LineWidth',2)
% ylabel('\bf\boldmath$c$ / molecule','interpreter','latex','Fontsize',18)
% yyaxis right
plot(t,c(:,strcmp(names,'P_N2')),'LineWidth',2)
xlabel('\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$c$ / molecule','interpreter','latex','Fontsize',18)
legend('[P_3]','[P_2]')
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(hf,'output/independent_Stricker.pdf','Resolution',300)

hf = figure;
hold on
plot(t,c(:,strcmp(names,'P_N6')),'LineWidth',2,'Color','#EDB120')
xlabel('\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$c$ / molecule','interpreter','latex','Fontsize',18)
legend('[P_6]')
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(hf,'output/independent_output_Stricker.pdf','Resolution',300)

%% Wavelet transform


% equidistant resampling
n_sample = 1e4;
newSimData = resample(simdata,linspace(t(1),t(end),n_sample));
[t,c,names] = getdata(newSimData);

% sampling frequency / Hz
Fs = 1/(t(2)-t(1))/60;

% plot wavelet transforms 
wavelet_plot(c(:,strcmp(names,'P_N2')),Fs)
exportgraphics(gcf,'output/AM_stricker_wavelet_N2.pdf','Resolution',300)

wavelet_plot(c(:,strcmp(names,'P_N3')),Fs)
exportgraphics(gcf,'output/AM_stricker_wavelet_N3.pdf','Resolution',300)

wavelet_plot(c(:,strcmp(names,'P_N6')),Fs)
exportgraphics(gcf,'output/AM_stricker_wavelet_N6.pdf','Resolution',300)


% calculating the Phase Locking Value (PLV)
cfs2 = cwt(c(:,strcmp(names,'P_N2')), 'amor',Fs); % 'amor' = analytic Morlet
[cfs3,freqs] = cwt(c(:,strcmp(names,'P_N3')), 'amor',Fs); % 'amor' = analytic Morlet

phase1 = angle(cfs2);
phase4 = angle(cfs3);

plv = abs(mean(exp(1i * (phase1- phase4)), 2));

figure
plot(freqs*1e3,plv,'LineWidth',2)
set(gca,'FontSize',16,'LineWidth',1.5)
xlabel('Frequency / mHz','FontSize',18)
ylabel('Phase Locking Value','FontSize',18)
set(gca,'YLim',[0,1])
exportgraphics(gcf,'output/AM_stricker_wavelet_PLV.pdf','Resolution',300)


function [] = wavelet_plot(c,Fs)
figure
cwt(c, 'amor',Fs)
set(gca,'FontSize',16,'LineWidth',1.5)
% cb = findall(gcf, 'Type', 'ColorBar');
cb = colorbar;
ylabel(cb,'Magnitude','FontSize',18)
cb.Ruler.Exponent = 0;
cb.Ruler.TickLabelFormat = '%.0e';  % or %g for compactset(cb,'FontSize',18,'LineWidth',1.5)
ylabel('Frequency / mHz','FontSize',18)
xlabel('Time / days','FontSize',18)
title ''

% shring th width a little
gcapos = get(gca,'Position'); % [left, bottom, width, height]
gcapos(3) = 0.925*gcapos(3);
% gcapos(4) = 0.975*gcapos(4);
set(gca, 'Position', gcapos);
end