%% independent oscillators
% interaction between an independent repressilator and an independent
% CRISPRlator through a not gate, FM

clc
clear
close all

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

%% create the model (examples/coupled_mixed_not.mat)

clean_up_GRN
Ecoli = Cell('CRISPR');
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-8);
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-6);
set(getconfigset(Ecoli.data.Mobj),'Stoptime',1e4);
Ecoli = Ecoli.add_node('N1','type1');
Ecoli = Ecoli.add_node('N2','type1');
Ecoli = Ecoli.add_node('N3','type1');
Ecoli = Ecoli.add_regulator('Repression_in','N2','NOHILL','sgRNA_N1');
Ecoli = Ecoli.add_regulator('Repression_in','N3','NOHILL','sgRNA_N2');
Ecoli = Ecoli.add_regulator('Repression_in','N1','NOHILL','sgRNA_N3');
Ecoli.set('sgRNA_N3','InitialAmount',1.000000e+01,'N3');
Ecoli = Ecoli.add_node('N4','Elowitz');
Ecoli = Ecoli.add_node('N5','Elowitz');
Ecoli = Ecoli.add_node('N6','Elowitz');
Ecoli = Ecoli.add_regulator('Repression_in_TF','N5','HILL','P_N4');
Ecoli = Ecoli.add_regulator('Repression_in_TF','N6','HILL','P_N5');
Ecoli = Ecoli.add_regulator('Repression_in_TF','N4','HILL','P_N6');
Ecoli.set('P_N5','InitialAmount',1.000000e+03,'N5');
Ecoli = Ecoli.add_node('N7','mixed');
Ecoli = Ecoli.add_regulator('Repression_in','N7','NOHILL','sgRNA_N1');
Ecoli = Ecoli.add_regulator('Repression_in_TF','N7','HILL','P_N4');
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N7'}];

%% simulation settings

% get the model
Mobj = Ecoli.get_model();

% speed up the repressilator by faster protein degradation
set(sbioselect(Mobj.Parameters,'Name','d_P'),'Value',0.5)

% start and end of the simulation
t_end = 1e5;
t_min = 1e3;
% number of data points
n = 1e6;

configsetobj = getconfigset(Mobj);
set(get(configsetobj,'RuntimeOptions'),'StatesToLog',{'P_N1','P_N4','P_N7'})
set(configsetobj,'SolverType','sundials','Stoptime',t_end)
set(get(configsetobj,'SolverOptions'),'AbsoluteTolerance',1e-14,...
    'RelativeTolerance',1e-12)

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
T1 = calc_timeperiod(t,c(:,strcmp(names,'P_N1'))); % CRISPR (slow)
T4 = calc_timeperiod(t,c(:,strcmp(names,'P_N4'))); % repressilator (fast)

f1 = 1/T1;
f4 = 1/T4;

% FFT
x = c(:,strcmp(names,'P_N7'));

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
set(gca,'XLim',[f4-f1*5,f4+f1*5])
xlabel('\bf\boldmath$\nu$ / min$^{-1}$','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$|\mathcal{F}([P_7])|$','interpreter','latex','Fontsize',18) % / (molecule)$^2$min
set(gca,'LineWidth',2,'Fontsize',16)
% legend('FFT','[P_4]')

% add extra lines where we would expect the peaks
ylim = get(gca,'YLim');
hold on
plot([f4,f4],ylim,'r--','LineWidth',2)
for i = 1:4
    plot([f4-i*f1,f4-i*f1],ylim,'k--','LineWidth',2)
    plot([f4+i*f1,f4+i*f1],ylim,'k--','LineWidth',2)
end

hl = legend('FFT','f_4','f_4\pm n\cdotf_1');
hl.Direction = 'reverse';

% reverse order, so the FFT fill be above the others
hc = get(gca,'Children');
set(gca,'Children',hc(end:-1:1))
exportgraphics(hf,'output/fft.pdf','Resolution',300)

%% example simulation

% simulation time
set(getconfigset(Mobj),'Stoptime',2e3)

% simulate
simdata = sbiosimulate(Mobj);
[t,c,names] = getdata(simdata);

hf = figure;
hold on
plot(t,c(:,strcmp(names,'P_N4')),'LineWidth',2)
plot(t,c(:,strcmp(names,'P_N1')),'LineWidth',2)
xlabel('\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$c$ / molecule','interpreter','latex','Fontsize',18)
legend('[P_4]','[P_1]')
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(hf,'output/independent.pdf','Resolution',300)

hf = figure;
hold on
plot(t,c(:,strcmp(names,'P_N7')),'LineWidth',2,'Color','#EDB120')
xlabel('\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$c$ / molecule','interpreter','latex','Fontsize',18)
legend('[P_7]')
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(hf,'output/independent_output.pdf','Resolution',300)