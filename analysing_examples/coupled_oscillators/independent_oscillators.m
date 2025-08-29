%% independent oscillators
% interaction between two independent repressilator and a not gate
% beat phenomenon

clc
clear
close all

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

%% create the model (examples/coupled_repressilators_not.mat)

clean_up_GRN
Ecoli = Cell('Tomazou');
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
Ecoli = Ecoli.add_protease('N1','PROT1','type1');
Ecoli = Ecoli.add_protease('N2','PROT1','type1');
Ecoli = Ecoli.add_protease('N3','PROT1','type1');
Ecoli = Ecoli.add_protease('N4','PROT2','type1');
Ecoli = Ecoli.add_protease('N5','PROT2','type1');
Ecoli = Ecoli.add_protease('N6','PROT2','type1');
Ecoli.data.Accelerate = true;
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1.000000e-08);
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1.000000e-10);
Ecoli = Ecoli.add_node('N7','type1');
Ecoli = Ecoli.add_regulator('Repression_in','N7','HILL','P_N1');
set(getconfigset(Ecoli.data.Mobj),'Stoptime',1.000000e+04);
Ecoli.set('P_N1','InitialAmount',1.000000e+03,'PROT1');
Ecoli.set('P_N4','InitialAmount',1.000000e+03,'PROT2');
Ecoli = Ecoli.add_regulator('Repression_in','N7','HILL','P_N4');
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N7'}];
Ecoli.set('PROT1','InitialAmount',1.640000e+02,'PROT1');

%% simulation settings

% get the model
Mobj = Ecoli.get_model();

% start and end of the simulation
t_end = 1e4;
t_min = 0.1*t_end;
% number of data points
n = 1e3;

% minimum number of peaks
n_minpeak = 3;

% number of workers
n_workers = 6;

% simulation settings:
configsetobj = getconfigset(Mobj);
set(configsetobj,'SolverType','sundials','Stoptime',t_end)
set(get(configsetobj,'SolverOptions'),'AbsoluteTolerance',1e-10,...
    'RelativeTolerance',1e-8)
% % output time points 
% configsetobj.SolverOptions.OutputTimes = linspace(t_min,t_end,n);

set(get(configsetobj,'RuntimeOptions'),'StatesToLog',{'P_N1','P_N4','P_N7'})

disp('accelerate')
tic
sbioaccelerate(Mobj);
toc

% prescribed 'PROT1' concentrations
c_prot = 10:5:200;
% c_prot should not be 'PROT2'
c_prot(c_prot==100) = [];

% store data: T1, err, T4, err, T7, err
data = zeros(length(c_prot),6);

%% simulate

% Get current parallel pool
p = gcp('nocreate');
% parallel computing
if isempty(p) || p.NumWorkers ~= n_workers
    delete(p);
    % start a new parallel pool
    parpool("Processes",n_workers);
end

tic
parfor i = 1:length(c_prot)

    t_end_par = t_end;
    t_min_par = t_min;
    n_par = n;

    % get the model
    Mobj_par = Mobj;

    %  change the 'PROT1' concentration
    set(sbioselect(Mobj_par.Species,'Name','PROT1'),'Value',c_prot(i))

    % number of peaks
    n_peaks = 0;

    % we will increase the length of the simulation if we have too few
    % peaks
    while n_peaks < n_minpeak || n_peaks == 0

        sbioaccelerate(Mobj_par);

        % [t,c,names] = sbiosimulate(Mobj);
        simdata = sbiosimulate(Mobj_par);
        [t,c,names] = getdata(simdata);
        % get rid of the beginning
        start_pos = find(t>=t_min,1,'first');
        t(1:start_pos) = [];
        c(1:start_pos,:) = [];

        % calculate the time period of the independent oscillators
        [T1,err1] = calc_timeperiod(t,c(:,strcmp(names,'P_N1')));
        [T4,err4] = calc_timeperiod(t,c(:,strcmp(names,'P_N4')));

        % resample th simulation
        newSimData = resample(simdata,linspace(t_min_par,t_end_par,n_par));
        dt = (t_end_par-t_min_par)/(n_par-1);
        [t_resampled,c_resampled,~] = getdata(simdata);
        % smooth over the time period of the fast oscillators
        a = smooth(c_resampled(:,strcmp(names,'P_N7')),100*ceil(max(T1,T4)/dt));
        % plot(t_resampled,a)
        [T7,err7,n_peaks] = calc_timeperiod(t_resampled,a,1/(abs(1/T1-1/T4))/4);

        % disp(1/abs(1/T1-1/T4))
        % disp(T7)
        % disp(n_peaks)

        % save data
        data(i,:) = [T1,err1,T4,err4,T7,err7];

        % set up a 10 times longer simulation, we will run it again if we
        % do not have enough peaks
        ratio = ceil(n_minpeak/n_peaks);
        if n_peaks == 0
            ratio = 10;
        end
        t_end_par = ratio*t_end_par;
        t_min_par = ratio*t_min_par;
        n_par = ratio*n_par;

        set(getconfigset(Mobj_par),'Stoptime',t_end_par)
    end
end
toc

%% plot the beats

!mkdir -p output

% base oscillator time period
T4 = mean(data(:,3));
% half position
pos = find(data(:,1)>T4,1,'last');

hf = figure;
fplot(@(x) 1./abs(1./T4-1./x),[min(data(:,1)),max(data(:,1))],'LineWidth',2)
hold on
errorbar(data(1:pos,1),data(1:pos,5),data(1:pos,6),'color',"#D95319",'LineWidth',2)
errorbar(data(pos+1:end,1),data(pos+1:end,5),data(pos+1:end,6),'color',"#D95319",'LineWidth',2)
set(gca,'LineWidth',2,'Fontsize',16)
% xlabel('\bf\boldmath$T$ / min','interpreter','latex','Fontsize',18)
xlabel('\bf Period of N1 / min','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$T_{beat}$ / min','interpreter','latex','Fontsize',18)
legend('Theoretical','Simulated')
exportgraphics(hf,'output/beat.pdf','Resolution',300)

%% example simulation

% change the 'PROT1' concentration
set(sbioselect(Mobj.Species,'Name','PROT1'),'Value',164)

set(get(configsetobj,'RuntimeOptions'),'StatesToLog',{'P_N1','P_N4','P_N7','protease_rate_PROT1','protease_rate_PROT2'})

% simulate
simdata = sbiosimulate(Mobj);
[t,c,names] = getdata(simdata);

%% plot

% example time series
tmin1 = 4500;
tmax1 = 5000;

tmin2 = 6750;
tmax2 = 7250;

hf = figure;
tlayout = tiledlayout(1,2,'TileSpacing','tight');
bgAx = axes(tlayout,'XTick',[],'YTick',[],'Box','off');
bgAx.Layout.TileSpan = [1 2];
ax1 = axes(tlayout);
hold on
plot(ax1,t,c(:,strcmp(names,'P_N4')),'LineWidth',2)
plot(ax1,t,c(:,strcmp(names,'P_N1')),'LineWidth',2)
xline(ax1,tmax1,'--','LineWidth',2);
ax1.Box = 'off';
xlim(ax1,[tmin1 tmax1])
set(ax1,'LineWidth',2,'Fontsize',16)
ytickformat('%.1f')

% gray area
yl = ylim;       % get current y-axis limits
fill([tmin1 tmax1 tmax1 tmin1], [yl(1) yl(1) yl(2) yl(2)], ...
     [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.3);  % gray with transparency

% Create second plot
ax2 = axes(tlayout);
ax2.Layout.Tile = 2;
hold on
plot(ax2,t,c(:,strcmp(names,'P_N4')),'LineWidth',2)
plot(ax2,t,c(:,strcmp(names,'P_N1')),'LineWidth',2)
xline(ax2,tmin2,'--','LineWidth',2);
ax2.YAxis.Visible = 'off';
ax2.Box = 'off';
xlim(ax2,[tmin2 tmax2])
set(ax2,'LineWidth',2,'Fontsize',16)

% Link the axes
linkaxes([ax1 ax2], 'y')

tlayout_children = get(tlayout,'children');
set(tlayout_children(end),'LineWidth',2,'Fontsize',16)

xlabel(tlayout,'\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel(tlayout,'\bf\boldmath$c$ / molecule','interpreter','latex','Fontsize',18)

% gray area
fill([tmin2 tmax2 tmax2 tmin2], [yl(1) yl(1) yl(2) yl(2)], ...
     [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.3);  % gray with transparency

legend('[P_4]','[P_1]')
set(gca,'ylim',yl)

% add axis break
% Create line
annotation(gcf,'line',[0.51 0.53],...
    [0.16 0.21],'LineWidth',2);

% Create line
annotation(gcf,'line',[0.52 0.54],...
    [0.16 0.21],'LineWidth',2);

exportgraphics(hf,'output/beat_oscillators.pdf','Resolution',300)

%%
hf = figure;
tlayout = tiledlayout(1,2,'TileSpacing','tight');
bgAx = axes(tlayout,'XTick',[],'YTick',[],'Box','off');
bgAx.Layout.TileSpan = [1 2];
ax1 = axes(tlayout);
hold on
plot(ax1,t,c(:,strcmp(names,'protease_rate_PROT1')),'LineWidth',2)
plot(ax1,t,c(:,strcmp(names,'protease_rate_PROT2')),'LineWidth',2)
xline(ax1,tmax1,'--','LineWidth',2);
ax1.Box = 'off';
xlim(ax1,[tmin1 tmax1])
set(ax1,'LineWidth',2,'Fontsize',16)
ytickformat('%.1f')

% gray area
yl = ylim;       % get current y-axis limits
fill([tmin1 tmax1 tmax1 tmin1], [yl(1) yl(1) yl(2) yl(2)], ...
     [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.3);  % gray with transparency

% Create second plot
ax2 = axes(tlayout);
ax2.Layout.Tile = 2;
hold on
plot(ax2,t,c(:,strcmp(names,'protease_rate_PROT1')),'LineWidth',2)
plot(ax2,t,c(:,strcmp(names,'protease_rate_PROT2')),'LineWidth',2)
xline(ax2,tmin2,'--','LineWidth',2);
ax2.YAxis.Visible = 'off';
ax2.Box = 'off';
xlim(ax2,[tmin2 tmax2])
set(ax2,'LineWidth',2,'Fontsize',16)

% Link the axes
linkaxes([ax1 ax2], 'y')

tlayout_children = get(tlayout,'children');
set(tlayout_children(end),'LineWidth',2,'Fontsize',16)

xlabel(tlayout,'\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel(tlayout,'\bf Protease activity\boldmath / min$^{-1}$','interpreter','latex','Fontsize',18)

% gray area
fill([tmin2 tmax2 tmax2 tmin2], [yl(1) yl(1) yl(2) yl(2)], ...
     [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.3);  % gray with transparency

legend('k_{PROT1}','k_{PROT2}')
set(gca,'ylim',yl)

% add axis break
% Create line
annotation(gcf,'line',[0.51 0.53],...
    [0.16 0.21],'LineWidth',2);

% Create line
annotation(gcf,'line',[0.52 0.54],...
    [0.16 0.21],'LineWidth',2);

exportgraphics(hf,'output/beat_oscillators_PROT.pdf','Resolution',300)
%%
% hf = figure;
% hold on
% plot(t,c(:,strcmp(names,'P_N4')),'LineWidth',2)
% plot(t,c(:,strcmp(names,'P_N1')),'LineWidth',2)
% xlabel('\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
% ylabel('\bf\boldmath$c$ / molecule','interpreter','latex','Fontsize',18)
% legend('[P_4]','[P_1]')
% set(gca,'LineWidth',2,'Fontsize',16)
% exportgraphics(hf,'output/beat_oscillators.pdf','Resolution',300)

% beat
hf = figure;
hold on
plot(t,c(:,strcmp(names,'P_N7')),'LineWidth',2,'Color','#EDB120')
xlabel('\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$c$ / molecule','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)

% gray area
yl = ylim;       % get current y-axis limits
fill([tmin1 tmax1 tmax1 tmin1], [yl(1) yl(1) yl(2) yl(2)], ...
     [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.3);  % gray with transparency
fill([tmin2 tmax2 tmax2 tmin2], [yl(1) yl(1) yl(2) yl(2)], ...
     [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.3);  % gray with transparency

% Create doublearrow
annotation(gcf,'doublearrow',[0.35 0.673214285714286],[0.599 0.6],'LineWidth',1.5);

% Create textbox
annotation(gcf,'textbox',...
    [0.455357142857141 0.579952383790705 0.166964281270547 0.119047616209303],...
    'String',{'$T_{beat}$'},...
    'FontSize',18,...
    'EdgeColor','none','Interpreter','latex');

legend('[P$_7$]','regions from ``b"','Location','northwest','interpreter','latex')

exportgraphics(hf,'output/beat_output.pdf','Resolution',300)

%% wavlet transform

% equidistant resampling
n_sample = 1e4;
newSimData = resample(simdata,linspace(t(1),t(end),n_sample));
[t,c,names] = getdata(newSimData);

% sampling frequency / Hz
Fs = 1/(t(2)-t(1))/60;

% plot wavelet transforms 
wavelet_plot(c(:,strcmp(names,'P_N1')),Fs)
exportgraphics(gcf,'output/beat_wavelet_N1.pdf','Resolution',300)

wavelet_plot(c(:,strcmp(names,'P_N4')),Fs)
exportgraphics(gcf,'output/beat_wavelet_N4.pdf','Resolution',300)

wavelet_plot(c(:,strcmp(names,'P_N7')),Fs)
exportgraphics(gcf,'output/beat_wavelet_N7.pdf','Resolution',300)


% calculating the Phase Locking Value (PLV)
cfs1 = cwt(c(:,strcmp(names,'P_N1')), 'amor',Fs); % 'amor' = analytic Morlet
[cfs4,freqs] = cwt(c(:,strcmp(names,'P_N4')), 'amor',Fs); % 'amor' = analytic Morlet

phase1 = angle(cfs1);
phase4 = angle(cfs4);

plv = abs(mean(exp(1i * (phase1- phase4)), 2));

figure
plot(freqs*1e3,plv,'LineWidth',2)
set(gca,'FontSize',16,'LineWidth',1.5)
xlabel('Frequency / mHz','FontSize',18)
ylabel('Phase Locking Value','FontSize',18)
set(gca,'YLim',[0,1])
exportgraphics(gcf,'output/beat_PLV.pdf','Resolution',300)

%% Hilbert transform

T = calc_timeperiod_avr_cross(t,c(:,strcmp(names,'P_N1')));

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
c1 = filtfilt(b,a,c(:,strcmp(names,'P_N1')));  % zero-phase filtering
c4 = filtfilt(b,a,c(:,strcmp(names,'P_N4')));  % zero-phase filtering

% phase calculation
phase1 = angle(hilbert(c1));
phase4 = angle(hilbert(c4));

% phase locking value
PLV = abs(mean(exp(1i * (phase1 - phase4))));
fprintf('PLV = %f\n',PLV);

% Instantaneous frequency
frequency1 = diff(unwrap(phase1)) / (2*pi*mean(diff(t)));
frequency4 = diff(unwrap(phase4)) / (2*pi*mean(diff(t)));

% frequency
figure
hold on
plot(t(1:end-1)/60/24,1e3*frequency1/60,'LineWidth',2)
plot(t(1:end-1)/60/24,1e3*frequency4/60,'LineWidth',2)
xlabel('\bf Time / days','interpreter','latex','Fontsize',18)
ylabel('\bf Frequency / mHz','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
ylim([0.1 0.15])
legend('N_1','N_4')
exportgraphics(gcf,'output/beat_hilbert_f.pdf','Resolution',300)

% phase difference
dphi = angle(exp(1i*(unwrap(phase1) - unwrap(phase4))));
figure
plot(t/60/24,dphi,'LineWidth',2)
xlabel('\bf Time / days','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$\Delta \phi$','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(gcf,'output/beat_hilbert_dphi.pdf','Resolution',300)

% Amplitude
amplitude1 = abs(hilbert(c1));
amplitude4 = abs(hilbert(c4));
figure
yyaxis left
plot(t/60/24,amplitude1,'LineWidth',2)
ylabel('\bf\boldmath Amplitude (N$_1$)','interpreter','latex','Fontsize',18)
yyaxis right
plot(t/60/24,amplitude4,'LineWidth',2)
xlabel('\bf Time / days','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath Amplitude (N$_4$)','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(gcf,'output/beat_hilbert_A.pdf','Resolution',300)

%% Functions
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