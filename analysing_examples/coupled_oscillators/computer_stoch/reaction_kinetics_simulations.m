%% Reaction kinetics simulations
clc
clear
close all

%% Basic settings

% time period of light (periodic forcing)
T_light = 18; %7.45 13 18 25

constlight = false; % whether the LED light is constant
I_const = 30;       % the intensity of the constant light / %

%% Model construction

% load fitted parameters
load('/code/model/parameters_avr_A.mat','p','T_scale','A_scale')

Mobj = model();
add_parameters(Mobj,p);


%% Simulation

tsim = 3e2/T_scale; % end of simulation

n_data = 1e4; % number of data points in the simulations

% tolerance values
reltol = 1e-8;
abstol = 1e-10;

% simulation settings:
configsetobj = getconfigset(Mobj);
set(configsetobj,'SolverType','sundials','Stoptime',tsim,'MaximumWallClock',inf)
set(get(configsetobj,'SolverOptions'),'AbsoluteTolerance',abstol,...
    'RelativeTolerance',reltol,'OutputTimes',linspace(0,tsim,n_data));

%% resonance frequency

% accelerate
% sbioaccelerate(Mobj)

% which species is interesting to us:
configsetobj.RuntimeOptions.StatesToLog = {'x','y','z','mx','my','mz'};
% initial concentrations according to the dark state
load('initial_concentrations_dark.mat','c_fix')

if constlight == false
    % sin wave
    set(sbioselect(Mobj.Parameters,'Name','IA'),'Value',100)
    set(sbioselect(Mobj.Parameters,'Name','I0'),'Value',0)
else
    % deactivate light forcing
    Mobj.Rule(1).active = false;
    % constlight
    set(sbioselect(Mobj.Parameters,'Name','I'),'Value',I_const)
    % set these for the stochastic simulations
    set(sbioselect(Mobj.Parameters,'Name','IA'),'Value',0)
    set(sbioselect(Mobj.Parameters,'Name','I0'),'Value',30)
end

% set the time period of the light source
set(sbioselect(Mobj.Parameters,'Name','T'),'Value',T_light)

% deteministic simulation
[t,c,names] = sbiosimulate(Mobj);

%% stochastic simulation
mex CFLAGS='$CFLAGS -std=c++20 -Wall -pedantic -largeArrayDims -o3 -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' '-output' reaction_sim /code/reaction_kinetics_simulations/reaction_sim.cpp

% set the parameters for the stochastic simulations
parameters = {'I0','IA','T','k1','k2','k1x','nI','KI','n'};
for i = 1:numel(parameters)
    p.(parameters{i}) = get(sbioselect(Mobj.Parameters,'Name',parameters{i}),'Value');
end
p.omega = 5;
p.tsim = tsim;
p.n_step = n_data*10; % number of steps during the solution
p.n_data = p.n_step/round(p.n_step/n_data); % output data number
p.n_sim = 1000; % number of the averaged simulations

% run stochastic simulation
[c_stoch,t_stoch] = reaction_sim(p,c_fix);

% rescale
c_stoch = A_scale*c_stoch;
t_stoch = T_scale*t_stoch;
c = A_scale*c;
t = T_scale*t;

figure
set(gcf,'Position',[956   276   316   235])
plot(t,c(:,strcmp(names,'y')),'LineWidth',2)
hold on
plot(t_stoch,c_stoch,'LineWidth',2)
xlabel('\bf\boldmath$t$ (h)','FontSize',14,'interpreter','latex')
ylabel('\bf\boldmath$y$ (a.u.)','FontSize',14,'interpreter','latex')
set(gca,'FontSize',14,'LineWidth',1.5)
box on
if constlight == false
    exportgraphics(gcf,['/results/trajectory_' num2str(T_light) '.png'],'resolution',300)
else
    exportgraphics(gcf,['/results/trajectory_const_' num2str(I_const) '.png'],'resolution',300)
end

%% calculate phase space

tsim = 1e4; % end of simulation
tmin = 1e3; % start of data collection

n_data = 1e5; % number of data points

% simulation settings:
set(configsetobj,'Stoptime',tsim)
set(get(configsetobj,'SolverOptions'),'OutputTimes',[]);

% which species is interesting to us:
configsetobj.RuntimeOptions.StatesToLog = {'x','y','z'};

% simulation
tic
simdata = sbiosimulate(Mobj);
toc

% get simulation data
newSimData = resample(simdata,linspace(tmin,tsim,n_data));
[t,c,names] = getdata(newSimData);

% rescale
c = A_scale*c;
t = T_scale*t;

% phase(sub)space
figure
set(gcf,'Position',[956   276   316   235])
plot3(c(:,strcmp(names,'x')),c(:,strcmp(names,'y')),c(:,strcmp(names,'z')),'LineWidth',1)
xlabel('\bf\boldmath$x$ (a.u.)','FontSize',14,'interpreter','latex')
ylabel('\bf\boldmath$y$ (a.u.)','FontSize',14,'interpreter','latex')
zlabel('\bf\boldmath$z$ (a.u.)','FontSize',14,'interpreter','latex')
set(gca,'FontSize',14,'LineWidth',1.5)
grid on
box on
view(120,25)
if constlight == false
    exportgraphics(gcf,['/results/phase_space_' num2str(T_light) '.png'],'resolution',300)
else
    exportgraphics(gcf,['/results/phase_space_const_' num2str(I_const) '.png'],'resolution',300)
end

%% First-return map of the maxima 
% maxima return map

figure
set(gcf,'Position',[956   276   316   235])
pks = findpeaks(c(:,strcmp(names,'y')));
scatter(pks(1:end-1),pks(2:end),'filled')
xlabel('\bf\boldmath$y_{max}(n)$ (a.u.)','FontSize',18,'interpreter','latex')
ylabel('\bf\boldmath$y_{max}(n+1)$ (a.u.)','FontSize',18,'interpreter','latex')
set(gca,'FontSize',14,'LineWidth',1.5)
axis([20,160,20,160]*A_scale)
box on
if constlight == false
    exportgraphics(gcf,['/results/maxima_return_' num2str(T_light) '.png'],'resolution',300)
else
    exportgraphics(gcf,['/results/maxima_return_const_' num2str(I_const) '.png'],'resolution',300)
end

%% Fourier spectra
dt = t(2)-t(1);
ffty = fft(c(:,strcmp(names,'y')));
fs = 1/dt;
f = (0:length(ffty)-1)*fs/length(ffty);

figure
set(gcf,'Position',[956   276   316   235])
plot(f,abs(ffty),'LineWidth',2)
xlabel('\bf\boldmath$\nu$ (h$^{-1}$)','FontSize',14,'interpreter','latex')
ylabel('\bf Magnitude (a.u.)','FontSize',14,'interpreter','latex')
axis([0 1 0 1e5])
set(gca,'FontSize',16,'LineWidth',1.5)
box on
if constlight == false
    exportgraphics(gcf,['/results/fourier_' num2str(T_light) '.png'],'resolution',300)
else
    exportgraphics(gcf,['/results/fourier_const_' num2str(I_const) '.png'],'resolution',300)
end

%% calculate stroboscopic map
% resample simulation data in given time points
newSimData = resample(simdata,tmin:T_light:tsim);
[~,c,names] = getdata(newSimData);

% rescale
c = A_scale*c;

figure
set(gcf,'Position',[956   276   316   235])
scatter3(c(:,strcmp(names,'x')),c(:,strcmp(names,'y')),c(:,strcmp(names,'z')),'filled')
xlabel('\bf\boldmath$x$ (a.u.)','FontSize',14,'interpreter','latex')
ylabel('\bf\boldmath$y$ (a.u.)','FontSize',14,'interpreter','latex')
zlabel('\bf\boldmath$z$ (a.u.)','FontSize',14,'interpreter','latex')
set(gca,'FontSize',14,'LineWidth',1.5)
grid on
view(120,25)
axis([0,160,0,160,0,160]*A_scale)
box on
if constlight == false
    exportgraphics(gcf,['/results/stroboscopic_' num2str(T_light) '.png'],'resolution',300)
else
    exportgraphics(gcf,['/results/stroboscopic_const_' num2str(I_const) '.png'],'resolution',300)
end
