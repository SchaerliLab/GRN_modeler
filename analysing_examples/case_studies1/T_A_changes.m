%% GNR simulation
clc
clear
close all

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

%% build the model
clean_up_GRN
Ecoli = Cell('CRISPR');

% use common parameters
set(sbioselect(Ecoli.data.node_models.type1.Parameters,'Name','a1'),'Notes','Common')

set(getconfigset(Ecoli.data.Mobj),'Stoptime',5e3);
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-8);
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-6);
Ecoli = Ecoli.add_node('N1','type1');
Ecoli = Ecoli.add_node('N2','type1');
Ecoli = Ecoli.add_node('N3','type1');
Ecoli = Ecoli.add_regulator('Repression_in','N2','NOHILL','sgRNA_N1');
Ecoli = Ecoli.add_regulator('Repression_in','N3','NOHILL','sgRNA_N2');
Ecoli = Ecoli.add_regulator('Repression_in','N1','NOHILL','sgRNA_N3');
Ecoli.set('sgRNA_N3','InitialAmount',1.000000e+01,'N3');
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];
% Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N2'}];
% Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N3'}];

Mobj = Ecoli.get_model();


% clear Mobj;
% const_dCas = false; % new model: [dCas]total can change 
% Mobj = model(modelnev,const_dCas);
% Mobj = set_ara(Mobj,n_edges,n_ara); % just for other settings
tmin = 1e4; % first time data in the outut
tsim = 1e4; % end of simulation / min
nt = 1e4; % number of time data in the output

% const_dCas = false; % new model: [dCas]total can change 
% if const_dCas==false
%     kbdCas = 1e-3;
%     set(sbioselect(Mobj.Parameters,'Name','kbdCas'),'Value',kbdCas)
%     % kfdCas is set by a rule at the beginning of the simulation
% %     set(sbioselect(Mobj.Parameters,'Name','kfdCas'),'Value',get(sbioselect(Mobj.Species,'Name','dCas'),'Value')*kbdCas)
% end

% simulation settings:
configsetobj = getconfigset(Mobj);
set(configsetobj,'SolverType','sundials','Stoptime',tmin+tsim,'MaximumWallClock',inf)
set(get(configsetobj,'SolverOptions'),'AbsoluteTolerance',1e-8,...
    'RelativeTolerance',1e-6)%,'OutputTimes',linspace(tmin,tmin+tsim,nt))
set(configsetobj.CompileOptions, 'DimensionalAnalysis', false);

configsetobj.RuntimeOptions.StatesToLog = {'P_N1'};

%% global sensitivity analysis

addobservable(Mobj,'Time period','calc_T2(P_N1,time,0.1)');%'calc_T(P_N1,time,20*60,0.1)');
addobservable(Mobj,'Amplitude','max(P_N1)-min(P_N1)');

modelParamNames = {'a1','dCas','k_P','d_P','krds','kfds','kfdsd'};
% bounds = [];
outputName = {'Time period','Amplitude'};
rng('default');

arany = 2;
bounds = zeros(length(modelParamNames),2);
for i=1:length(modelParamNames)
    par = get(sbioselect(Mobj,'Name',modelParamNames{i}),'Value');
    bounds(i,:) = [par/arany,par*arany];
end

%,'Bounds',bounds
tic
sobolResults = sbiosobol(Mobj,modelParamNames,outputName,'Bounds',bounds,'ShowWaitBar',true,'NumberSamples',1000,'Useparallel',false,'Accelerate',true,'OutputTimes',linspace(tmin,tmin+tsim,nt));
toc

!mkdir -p output
% Show the mean model response, the simulation results, and a shaded region covering 90% of the simulation results.
plotData(sobolResults,ShowMedian=true,ShowMean=false);
exportgraphics(gcf,'output/sobol1.pdf','Resolution',300)

% Plot the time course of the first- and total-order Sobol indices.
h = plot(sobolResults);
% Resize the figure.
h.Position(:) = [100 100 1280 800];
exportgraphics(gcf,'output/sobol2.pdf','Resolution',300)

% Display the magnitudes of the sensitivities in a bar plot. Darker colors mean that those values occur more often over the whole time course.
bar(sobolResults);
% a = get(get(gcf,'Children'),'Children');
% set(a(4),'XLim',[-.1,1.4])
ht = get(get(gcf,'children'),'children');
set(ht(2),'LineWidth',2,'FontSize',16)
yticklabels(ht(2),[])
ylabel(ht(2),[])
set(ht(4),'LineWidth',2,'FontSize',16)
ylabel(ht(4),'\bf Sensitivity Input')
set(ht(1),'Location','southeast','FontSize',12)
set(ht(3),'Location','southeast','FontSize',12)
title(get(gcf,'children'),'\bf Sensitivity Output','FontSize',18)
title(ht(4),'Time period','FontSize',18)
title(ht(2),'Amplitude','FontSize',18)
set(get(gcf,'children'),'TileSpacing','compact')
xlabel(ht(4),'\bf Sobol Index')
xlabel(ht(2),'\bf Sobol Index')
exportgraphics(gcf,'output/sobol3.pdf','Resolution',300)

%% Simulation

% tmin = 1e4; % first time data in the outut
% tsim = 1e4; % end of simulation / min
% nt = 1e4; % number of time data in the output
% T_approx = 60*20; % for autocorrelation function lag time
% 
% % simulation settings:
% configsetobj = getconfigset(Mobj);
% set(configsetobj,'SolverType','sundials','Stoptime',tmin+tsim,'MaximumWallClock',inf)
% set(get(configsetobj,'SolverOptions'),'AbsoluteTolerance',1e-8,... % 'AbsoluteToleranceScaling',false
%     'RelativeTolerance',1e-6) % ,'OutputTimes',linspace(tmin,tmin+tsim,nt))

% configsetobj.RuntimeOptions.StatesToLog = {'P_N1'};


disp('accelerating...')
tic
sbioaccelerate(Mobj)
% accelerate(Mobjexp)
toc

ncas = 50;
nkd = 50;
cas0 = get(sbioselect(Mobj,'Name','dCas'),'Value');
k_P = get(sbioselect(Mobj,'Name','k_P'),'Value');
[CAS,KP] = meshgrid(linspace(cas0*0.25,cas0,ncas),linspace(k_P*0.5,k_P*2,nkd));
T = zeros(nkd,ncas);
A = zeros(nkd,ncas);

tic
for i = 1:nkd
    for j = 1:ncas
        set(sbioselect(Mobj,'Name','dCas'),'Value',CAS(i,j))
        set(sbioselect(Mobj,'Name','k_P'),'Value',KP(i,j))
        [t,c,names] = sbiosimulate(Mobj);
        T(i,j) = calc_timeperiod_avr_cross(t,c(:,strcmp(names,'P_N1')));% calc_T(c(:,strcmp(names,'P_N1')),t,T_approx,0.1);
        A(i,j) = (max(c(:,strcmp(names,'P_N1')))-min(c(:,strcmp(names,'P_N1'))));%/max(c(:,strcmp(names,'P_N1')));
    end
end
toc

T = T/60; % h

figure;
surf(CAS,KP,T)
view(2)
% colormap winter
shading interp
% shading flat
hc = colorbar;
% title(['t = ' num2str(t,'%2.1f')])
% axis([X(1,1) X(1,end) Y(1,1) Y(end,1)])
% clim([minc maxc])
% hc.Ruler.TickLabelFormat = '%.2f';
ylabel(hc,'\bf\boldmath$T$ / h','FontSize',18,'interpreter','latex')
% xtickformat('%2.0f')
% ytickformat('%2.0f')
set(gca,'FontSize',16)
set(hc,'FontSize',16)
set(gca,'LineWidth',1.5)
set(hc,'LineWidth',1.5)
ylabel('\bf\boldmath$k_{P}$  / (molecule/min)','FontSize',18,'interpreter','latex')
xlabel('\bf\boldmath$c_{cas}$ / nM','FontSize',18,'interpreter','latex')
% Show tick marks above plot
set(gca,'Layer','top')
grid off
box on
axis([min(CAS,[],'all') max(CAS,[],'all') min(KP,[],'all') max(KP,[],'all')])
% set(gca,'xaxisLocation','top')
% set(gca,'YScale','log')
exportgraphics(gcf,'output/T.pdf','Resolution',300)

figure;
surf(CAS,KP,A/1e3)
view(2)
% colormap winter
shading interp
% shading flat
hc = colorbar;
% title(['t = ' num2str(t,'%2.1f')])
% axis([X(1,1) X(1,end) Y(1,1) Y(end,1)])
% clim([minc maxc])
% hc.Ruler.TickLabelFormat = '%.2f';
ylabel(hc,'\bf \boldmath$A$ / $10^3$ molecule','FontSize',18,'interpreter','latex')
% xtickformat('%2.0f')
% ytickformat('%2.0f')
set(gca,'FontSize',16)
set(hc,'FontSize',16)
set(gca,'LineWidth',1.5)
set(hc,'LineWidth',1.5)
ylabel('\bf\boldmath$k_{P}$  / (molecule/min)','FontSize',18,'interpreter','latex')
xlabel('\bf\boldmath$c_{cas}$ / nM','FontSize',18,'interpreter','latex')
% Show tick marks above plot
set(gca,'Layer','top')
grid off
box on
axis([min(CAS,[],'all') max(CAS,[],'all') min(KP,[],'all') max(KP,[],'all')])
% set(gca,'xaxisLocation','top')
% set(gca,'YScale','log')
exportgraphics(gcf,'output/A.pdf','Resolution',300)

