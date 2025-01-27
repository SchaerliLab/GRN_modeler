%% The Tomazou circuit
clc
close all
clear

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

%% build the model
clean_up_GRN

Ecoli = Cell('Tomazou');

%% simulation settings

set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-12);
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-10);
set(getconfigset(Ecoli.data.Mobj),'Stoptime',1e4);

%% Add nodes

Ecoli = Ecoli.add_node('R1');
Ecoli = Ecoli.add_node('R2');
Ecoli = Ecoli.add_node('R3');
Ecoli = Ecoli.add_node('G');

%% Add protease

Ecoli = Ecoli.add_protease('R1','C');
Ecoli = Ecoli.add_protease('R2','C');
Ecoli = Ecoli.add_protease('R3','C');
Ecoli = Ecoli.add_protease('G','L');


%% Add regulators

Ecoli = Ecoli.add_regulator('Repression_in','R1','R3');
Ecoli = Ecoli.add_regulator('Repression_in','R2','R1');
Ecoli = Ecoli.add_regulator('Repression_in','R3','R2');
Ecoli = Ecoli.add_regulator('Repression_in','G','R3');
Ecoli = Ecoli.add_regulator('Activation_in','R2','Y');
Ecoli.set('Y','Constant',true,'R2','Y');
Ecoli = Ecoli.add_regulator('Activation_out','R2','Y','I2');
Ecoli = Ecoli.add_regulator('Activation_in','G','U');
Ecoli.set('U','Constant',true,'G','U')
Ecoli = Ecoli.add_regulator('Activation_out','G','U','I1');


%% Initial conentrations


Ecoli.set('P_R2','InitialAmount',100,'R2')
Ecoli.set('C','InitialAmount',100,'C')
Ecoli.set('U','Value',100,'G','U')
Ecoli.set('Y','Value',100,'R2','Y')
Ecoli.set('L','InitialAmount',300,'L')

%% Parameters

% G copy number 
Ecoli.set('n_copy_G','Value',50,'G')
Ecoli.set('a1_G','Value',1000,'G')
Ecoli.set('K_molecule_G<-U','Value',50,'G','U')
Ecoli.set('n_copy_R2','Value',50,'R2')
Ecoli.set('a1_R2','Value',10000,'R2')
Ecoli.set('K_molecule_R2<-Y','Value',50,'R2','Y')

%% test the model

Mobj = Ecoli.get_model();

I1 = linspace(20,100,10); % uM
I2 = linspace(20,100,10); % uM

% minimum time when the concentration will be evaluated (get rid of
% transient)
t_min = 1e3;

% szimulacios beallitasok:
configsetobj = getconfigset(Mobj);
set(configsetobj,'SolverType','ode15s')

set(get(configsetobj,'RuntimeOptions'),'StatesToLog',{'P_G'})

sbioaccelerate(Mobj)

T = zeros(numel(I1),numel(I2));
A = zeros(numel(I1),numel(I2));

tic
for i = 1:numel(I1)
    for j = 1:numel(I2)

        set(sbioselect(Mobj.Species,'Name','I1'),'InitialAmount',I1(i));
        set(sbioselect(Mobj.Species,'Name','I2'),'InitialAmount',I2(j));
        
        % simulation
        [t,c,names] = sbiosimulate(Mobj);

        t_pos = find(t>t_min,1,'first');

        % calculate amplitude, time period
        T(i,j) = calc_timeperiod(t(t_pos:end),c(t_pos:end));
        A(i,j) = max(c(t_pos:end))-min(c(t_pos:end));
    end
end
toc


%% plot

!mkdir -p output
[X,Y] = meshgrid(I2,I1);

hf1 = figure;
surf(X,Y,A)
view(2)
% colormap winter
shading interp
% shading flat
hc = colorbar;
axis([X(1,1) X(1,end) Y(1,1) Y(end,1)])
% clim([minc maxc])
hc.Ruler.TickLabelFormat = '%.1f';
ylabel(hc,'\bf\boldmath$A$ / molecule','FontSize',18,'interpreter','latex')
% xtickformat('%2.0f')
% ytickformat('%2.0f')
set(gca,'FontSize',16)
set(hc,'FontSize',16)
set(gca,'LineWidth',1.5)
set(hc,'LineWidth',1.5)
xlabel('\bf\boldmath$I_{2}$ / $\mu$M','FontSize',18,'interpreter','latex')
ylabel('\bf\boldmath$I_{1}$ / $\mu$M','FontSize',18,'interpreter','latex')
% Show tick marks above plot
set(gca,'Layer','top')
grid off
box on
exportgraphics(hf1,'output/Tomazou_A.pdf','Resolution',300)

hf2 = figure;
surf(X,Y,T)
view(2)
shading interp
% shading flat
hc = colorbar;
axis([X(1,1) X(1,end) Y(1,1) Y(end,1)])
% clim([minc maxc])
% hc.Ruler.TickLabelFormat = '%.2f';
ylabel(hc,'\bf\boldmath$T$ / min','FontSize',18,'interpreter','latex')
% xtickformat('%2.0f')
% ytickformat('%2.0f')
set(gca,'FontSize',16)
set(hc,'FontSize',16)
set(gca,'LineWidth',1.5)
set(hc,'LineWidth',1.5)
xlabel('\bf\boldmath$I_{2}$ / $\mu$M','FontSize',18,'interpreter','latex')
ylabel('\bf\boldmath$I_{1}$ / $\mu$M','FontSize',18,'interpreter','latex')
% Show tick marks above plot
set(gca,'Layer','top')
grid off
box on
exportgraphics(hf2,'output/Tomazou_T.pdf','Resolution',300)

%% make graph
figure
Ecoli.data.inner_regulator_size = 6;
Ecoli.data.node_fontsize = 10;
Ecoli.make_graph();
axis off
exportgraphics(gcf,'output/tomazou_graph.pdf','Resolution',300)