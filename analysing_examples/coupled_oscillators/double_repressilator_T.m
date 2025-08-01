%% calculate the time period of Elowitz repressilator in the function of nodes
clc
clear
close all

% we evaluate the simulations only after this time / min
start_time = 1e3;
stop_time = 1e4;

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

% number of nodes
% first loop
N1 = 1:2:15;
% second loop
N2 = 1:2:15;

% Time period / min
T = zeros(numel(N1,N2));

%% repressilator

tic
for i = 1:numel(N1)
    for j = 1:numel(N2)

        % we need nodes in the first loop
        if N1(i) >= 3
            % build the model
            Ecoli = build_double_repressilator('Elowitz',N1(i),N2(j));

            % simulation settings
            Ecoli.set('P_N1','InitialAmount',1.000000e+02,'N1');
            set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-8)
            set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-6)
            set(getconfigset(Ecoli.data.Mobj),'Stoptime',stop_time);
            Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];

            % run simulation
            Mobj = Ecoli.get_model();
            [t,c,~] = run_simulation(Mobj,'ode15s');

            % get rid of the initial transient part
            t_min = find(t>start_time,1,'first');
            t(1:t_min) = [];
            c(1:t_min,:) = [];

            % calculate time period
            T(i,j) = calc_timeperiod(t,c);
        else
            T(i,j) = nan;
        end

    end
    disp([int2str(i) '/' int2str(numel(N1))])
end
toc

% it is symmetric
if numel(N1)==numel(N2) && all(N1 == N2) && N1(1) == 1
    T(1,:) = T(:,1);
end

%% plot the data

!mkdir -p output

% [X,Y] = meshgrid(N2,N1);

% time period
figure
% hold on
% surf(X,Y,T)
h = imagesc(N2,N1,T);
% nan value should be transparent:
h.AlphaData = ~isnan(T);
axis equal
axis tight
axis xy
xticks(N2)
yticks(N1)
% view(2)
grid off
% box on
hc = colorbar;
% plot(fitobject)
xlabel('\bf\boldmath\#nodes$_1$','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath\#nodes$_2$','interpreter','latex','Fontsize',18)
zlabel('\bf\boldmath$T$ / min','interpreter','latex','Fontsize',18)
ylabel(hc,'\bf\boldmath$T$ / min','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(gcf,'output/T_double_repressilator.pdf')

%% 3D 
figure
% hold on
[X,Y] = meshgrid(N2,N1);
surf(X,Y,T)
xticks(N2)
yticks(N1)
hc = colorbar;
xlabel('\bf\boldmath\#nodes$_1$','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath\#nodes$_2$','interpreter','latex','Fontsize',18)
zlabel('\bf\boldmath$T$ / min','interpreter','latex','Fontsize',18)
ylabel(hc,'\bf\boldmath$T$ / min','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
drawnow
exportgraphics(gcf,'output/T_double_repressilator_3D.pdf')

%% plot linear functions

figure
hold on
plot(T(:,1),'o-','LineWidth',2)
plot(diag(T),'o-','LineWidth',2)
set(gca,'LineWidth',2,'Fontsize',16)
xlabel('\bf\boldmath\#nodes','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$T$ / min','interpreter','latex','Fontsize',18)
legend('single repressilator','double repressilator','Location','northwest')
exportgraphics(gcf,'output/T_double_repressilator_plot.pdf')