%% deeply coupled oscillators
% 4 node circuit containing two repressilator

clc
clear
close all

%% create the model (examples/node4_deeply_coupled.mat) with an extra protease

clean_up_GRN
Ecoli = Cell('Tomazou');
Ecoli = Ecoli.add_node('N1','type1');
Ecoli = Ecoli.add_node('N2','type1');
Ecoli = Ecoli.add_node('N3','type1');
Ecoli = Ecoli.add_node('N4','type1');
Ecoli = Ecoli.add_regulator('Repression_in','N2','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Repression_in','N3','HILL','P_N2');
Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL','P_N3');
Ecoli = Ecoli.add_regulator('Repression_in','N4','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Repression_in','N3','HILL','P_N4');
Ecoli = Ecoli.add_protease('N1','PROT1','type1');
Ecoli = Ecoli.add_protease('N2','PROT1','type1');
Ecoli = Ecoli.add_protease('N3','PROT1','type1');
Ecoli = Ecoli.add_protease('N4','PROT1','type1');
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-10);
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-8);
set(getconfigset(Ecoli.data.Mobj),'Stoptime',1e4);
Ecoli.data.Accelerate = true;
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N4'}];
Ecoli.set('P_N1','InitialAmount',1.000000e+03,'N1'); % Tomazou
Ecoli = Ecoli.add_protease('N4','PROT2','type1');
Ecoli.set('PROT2','InitialAmount',0.000000e+00,'PROT2');
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N2'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N3'}];

%% simulation settings

% get the model
Mobj = Ecoli.get_model();

% start and end of the simulation
t_end = 1e4;
t_min = 0.1*t_end;
% number of data points
n = 1e3;

% simulation settings:
configsetobj = getconfigset(Mobj);
set(configsetobj,'SolverType','sundials','Stoptime',t_end)
set(get(configsetobj,'SolverOptions'),'AbsoluteTolerance',1e-10,...
    'RelativeTolerance',1e-8)

disp('accelerate')
tic
sbioaccelerate(Mobj);
toc

% prescribed 'PROT1' concentrations
c_prot = [0:2:50,60:10:200];

% store data: T1, err, T4, err, T7, err. Amplitude in P_N1
data = zeros(length(c_prot),5);

%% simulate: the effect of 'PROT2'

tic
for i = 1:numel(c_prot)

    %  change the 'PROT2' concentration
    set(sbioselect(Mobj.Species,'Name','PROT2'),'Value',c_prot(i))
    
    % simulate
    simdata = sbiosimulate(Mobj);
    [t,c,names] = getdata(simdata);
    % get rid of the beginning
    start_pos = find(t>=t_min,1,'first');
    t(1:start_pos) = [];
    c(1:start_pos,:) = [];

    % calculate the time period of the independent oscillators
    [T1,err1] = calc_timeperiod(t,c(:,strcmp(names,'P_N1')));
    [T4,err4] = calc_timeperiod(t,c(:,strcmp(names,'P_N4')));

    % amplitude in P_N1
    A = max(c(:,strcmp(names,'P_N2')))-min(c(:,strcmp(names,'P_N2')));

    % save data
    data(i,:) = [T1,err1,T4,err4,A];
end
toc

%% plot the effect of 'PROT2' on the time period

!mkdir -p output

hf = figure;
plot(c_prot,data(:,3),'LineWidth',2)
set(gca,'LineWidth',2,'Fontsize',16)
ylabel('\bf\boldmath$T$ / min','interpreter','latex','Fontsize',18)
xlabel('\bf\boldmath$\left[\textrm{PROT}_2\right]$ / molecule','interpreter','latex','Fontsize',18)
exportgraphics(hf,'output/deeply_coupled_PROT_T.pdf','Resolution',300)

figure
plot(c_prot,data(:,5),'LineWidth',2)
set(gca,'LineWidth',2,'Fontsize',16)
ylabel('\bf\boldmath$A\left(\left[P_{N2}\right]\right)$ / molecule','interpreter','latex','Fontsize',18)
xlabel('\bf\boldmath$\left[\textrm{PROT}_2\right]$ / molecule','interpreter','latex','Fontsize',18)
exportgraphics(hf,'output/deeply_coupled_PROT_A.pdf','Resolution',300)

%% example simulation

c_prot = [0,50,500];

% start and end of the simulation
t_end = 2e3;
t_min = 1e3;

configsetobj = getconfigset(Mobj);
set(configsetobj,'Stoptime',t_end)

hf = figure;
ht = tiledlayout(3,1,'TileSpacing','tight');

for i = 1:numel(c_prot)

    % change the 'PROT2' concentration
    set(sbioselect(Mobj.Species,'Name','PROT2'),'Value',c_prot(i))

    % simulate
    simdata = sbiosimulate(Mobj);
    [t,c,names] = getdata(simdata);
    % get rid of the beginning
    start_pos = find(t>=t_min,1,'first');
    t(1:start_pos) = [];
    c(1:start_pos,:) = [];

    % plot the stuff
    nexttile
    hold on

    plot(t,c(:,strcmp(names,'P_N1'))/1e5,'LineWidth',2)
    plot(t,c(:,strcmp(names,'P_N2'))/1e5,'LineWidth',2)
    plot(t,c(:,strcmp(names,'P_N3'))/1e5,'LineWidth',2)
    plot(t,c(:,strcmp(names,'P_N4'))/1e5,'LineWidth',2)

    % title(['PROT2 = ' num2str(c_prot(i)) ' molecule'],'Fontsize',16)
    % no labels on the x axis
    if i ~= 3
        set(gca,'XTickLabel',[])
    end
    set(gca,'LineWidth',2,'Fontsize',16)

    if i == 1
        % to create a legend at the end so it will be on the top
        hax1 = gca;
    end

end
xlabel(ht,'\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel(ht,'\bf\boldmath$c$ / $10^5$ molecule','interpreter','latex','Fontsize',18)
legend(hax1,'[P_1]','[P_2]','[P_3]','[P_4]')
exportgraphics(hf,'output/deeply_coupled_PROT.pdf','Resolution',300)
