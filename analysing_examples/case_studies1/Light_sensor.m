%% Light and Ara inducible circuit
clc
close all
clear

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

clean_up_GRN

% open, closed or 1 node circuit
circuit = 'closed'; % 'open' 'closed' 'N1'
% relative leakage strength
leak = 0.01;
% amplitude of light (in %)
A = 100;
% arabinose concentration
ara = 0*1e-3;

%% build the circuit
Ecoli = Cell('CRISPR');
solver = 'ode15s';
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-8)
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-6)
set(getconfigset(Ecoli.data.Mobj),'MaximumWallClock',10)
Ecoli = Ecoli.add_node('N1','type1');
Ecoli = Ecoli.add_regulator('Light_Ara','N1','HILL','R1');
Ecoli = Ecoli.add_regulator('Act_direct','N1','HILL','R1','HILL_light','Light');
Ecoli = Ecoli.add_regulator('Act_direct','N1','HILL','R1','HILL_ara','Ara');
set(sbioselect(Ecoli.data.Mobj,'Name','Light'),'Constant',false);
addrule(Ecoli.data.Mobj,'Light = A*(1+sign(sin(2*pi/T*time)))/2','repeatedAssignment','Name','SquareWave');
addparameter(Ecoli.data.Mobj,'A','Value',A,'Units','micromolarity','Notes','Individual');
addparameter(Ecoli.data.Mobj,'T','Value',1.440000e+03,'Units','minute','Notes','Individual');
set(get(getconfigset(Ecoli.data.Mobj),'CompileOptions'),'Dimensionalanalysis',false)
set(getconfigset(Ecoli.data.Mobj),'StopTime',1e3);
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];
switch circuit
    case {'open','closed'}
        Ecoli = Ecoli.add_node('N2','type1');
        Ecoli = Ecoli.add_regulator('Repression_in','N2','NOHILL','sgRNA_N1');
        Ecoli = Ecoli.add_node('N3','type1');
        Ecoli = Ecoli.add_regulator('Repression_in','N3','NOHILL','sgRNA_N2');
        Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N2'}];
        Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N3'}];
end
switch circuit
    case 'closed'
        Ecoli = Ecoli.add_regulator('Repression_in','N1','NOHILL','sgRNA_N3');
end

% follow the light:
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'Light'}];

set(getconfigset(Ecoli.data.Mobj),'StopTime',5.760e3);
Ecoli.set('sgRNA_N1','InitialAmount',1.000000e+01,'N1');

% Ecoli.make_graph();

% changed parameters against oscillation:
set(sbioselect(Ecoli.data.Mobj,'Name','dilution'),'Value',5.000000e-03); % dilution is slover on solid surface
set(sbioselect(Ecoli.data.Mobj,'Name','krdsd'),'Value',0.7762); % Javier

% leakage
set(sbioselect(Ecoli.data.Mobj,'Name','leak'),'Value',leak);

% arabionose
Ecoli.set('Ara','InitialAmount',ara,'N1','HILL','<-','R1','HILL_ara','<-','Ara');

% strength of the light inducible node: 2.36
Ecoli.set('a1_N1','Value',2.360000e+00,'N1');

%% Run simulation
Mobj = Ecoli.get_model();
set(getconfigset(Mobj),'SolverType',solver);
[t,c,names] = run_simulation(Mobj,solver);

hf = figure;
yyaxis left
hold on
switch circuit
    case 'N1'
        plot(t,c(:,1),'LineWidth',2,'Color',"#77AC30")
    case {'open','closed'}
        plot(t,c(:,1),'LineWidth',2,'Color',"#A2142F")
        plot(t,c(:,2),'LineWidth',2,'Color',"#0072BD")
        plot(t,c(:,3),'LineWidth',2,'Color',"#77AC30")
end
xlabel('\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$c$ / molecule','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
yyaxis right
plot(t,c(:,4),'LineWidth',2,'Color',"#EDB120")
ax = gca;
ax.YAxis(2).Color = "#EDB120";
legend(Ecoli.data.StatesToLog,'Interpreter','none')
ylabel('\bf\boldmath$I$ / \%','interpreter','latex','Fontsize',18)
!mkdir -p output
exportgraphics(hf,['output/leak_' num2str(leak) '_' circuit '_A_' num2str(A) '_Ara_' num2str(ara) '.pdf'],'Resolution',300)


% rotate, 2D plot
% color limit: between zero and max after t_init min
t_init = 1e3;
pos_init = find(t>t_init,1,'first');
switch circuit
    case 'N1'
        hf = figure;
        hc = create_2D_plot(t,c(:,1),1e3,"#77AC30",[0,max(c(pos_init:end,1))]);
        ylabel(hc,'\bf\boldmath$c$ / molecule','FontSize',18,'interpreter','latex')
        exportgraphics(hf,['output/leak_' num2str(leak) '_' circuit '_A_' num2str(A) '_Ara_' num2str(ara)  '_2D.pdf'],'Resolution',300)
    case {'open','closed'}
        % set(hf,'Position',[61         522        1762         335])
        % tiledlayout(1,4,"TileSpacing","tight")
        % nexttile
        hf = figure;
        hc = create_2D_plot(t,c(:,1),1e3,"#A2142F",[0,max(c(pos_init:end,1))]);
        ylabel(hc,'\bf\boldmath$c$ / molecule','FontSize',18,'interpreter','latex')
        exportgraphics(hf,['output/leak_' num2str(leak) '_' circuit '_A_' num2str(A) '_Ara_' num2str(ara)  '_red_2D.pdf'],'Resolution',300)
        % nexttile
        hf = figure;
        hc = create_2D_plot(t,c(:,2),1e3,"#0072BD",[0,max(c(pos_init:end,2))]);
        ylabel(hc,'\bf\boldmath$c$ / molecule','FontSize',18,'interpreter','latex')
        exportgraphics(hf,['output/leak_' num2str(leak) '_' circuit '_A_' num2str(A) '_Ara_' num2str(ara)  '_blue_2D.pdf'],'Resolution',300)
        % nexttile
        hf = figure;
        hc = create_2D_plot(t,c(:,3),1e3,"#77AC30",[0,max(c(pos_init:end,3))]);
        ylabel(hc,'\bf\boldmath$c$ / molecule','FontSize',18,'interpreter','latex')
        exportgraphics(hf,['output/leak_' num2str(leak) '_' circuit '_A_' num2str(A) '_Ara_' num2str(ara)  '_green_2D.pdf'],'Resolution',300)
end

switch circuit
    case 'closed'
        save(['output/leak_' num2str(leak) '_' circuit '_A_' num2str(A) '_Ara_' num2str(ara) '.mat'],'t','c')
end

%% plot the circuit
figure
Ecoli.data.original_graph = false;
Ecoli.data.layout_type = 'layered';
h = Ecoli.make_graph();
axis off

exportgraphics(gcf,'output/N3_Jung_make_graph.pdf','Resolution',300)