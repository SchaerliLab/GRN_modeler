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

Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N2'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N3'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N4'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N5'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N6'}];

Ecoli.set('dCas','InitialAmount',2.00000e+02,'N6'); % change dCas here
Ecoli.set('a1_N6','Value',1.1,'N6');
Ecoli.set('a1_N2','Value',1.1,'N2');

set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-10);
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-8);
set(getconfigset(Ecoli.data.Mobj),'Stoptime',1e5);
Ecoli.data.Accelerate = true;

solver = 'sundials';

Mobj = Ecoli.get_model();

% output time
configsetobj = getconfigset(Mobj);
set(configsetobj,'SolverType',solver);
configsetobj.SolverOptions.OutputTimes = linspace(0.5e5,1e5,1e5);

sbioaccelerate(Mobj);

%% And gate demonstration

% if dCas > 2e2 we save it with a different name
if get(sbioselect(Mobj,'Name','dCas'),'Value') > 2e2
    suffix = '_SI';
else
    suffix = '';
end

% a1 values for the N1 and N4 nodes
a1_values = [1,1;1.2,1;1,1.2;1.2,1.2];

tic
for i = 1:size(a1_values,1)

    set(sbioselect(Mobj.Parameters,'Name','a1_N1'),'Value',a1_values(i,1))
    set(sbioselect(Mobj.Parameters,'Name','a1_N4'),'Value',a1_values(i,2))

    [t,c,names] = run_simulation(Mobj,solver);

    % calculate the correlation matrix between the time series
    M = corr(c);

    % figure
    % hold on
    % for j = 1:numel(Ecoli.data.StatesToLog)
    %     plot(t,c(:,strcmp(names,Ecoli.data.StatesToLog{j})))
    % end
    % xlabel('time')
    % ylabel('concentration')
    % legend(Ecoli.data.StatesToLog,'Interpreter','none')
    
    t_endpos = find(t>=t(end)*0.997,1,'first');
    figure
    set(gcf,'Position',[581   543   570   255])
    hold on
    % linewidth according to correlation
    plot(t(t_endpos:end),c(t_endpos:end,strcmp(names,'P_N1')),'LineWidth',1+heaviside(M(1,1))*M(1,1)*3)
    plot(t(t_endpos:end),c(t_endpos:end,strcmp(names,'P_N4')),'LineWidth',1+heaviside(M(4,1))*M(4,1)*3)
    plot(t(t_endpos:end),c(t_endpos:end,strcmp(names,'P_N5')),'LineWidth',1+heaviside(M(5,1))*M(5,1)*3)
    plot(t(t_endpos:end),c(t_endpos:end,strcmp(names,'P_N6')),'LineWidth',1+heaviside(M(6,1))*M(6,1)*3)
    legend({'P_N1','P_N4','P_N5','P_N6'},'Interpreter','none')
    xlabel('\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
    ylabel('\bf\boldmath$c$ / molecule','interpreter','latex','Fontsize',18)
    set(gca,'LineWidth',2,'Fontsize',16)
    % draw the closest with wider line
    % child = get(gca,'children');
    % line_number = find(M(:,1)==max(M(2:end,1)));
    % switch line_number
    %     case 4
    %         set(child(end-1),'LineWidth',4)
    %     case 5
    %         set(child(end-2),'LineWidth',4)
    %     case 6
    %         set(child(end-3),'LineWidth',4)
    % end
    drawnow
    exportgraphics(gcf,['output/computer_' int2str(i) suffix '.pdf'],'Resolution',300)

    % figure
    % imagesc(M)
    % colorbar
    % axis square
    % clim([-1 1])

    figure
    variableNames = {'N_1', 'N_2', 'N_3', 'N_4', 'N_5', 'N_6'};
    h = heatmap(variableNames, variableNames, M, 'Colormap', winter);
    h.ColorLimits = [-1 1]; % Set color limits for consistent scaling
    set(gca,'FontSize',16)
    drawnow
    exportgraphics(gcf,['output/computer_corr_' int2str(i) suffix '.pdf'],'Resolution',300)

    % calculate the correlation matrix between the time series
    M = corr(c(:,[1,4,5,6]));

    figure
    variableNames = {'N_1', 'N_4', 'N_5', 'N_6'};
    h = heatmap(variableNames, variableNames, M, 'Colormap', winter);
    h.ColorLimits = [-1 1]; % Set color limits for consistent scaling
    % Get the heatmap's underlying axes for customization
    ax = gca;
    % Add a wide red rectangle around the first column
    % hold(ax, 'on'); % Keep the heatmap displayed while adding the rectangle
    xPos = ax.Position(1); % X-coordinate of the heatmap position
    yPos = ax.Position(2); % Y-coordinate of the heatmap position
    width = ax.Position(3) / size(M, 2); % Width of one column
    height = ax.Position(4); % Full height of the heatmap
    % Create the rectangle
    annotation('rectangle', ...
        [xPos, yPos, width, height], ... % Normalized coordinates
        'EdgeColor', 'red', ...
        'LineWidth', 10);
    set(gca,'FontSize',16)
    drawnow
    exportgraphics(gcf,['output/computer_corr_' int2str(i) suffix '_simple.pdf'],'Resolution',300)


    % microfluidics representation (phase space)
    figure
    plot(c(:,strcmp(names,'P_N1')),c(:,strcmp(names,'P_N4')),'Marker','.')
    xlabel('\bf\boldmath$\left[P_1\right]$ / molecule','interpreter','latex','Fontsize',18)
    ylabel('\bf\boldmath$\left[P_4\right]$ / molecule','interpreter','latex','Fontsize',18)
    set(gca,'FontSize',16,'LineWidth',2)
    exportgraphics(gcf,['output/computer_phasespace_' int2str(i) suffix '.pdf'],'Resolution',300)

end
toc

%% create nonlinear output

% output time
configsetobj = getconfigset(Mobj);
set(configsetobj,'SolverType',solver);
configsetobj.SolverOptions.OutputTimes = linspace(0.5e4,1e4,1e4);

sbioaccelerate(Mobj);

% number of divisions
n = 100;

% a1 values for the N1 and N4 nodes
set(sbioselect(Mobj.Parameters,'Name','a1_N1'),'Value',1.2)
a1_values = linspace(1.0,1.3,n);

% store the output vales
output = zeros(n,1);

tic
for i = 1:n

    set(sbioselect(Mobj.Parameters,'Name','a1_N4'),'Value',a1_values(i))

    [t,c,names] = run_simulation(Mobj,solver);


    % calculate the correlation matrix between the time series
    M = corr(c(:,[1,4]));

    output(i) = M(1,2);
end
toc

figure
plot(a1_values,output,'LineWidth',2)
xlabel('\bf\boldmath$a_{1,N4}$ / min$^{-1}$','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath $\mathrm{corr}(c_{N1}, c_{N4})$ / -','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
exportgraphics(gcf,'output/computer_nonlinear.pdf','Resolution',300)