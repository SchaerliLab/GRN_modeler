%% phase space graph for the togglelator

clc
clear
close all

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

%% create protease oscillator model (exapmles/protease_oscillator.mat)

clean_up_GRN
app.Ecoli = Cell('Tomazou');
app.Ecoli = app.Ecoli.add_node('N1','type1');
app.Ecoli = app.Ecoli.add_regulator('Repression_in','N1','HILL','P_N1');
app.Ecoli = app.Ecoli.add_protease('N1','PROT1','type1');
app.Ecoli.set('PROT1','InitialAmount',5.000000e+01,'PROT1');
set(getconfigset(app.Ecoli.data.Mobj),'Stoptime',1e2);
set(get(getconfigset(app.Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-10);
set(get(getconfigset(app.Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-8);
app.Ecoli.data.StatesToLog = [app.Ecoli.data.StatesToLog, {'uP_N1'}];
app.Ecoli.data.StatesToLog = [app.Ecoli.data.StatesToLog, {'P_N1'}];
app.Ecoli.data.StatesToLog = [app.Ecoli.data.StatesToLog, {'mRNA_N1'}];
app.Ecoli.data.StatesToLog(strcmp(app.Ecoli.data.StatesToLog,'mRNA_N1')) = [];
app.Ecoli.data.StatesToLog(strcmp(app.Ecoli.data.StatesToLog,'uP_N1')) = [];
app.Ecoli.data.StatesToLog(strcmp(app.Ecoli.data.StatesToLog,'P_N1')) = [];
app.Ecoli.data.StatesToLog = [app.Ecoli.data.StatesToLog, {'mRNA_N1'}];
app.Ecoli.data.StatesToLog = [app.Ecoli.data.StatesToLog, {'uP_N1'}];
app.Ecoli.data.StatesToLog = [app.Ecoli.data.StatesToLog, {'P_N1'}];

% get the model
Mobj_prot = app.Ecoli.get_model();
save 'protease_oscillator_model' Mobj_prot
clear Mobj_prot

%% create the model (examples/togglelator.mat)

clean_up_GRN
app.Ecoli = Cell();
app.Ecoli = app.Ecoli.add_node('N1','type1');
app.Ecoli = app.Ecoli.add_node('N2','type1');
app.Ecoli = app.Ecoli.add_regulator('Repression_in','N2','HILL','P_N1');
app.Ecoli = app.Ecoli.add_regulator('Repression_in','N1','HILL','P_N2');
app.Ecoli = app.Ecoli.add_protease('N1','PROT1','type1');
app.Ecoli = app.Ecoli.add_protease('N2','PROT1','type1');
app.Ecoli.set('P_N2','InitialAmount',1.000000e+02,'N2');
app.Ecoli.set('uP_N2','InitialAmount',1.000000e+01,'N2');
app.Ecoli.data.StatesToLog = [app.Ecoli.data.StatesToLog, {'P_N1'}];
app.Ecoli.data.StatesToLog = [app.Ecoli.data.StatesToLog, {'P_N2'}];
set(getconfigset(app.Ecoli.data.Mobj),'Stoptime',1e3);
set(getconfigset(app.Ecoli.data.Mobj),'Stoptime',1e2);

%% simulation settings

% get the model
Mobj = app.Ecoli.get_model();

% get the protease oscillator model
load('protease_oscillator_model.mat','Mobj_prot')
% change the parameters to be in accordance with the togglelator
set(sbioselect(Mobj_prot.Parameters,'Name','K_protease'),'Value',get(sbioselect(Mobj_prot.Parameters,'Name','K_protease'),'Value')/2)

% end of the simulation
t_end = 1e4;

% simulation settings:
configsetobj = getconfigset(Mobj);
set(get(configsetobj,'RuntimeOptions'),'StatesToLog',{'P_N1','P_N2','uP_N1'}) 
set(configsetobj,'SolverType','sundials','Stoptime',t_end)
set(get(configsetobj,'SolverOptions'),'AbsoluteTolerance',1e-10,...
    'RelativeTolerance',1e-8)
%the same for the protease oscillator
configsetobj = getconfigset(Mobj_prot);
set(get(configsetobj,'RuntimeOptions'),'StatesToLog',{'mRNA_N1','uP_N1','P_N1'}) 
set(configsetobj,'SolverType','sundials','Stoptime',t_end)
set(get(configsetobj,'SolverOptions'),'AbsoluteTolerance',1e-10,...
    'RelativeTolerance',1e-8)

disp('accelerate')
tic
sbioaccelerate(Mobj);
sbioaccelerate(Mobj_prot);
toc

%% fixed points

% protease concentration
c_prot_act = 100;

% fixed points
syms c1 c2 c3 c4 c5 c6 c7
assume([c1,c2,c3,c4,c5,c6,c7] >= 0)

% get ode from simbiology model
[f_model,f] = simbio2ode(Mobj);
% convert to symboic function
f = sym(f);
% model parameters
p = parameters2struct(Mobj);

% set the protease concentration
f_prot = subs(f,sym('c7'),c_prot_act);
% find the fixed points
sol = solve(f_prot);

% fixed points
P1 = log10(double(vpa(sol.c3)));
P2 = log10(double(vpa(sol.c6)));
uP1 = log10(double(vpa(sol.c2)));

%% plot

% we need every species to set the initial concentrations later
set(get(getconfigset(Mobj_prot),'RuntimeOptions'),'StatesToLog',{Mobj_prot.Species.Name})

% find the periodic orbit with the protease oscillator
Mobj_prot.Species(4).Value = c_prot_act/2;
[~,c_per,~] = sbiosimulate(Mobj_prot);
% set the initial concentration according to the periodic orbit
for j = 1:3
    Mobj.Species(j).Value = c_per(end,j);
    Mobj.Species(j+3).Value = c_per(end,j);
end
Mobj.Species(7).Value = c_prot_act;


% simulate
[~,c,names] = sbiosimulate(Mobj);


% plot the result
hf = figure;
set(gcf,'Position',[568   280   481   418])
hold on

c = log10(c);
plot3(c(:,strcmp(names,'P_N1')),c(:,strcmp(names,'P_N2')),c(:,strcmp(names,'uP_N1')),'LineWidth',4,'Color',"#EDB120")


linewidth = 1;
% spiraling out
for j = 1:3
    Mobj.Species(j).Value = c_per(end,j)/4;
    Mobj.Species(j+3).Value = c_per(end,j)*4;
end
% simulate
[~,c,names] = sbiosimulate(Mobj);
c = log10(c);
plot3(c(:,strcmp(names,'P_N1')),c(:,strcmp(names,'P_N2')),c(:,strcmp(names,'uP_N1')),'LineWidth',linewidth,'Color',"#A2142F")

% spiraling in
for j = 1:3
    Mobj.Species(j).Value = c_per(end,j)/3.5;
    Mobj.Species(j+3).Value = c_per(end,j)*3.5;
end
% simulate
[~,c,names] = sbiosimulate(Mobj);
c = log10(c);
plot3(c(:,strcmp(names,'P_N1')),c(:,strcmp(names,'P_N2')),c(:,strcmp(names,'uP_N1')),'LineWidth',linewidth,'Color',"#EDB120")

% spiraling out
for j = 1:3
    Mobj.Species(j).Value = c_per(end,j)*4;
    Mobj.Species(j+3).Value = c_per(end,j)/4;
end
% simulate
[~,c,names] = sbiosimulate(Mobj);
c = log10(c);
plot3(c(:,strcmp(names,'P_N1')),c(:,strcmp(names,'P_N2')),c(:,strcmp(names,'uP_N1')),'LineWidth',linewidth,'Color',"#0072BD")

% spiraling in
for j = 1:3
    Mobj.Species(j).Value = c_per(end,j)*3.5;
    Mobj.Species(j+3).Value = c_per(end,j)/3.5;
end
% simulate
[~,c,names] = sbiosimulate(Mobj);
c = log10(c);
plot3(c(:,strcmp(names,'P_N1')),c(:,strcmp(names,'P_N2')),c(:,strcmp(names,'uP_N1')),'LineWidth',linewidth,'Color',"#EDB120")

xlabel('\bf\boldmath$\log_{10}{\left(\left[P_{N1}\right] / \mathrm{molecule}\right)}$','interpreter','latex','Fontsize',14)
ylabel('\bf\boldmath$\log_{10}{\left(\left[P_{N2}\right] / \mathrm{molecule}\right)}$','interpreter','latex','Fontsize',14)
zlabel('\bf\boldmath$\log_{10}{\left(\left[uP_{N1}\right] / \mathrm{molecule}\right)}$','interpreter','latex','Fontsize',14)
set(gca,'LineWidth',2,'Fontsize',16)
grid on
view(-64,22)
axis equal

% plane 
ax = axis;
[X, Y, Z] = equalXYPlane(0,ax(1:2),ax(5:6)); % X = Y, Z ranges from 0 to 5
surf(X, Y, Z, 'FaceColor', [0.5, 0.5, 0.5], 'FaceAlpha', 0.4, 'EdgeColor', 'none'); 

% box all around
box on
ax = gca;
ax.BoxStyle = 'full';

% fixed points
scatter3(P1(2),P2(2),uP1(2),50,[0, 114, 189]/255,'filled') % stable1
scatter3(P1(3),P2(3),uP1(3),50,[162, 20, 47]/255,'filled') % stable2
scatter3(P1(1),P2(1),uP1(1),50,[126, 47, 142]/255,'filled') % unstalble

% Create arrow
annotation(gcf,'arrow',[0.393729231428877 0.327062564762211],...
    [0.4808497423629 0.433002852410747],'Color',"#A2142F",'LineWidth',4);
% Create arrow
annotation(gcf,'arrow',[0.403848279047925 0.463497401854942],...
    [0.488072344322345 0.531134545279283],'Color',"#EDB120",'LineWidth',4);
% Create arrow
annotation(gcf,'arrow',[0.603069537387256 0.560807632625352],...
    [0.66668956043956 0.595728065793856] ,'Color',"#EDB120",'LineWidth',4);
% Create arrow
annotation(gcf,'arrow',[0.610651629072682 0.647493734335839],...
    [0.679425837320574 0.755980861244019],'Color',"#0072BD",'LineWidth',4);

exportgraphics(gcf,'output/toggleswitch_phasespace.pdf','Resolution',300)

function [X, Y, Z] = equalXYPlane(value, range, rangeZ)
% EQUALXYPLANE Generates a 3D plane where X = Y.
% 
% INPUTS:
%   value   - The constant difference or offset for X = Y (e.g., value = 0 gives x = y).
%   range   - Range for X (or Y) as [min, max].
%   rangeZ  - Range for Z (vertical dimension) as [zmin, zmax].
%
% OUTPUTS:
%   X, Y, Z - Meshgrid coordinates for the plane where X = Y.
%
% Example:
%   [X, Y, Z] = equalXYPlane(0, [0, 10], [0, 5]);
%   surf(X, Y, Z);

    % Generate the line where X = Y
    diagRange = range(1):0.1:range(2);
    [X, Z] = meshgrid(diagRange, rangeZ(1):0.1:rangeZ(2));
    Y = X + value; % Enforce X = Y (with optional offset)
end
