%% goodwin (prolator)
% the toggle switch can oscillate

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

%% find fix points

% get the model
Mobj = app.Ecoli.get_model();

% prescribed protease concentration
c_prot = logspace(0,4,1000);
n_prot = numel(c_prot);

syms c1 c2 c3 c4
assume(c1 >= 0)
assume(c2 >= 0)
assume(c3 >= 0)
assume(c4 >= 0)
% assume([c1,c2,c3,c4],'real')

% get ode from simbiology model
[f_model,~] = simbio2ode(Mobj);
% convert to symboic function
% f = sym(f);
% model parameters
p = parameters2struct(Mobj);
% calculate the jacobian
% the protease (the last one is constant)
f_jac = calc_jacobian_from_ode(f_model,p);

% fix points c_fix: nvar x nprot
c_fix = zeros(4,n_prot);
% stability of the fix point nprot
stable = zeros(n_prot,1);

% optiond for searching the fixed points
options = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','off','OptimalityTolerance',1e-6,'FunctionTolerance',1e-6,'StepTolerance',1e-6,MaxFunctionEvaluations=1e5);

tic
for np = 1:n_prot
    % % set the protease concentration
    % f_prot = subs(f,sym('c4'),c_prot(np));
    % % find the fix point
    % sol = vpasolve(f_prot);
    % 
    % if numel(sol.c1)>1
    %     error('We have more than 1 fix points!')
    % else
    %     c_fix(1,np) = sol.c1;
    %     c_fix(2,np) = sol.c2;
    %     c_fix(3,np) = sol.c3;
    %     c_fix(4,np) = c_prot(np);
    % end

    % start from the previous fixed point
    % except if we are changing the parameters radically
    fix_point = [1,1,1];

    % search for the solution
    [fix_point,~,exitflag] = fsolve(@(x)f_model([10.^x(1),10.^x(2),10.^x(3),c_prot(np)]),log10(fix_point),options);
    c_fix(:,np) = [10.^fix_point,c_prot(np)];

    % evaluate the jacobian
    J = f_jac(c_fix(:,np));
    [~,D] = eig(J);
    % check stability
    stable(np) = all(real(diag(D(1:3,1:3)))<0);
end
toc

%% stability region

syms k_mat

% get ode from simbiology model
[f_model,f] = simbio2ode(Mobj,'k_mat');
% convert to symboic function
f = sym(f);
% model parameters
p = parameters2struct(Mobj);
% calculate the jacobian
% k_mat is an extra input
f_jac = calc_jacobian_from_ode(f_model,p,'k_mat');

% number of divisions:
n_prot = 200;
n_kmat = 200;

% protease
c_prots = logspace(0,4,n_prot);
k_mats = logspace(-6,2,n_kmat);

% type of the fix point: 0: complex unstable, 1: stable, 2: unstable real
fix_type = zeros(n_kmat,n_prot);

options = optimoptions('fsolve','Algorithm','levenberg-marquardt','Display','off','OptimalityTolerance',1e-6,'FunctionTolerance',1e-6,'StepTolerance',1e-6,MaxFunctionEvaluations=1e5); % ,'Display','off'

tic

for np = 1:n_prot
    for nk = 1:n_kmat
        % % % set the protease concentration and k_mat
        % f_act = subs(f,sym('c4'),c_prots(np));
        % f_act = subs(f_act,k_mat,k_mats(nk));
        % % find the fix point
        % sol = vpasolve(f_act);

        % if numel(sol.c1)>1
        %     error('We have more than 1 fix points!')
        % else
        %     fix_point = [sol.c1,sol.c2,sol.c3,c_prots(np)];
        % end

        % start from the previous fixed point
        % except if we are changing the parameters radically
        if np==1
            fix_point = [1,1,1];
        end
        
        % search for the solution
        [fix_point,~,exitflag] = fsolve(@(x)f_model([10.^x(1),10.^x(2),10.^x(3),c_prots(np)],k_mats(nk)),log10(fix_point),options);
        fix_point = 10.^(fix_point); % real

        % evaluate the jacobian
        J = f_jac([fix_point,c_prots(np)],k_mats(nk));
        % get rid of the constant PROT variable
        J = J(1:3,1:3);
        [~,D] = eig(J);
        % check stability
        if all(real(diag(D))<0)
            fix_type(nk,np) = 1; % stable
        elseif isreal(D)
            fix_type(nk,np) = 2; % unstable real
        else
            fix_type(nk,np) = 0; % unstable complex
        end

    end
end
toc

%% plot stability region

[C_PROTS,K_MATS] = meshgrid(c_prots,k_mats);

hf = figure;
surf(C_PROTS,K_MATS,fix_type)
shading flat
view(2)
set(gca,'Layer','top')
set(gca,'YScale','log')
set(gca,'XScale','log')
xlabel('\bf\boldmath$\left[\textrm{PROT}_1\right]$ / molecule','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$k_{mat}$ / min$^{-1}$','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
box on
grid off
set(gca,'Ylim',[min(k_mats),max(k_mats)])

% Create dummy patches for legend
hold on
h1 = patch(NaN, NaN, 'y');  % yellow
h2 = patch(NaN, NaN, 'b');  % blue
% Add legend
legend([h1, h2], {'Stable fixed point', 'Oscillation'})

exportgraphics(hf,'output/goodwin_stability.pdf','Resolution',300)


%% find min and max values for oscillation

% start and end of the simulation
t_end = 2e3;
t_min = 1e3;

n_prot = numel(stable);

% simulation settings:
configsetobj = getconfigset(Mobj);
set(get(configsetobj,'RuntimeOptions'),'StatesToLog',{'P_N1'})
set(configsetobj,'SolverType','sundials','Stoptime',t_end)
set(get(configsetobj,'SolverOptions'),'AbsoluteTolerance',1e-10,...
    'RelativeTolerance',1e-8)

disp('accelerate')
tic
sbioaccelerate(Mobj);
toc

% maximum and minimum values in the oscillator
c_osc = zeros(2,n_prot);

tic
% try to find periodic attractors, calculate the amplitude
for i = 1:n_prot

    % check just the unstable one
    if stable(i)==true
        continue;
    end

    % start from the periodic orbit from the corresponding protease
    % oscillator
    % we have half amount protein, we are using half amount protease
    Mobj.Species(4).Value = c_prot(i);
    [t,c,names] = sbiosimulate(Mobj);

    % get rid of the beginning
    start_pos = find(t>=t_min,1,'first');

    % if it is constant, we might not have enogh data points
    if size(c,1)-start_pos > 3
        pks = findpeaks(c(start_pos:end,strcmp(names,'P_N1')),t(start_pos:end),'MinPeakProminence',1);
        % if it is not oscillating, we do not have peaks
        if ~isempty(pks)
            % maximum peaks
            c_osc(1,i) = max(pks);
            % minimum peaks
            pks = findpeaks(-c(start_pos:end,strcmp(names,'P_N1')),t(start_pos:end),'MinPeakProminence',1);
            c_osc(2,i) = min(-pks);
        end
    end

end
toc

%% plot the bifurcation diagram

% end of the first stable region
stable1 = find(stable==false,1,'first')-1;
% beginning of the second stable region
stable2 = find(stable==false,1,'last')+1;

hf = figure;
hold on
plot(c_prot(1:stable1),c_fix(3,1:stable1),'Color','#0072BD','LineWidth',2)
plot(c_prot(stable2:end),c_fix(3,stable2:end),'Color','#0072BD','LineWidth',2)
plot(c_prot(stable==false),c_fix(3,stable==false),'--','Color','#0072BD','LineWidth',2) % unstalble
plot(c_prot(stable==false),c_osc(1,stable==false),'Color','#EDB120','LineWidth',1) % oscillation max
plot(c_prot(stable==false),c_osc(2,stable==false),'Color','#EDB120','LineWidth',1) % oscillation min
set(gca,'YScale','log')
set(gca,'XScale','log')
xlabel('\bf\boldmath$\left[\textrm{PROT}_1\right]$ / molecule','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$\left[P_{N1}\right]$ / molecule','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)

% size of the dots in the scatter
scattersize = 100;
% find the bifurcation point in the diagram
HOPF_start = find(stable==false,1,'first');
HOPF_end = find(stable==false,1,'last');
% bifurcation points
% HOPF
scatter(c_prot(1,[HOPF_start,HOPF_end]),c_fix(3,[HOPF_start,HOPF_end]),scattersize,[0.9290 0.6940 0.1250],'filled')

hl = legend('','','','','','Supercritical HOPF','Location','southwest','FontSize',14,'Box','off','interpreter','latex');
title(hl,'\bf Bifurcations \hspace{1 cm}','interpreter','latex')

exportgraphics(hf,'output/goodwin_bif.pdf','Resolution',300)