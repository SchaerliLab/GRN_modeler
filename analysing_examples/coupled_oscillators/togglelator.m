%% togglelator
% the toggle switch can oscillate

clc
clear
close all

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))

%% create protease oscillator model (exapmles/protease_oscillator.mat)

clean_up_GRN
Ecoli = Cell('Tomazou');
Ecoli = Ecoli.add_node('N1','type1');
Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL','P_N1');
Ecoli = Ecoli.add_protease('N1','PROT1','type1');
Ecoli.set('PROT1','InitialAmount',5.000000e+01,'PROT1');
set(getconfigset(Ecoli.data.Mobj),'Stoptime',1e2);
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-10);
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-8);
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'uP_N1'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'mRNA_N1'}];
Ecoli.data.StatesToLog(strcmp(Ecoli.data.StatesToLog,'mRNA_N1')) = [];
Ecoli.data.StatesToLog(strcmp(Ecoli.data.StatesToLog,'uP_N1')) = [];
Ecoli.data.StatesToLog(strcmp(Ecoli.data.StatesToLog,'P_N1')) = [];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'mRNA_N1'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'uP_N1'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];

% get the model
Mobj_prot = Ecoli.get_model();
save 'protease_oscillator_model' Mobj_prot
clear Mobj_prot

%% create the model (examples/togglelator.mat)

clean_up_GRN
Ecoli = Cell();
Ecoli = Ecoli.add_node('N1','type1');
Ecoli = Ecoli.add_node('N2','type1');
Ecoli = Ecoli.add_regulator('Repression_in','N2','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL','P_N2');
Ecoli = Ecoli.add_protease('N1','PROT1','type1');
Ecoli = Ecoli.add_protease('N2','PROT1','type1');
Ecoli.set('P_N2','InitialAmount',1.000000e+02,'N2');
Ecoli.set('uP_N2','InitialAmount',1.000000e+01,'N2');
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N2'}];
set(getconfigset(Ecoli.data.Mobj),'Stoptime',1e3);
set(getconfigset(Ecoli.data.Mobj),'Stoptime',1e2);

%% find fix points

% get the model
Mobj = Ecoli.get_model();

% Unstable Supercritical Hopf Bifurcation

% prescribed protease concentration
c_prot = [logspace(0,3.3,5),logspace(3.3,3.7,20),logspace(3.7,4,5)]; % pitchfok
% c_prot = [logspace(0,3.3,5),logspace(3.3,3.69,5),linspace(10^3.69,10^3.7,10),logspace(3.7,4,5)]; % pitchfok
% c_prot = [logspace(0,3.69,10),linspace(10^3.69,10^3.7,10)]; % pitchfok
% c_prot = logspace(-12,-6,10); % nothing there
n_prot = numel(c_prot);

syms c1 c2 c3 c4 c5 c6 c7
assume(c1 >= 0)
assume(c2 >= 0)
assume(c3 >= 0)
assume(c4 >= 0)
assume(c5 >= 0)
assume(c6 >= 0)
assume(c7 >= 0)
% assume([c1,c2,c3,c4,c5,c6,c7],'real')

% get ode from simbiology model
[f_model,f] = simbio2ode(Mobj);
% convert to symboic function
f = sym(f);
% model parameters
p = parameters2struct(Mobj);
% calculate the jacobian
% the protease (the last one is constant)
f_jac = calc_jacobian_from_ode(f_model,p);

% 3 fix points c_fix: nvar x nprot x nfix
c_fix = zeros(7,n_prot,3);
% stability of the fix point nfix x nprot
stable = zeros(3,n_prot);

tic
for np = 1:n_prot
    % set the protease concentration
    f_prot = subs(f,sym('c7'),c_prot(np));
    % find the fixed points
    sol = solve(f_prot);

    if numel(sol.c1)>3
        warning('We have more than 3 fixed points!')
    end

    for n_fix = 1:numel(sol.c1)
        % fix protease concentration
        c_fix(7,np,n_fix) = c_prot(np);
        for i = 1:6
            c_fix(i,np,n_fix) = vpa(sol.(['c' int2str(i)])(n_fix));
        end
        % evaluate the jacobian
        J = f_jac(c_fix(:,np,n_fix));
        [~,D] = eig(J);
        % check stability
        stable(n_fix,np) = all(real(diag(D(1:6,1:6)))<0);
    end

    % if we have less than 3 fix points
    nfix = numel(sol.c1);
    if nfix < 3
        stable(nfix+1:end,np) = nan(3-nfix,1);
        c_fix(1:6,np,numel(sol.c1)+1:end) = nan(6,1,3-nfix);
    end
end
toc

%% sort the fix points

% the fix points for P1
c_fix_P1 = reshape(c_fix(3,:,:),n_prot,3).';

% stable points first
[stable,I] = sort(stable,1,'descend','MissingPlacement','last');
for i = 1:n_prot
    c_fix_P1(:,i) = c_fix_P1(I(:,i),i);
end

% we start with the smallest stable fix point
for i = 1:n_prot
    [c_fix_P1(1:sum(stable(:,i)==1),i),I] = sort(c_fix_P1(1:sum(stable(:,i)==1),i),1,'ascend');
end

%% find min and max values for oscillation

% get more points here
n_prot_dense = 10*n_prot;
c_prot_dense = logspace(log10(c_prot(1)),log10(c_prot(end)),n_prot_dense);

% get the protease oscillator model
load('protease_oscillator_model.mat','Mobj_prot')
% change the parameters to be in accordance with the togglelator
set(sbioselect(Mobj_prot.Parameters,'Name','K_protease'),'Value',get(sbioselect(Mobj_prot.Parameters,'Name','K_protease'),'Value')/2)

% start and end of the simulation
t_min = 1e3;
t_end = t_min+1e3;

% simulation settings:
% configsetobj = getconfigset(Mobj);
% set(get(configsetobj,'RuntimeOptions'),'StatesToLog',{'P_N1'}) 
% set(configsetobj,'SolverType','sundials','Stoptime',t_end)
% set(get(configsetobj,'SolverOptions'),'AbsoluteTolerance',1e-14,...
%     'RelativeTolerance',1e-12)
% the same for the protease oscillator
configsetobj = getconfigset(Mobj_prot);
set(get(configsetobj,'RuntimeOptions'),'StatesToLog',{'mRNA_N1','uP_N1','P_N1'}) % to be able to continue the simulations
set(configsetobj,'SolverType','sundials','Stoptime',t_end)
set(get(configsetobj,'SolverOptions'),'AbsoluteTolerance',1e-10,...
    'RelativeTolerance',1e-8)

disp('accelerate')
tic
sbioaccelerate(Mobj);
sbioaccelerate(Mobj_prot);
toc

% maximum and minimum values in the oscillator
c_osc = zeros(2,n_prot_dense);

tic
% try to find periodic attractors, calculate the amplitude
for i = n_prot_dense:-1:1

%     if sum(stable(:,i)==0) ~= 1
%         warning('We have multiple unstable fix points.')
%     end

    % start from the periodic orbit from the corresponding protease
    % oscillator
    % we have half amount protein, we are using half amount protease
    Mobj_prot.Species(4).Value = c_prot_dense(i)/2;
    [t,c,names] = sbiosimulate(Mobj_prot);

    % % set the concentrations according to the goodwin oscillator
    % set(sbioselect(Mobj.Species,'Name','mRNA_N1'),'Value',c(end,strcmp(names,'mRNA_N1')))
    % set(sbioselect(Mobj.Species,'Name','mRNA_N2'),'Value',c(end,strcmp(names,'mRNA_N1')))
    % set(sbioselect(Mobj.Species,'Name','uP_N1'),'Value',c(end,strcmp(names,'uP_N1')))
    % set(sbioselect(Mobj.Species,'Name','uP_N2'),'Value',c(end,strcmp(names,'uP_N1')))
    % set(sbioselect(Mobj.Species,'Name','P_N1'),'Value',c(end,strcmp(names,'P_N1')))
    % set(sbioselect(Mobj.Species,'Name','P_N2'),'Value',c(end,strcmp(names,'P_N1')))
    % set(sbioselect(Mobj.Species,'Name','PROT1'),'Value',c_prot_dense(i))
    % [t,c,names] = sbiosimulate(Mobj);

    % get rid of the beginning
    start_pos = find(t>=t_min,1,'first');

    % if it is constant, we might not have enogh data points
    if size(c,1)-start_pos > 3
        pks = findpeaks(c(start_pos:end,strcmp(names,'P_N1')),t(start_pos:end),'MinPeakProminence',1);
        % if it is not oscillating, we do not have peaks
        if ~isempty(pks)
            % maximum paks
            c_osc(1,i) = max(pks);
            % minimum peaks
            pks = findpeaks(-c(start_pos:end,strcmp(names,'P_N1')),t(start_pos:end),'MinPeakProminence',1);
            c_osc(2,i) = min(-pks);
        end
    end

end
toc

%% plot the bifurcation diagram

% load the unstable orbits calculated by continuation_togglelator
load('output/c0_UPO','extremum_UPO','c0_UPO')

% beginning of the monostable state
firtsnan = find(isnan(stable(2,:)),1,'first');
if isempty(firtsnan)
    firtsnan = n_prot;
end
% monostable solution, copy the first element to every other solution
c_fix_P1(2,firtsnan) = c_fix_P1(1,firtsnan);
c_fix_P1(3,firtsnan) = c_fix_P1(1,firtsnan);

% % the state where we have 3 fix point + the first point from the monostable
% c_prot_3fix = [c_prot(~isnan(c_fix_P1(2,:))),

% periodic orbit stability region
cstab_min_pos = find(c_prot_dense>c0_UPO(1,end),1,'first');
cstab_max_pos = find(c_prot_dense<c0_UPO(end,end),1,'last');
%%
hf = figure;
hold on
% fixed points
plot(c_prot(1:firtsnan),c_fix_P1(2,1:firtsnan),'Color',"#0072BD",'LineWidth',3) % stable1
plot(c_prot(1:firtsnan),c_fix_P1(1,1:firtsnan),'Color',"#A2142F",'LineWidth',3) % stable2
plot(c_prot(1:firtsnan),c_fix_P1(3,1:firtsnan),'--','Color','#7E2F8E','LineWidth',3) % unstalble
% periodic orbits
plot(c_prot_dense(1,1:cstab_min_pos),c_osc(1,1:cstab_min_pos),'--','Color','#EDB120','LineWidth',1) % unstable oscillation max
plot(c_prot_dense(1,1:cstab_min_pos),c_osc(2,1:cstab_min_pos),'--','Color','#EDB120','LineWidth',1) % unstable oscillation min
plot(c_prot_dense(1,cstab_min_pos:cstab_max_pos),c_osc(1,cstab_min_pos:cstab_max_pos),'Color','#EDB120','LineWidth',1) % stable oscillation max
plot(c_prot_dense(1,cstab_min_pos:cstab_max_pos),c_osc(2,cstab_min_pos:cstab_max_pos),'Color','#EDB120','LineWidth',1) % stable oscillation min
plot(c_prot_dense(1,cstab_max_pos:end),c_osc(1,cstab_max_pos:end),'--','Color','#EDB120','LineWidth',1) % unstable oscillation max
plot(c_prot_dense(1,cstab_max_pos:end),c_osc(2,cstab_max_pos:end),'--','Color','#EDB120','LineWidth',1) % unstable oscillation min
if firtsnan ~= n_prot % monostable
    plot(c_prot(firtsnan:end),c_fix_P1(1,firtsnan:end),'Color','#7E2F8E','LineWidth',3) % monostable
end
% unstable orbits with continuation method
plot(c0_UPO(:,end),extremum_UPO(:,1),'--','Color',"#D95319",'LineWidth',0.5)
plot(c0_UPO(:,end),extremum_UPO(:,2),'--','Color',"#D95319",'LineWidth',0.5)
plot(c0_UPO(:,end),extremum_UPO(:,3),'--','Color',"#D95319",'LineWidth',0.5)
plot(c0_UPO(:,end),extremum_UPO(:,4),'--','Color',"#D95319",'LineWidth',0.5)
set(gca,'YScale','log')
set(gca,'XScale','log')
xlabel('\bf\boldmath$\left[\textrm{PROT}_1\right]$ / molecule','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$\left[P_{N1}\right]$ / molecule','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)

% size of the dots in the scatter
scattersize = 100;
% find the bifurcation point in the diagram
subk_HOPF_start = find(c_osc(1,:)>0,1,'first');
subk_HOPF_end = find(c_osc(1,:)>0,1,'last');
% bifurcation points
% pitchfork
scatter(c_prot(1,firtsnan),c_fix_P1(3,firtsnan),scattersize,[0.4660 0.6740 0.1880],'filled')
% subcritical HOPF
scatter(c_prot_dense(1,[subk_HOPF_start,subk_HOPF_end]),...
    (c_osc(1,[subk_HOPF_start,subk_HOPF_end])+c_osc(2,[subk_HOPF_start,subk_HOPF_end]))/2,scattersize,[0.9290 0.6940 0.1250],'filled')
% subcritical pitchfork of POs
scatter([c_prot_dense(1,cstab_min_pos),c_prot_dense(1,cstab_min_pos),c_prot_dense(1,cstab_max_pos),c_prot_dense(1,cstab_max_pos)],...
    [c_osc(1,cstab_min_pos),c_osc(2,cstab_min_pos),c_osc(1,cstab_max_pos),c_osc(2,cstab_max_pos)],scattersize,[0.6350 0.0780 0.1840],'filled')

hl = legend('','','','','','','','','','','','','','','Pitchfork','Subcritical HOPF','Subcritical pitchfork of periodic orbits','Location','southwest','FontSize',14,'Box','off','interpreter','latex');
title(hl,'\bf Bifurcations \hspace{4 cm}','interpreter','latex')
exportgraphics(hf,'output/togglelator_bif.pdf','Resolution',300)

%% compare the togglelator and the protease oscillator
c_prot_act = 100;

% we need every species to set the initial concentrations later
set(get(getconfigset(Mobj_prot),'RuntimeOptions'),'StatesToLog',{Mobj_prot.Species.Name})

% find the periodic orbit with the protease oscillator
Mobj_prot.Species(4).Value = c_prot_act/2;
[t,c,names] = sbiosimulate(Mobj_prot);
% set the initial concentration according to the periodic orbit
for j = 1:3
    Mobj.Species(j).Value = c(end,j)*(1+0.1*randn);
    Mobj.Species(j+3).Value = c(end,j)*(1+0.1*randn);
    Mobj_prot.Species(j).Value = c(end,j);
end
Mobj.Species(7).Value = c_prot_act;

set(get(getconfigset(Mobj_prot),'RuntimeOptions'),'StatesToLog',{'P_N1'})
set(getconfigset(Mobj_prot),'Stoptime',1e2)
% set(getconfigset(Mobj),'Stoptime',1e3)
set(getconfigset(Mobj),'Stoptime',1e2)


% compare the simulations
[t1,c1,names1] = sbiosimulate(Mobj);
[t2,c2,names2] = sbiosimulate(Mobj_prot);

% plot the result
hf = figure;
hold on
plot(t1,c1(:,strcmp(names1,'P_N1')),'LineWidth',2)
plot(t2,c2(:,strcmp(names2,'P_N1')),'LineWidth',2)
xlabel('\bf\boldmath$t$ / min','interpreter','latex','Fontsize',18)
ylabel('\bf\boldmath$\left[P_{N1}\right]$ / molecule','interpreter','latex','Fontsize',18)
set(gca,'LineWidth',2,'Fontsize',16)
legend('protease oscillator','toggle switch')
exportgraphics(hf,'output/protosc_togglelator_comparision.pdf','Resolution',300)