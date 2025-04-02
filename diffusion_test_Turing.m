%% Diffusion test
clc
close all
clear

% add to path
path = strsplit(pwd,'GRN');
addpath(genpath([path{1},'GRN']))


%% build reaction kinetics model

clean_up_GRN

Ecoli = Cell('Elowitz');
solver = 'sundials';
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'AbsoluteTolerance',1e-8)
set(get(getconfigset(Ecoli.data.Mobj),'SolverOptions'),'RelativeTolerance',1e-6)
set(getconfigset(Ecoli.data.Mobj),'MaximumWallClock',10)
Ecoli = Ecoli.add_node('N1','type1');
Ecoli = Ecoli.add_node('N2','type1');
Ecoli = Ecoli.add_regulator('Activation_in','N1','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Activation_in','N2','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL','P_N2');
Ecoli = Ecoli.add_regulator('Repression_in','N2','HILL','P_N2');
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N2'}];
% higher nonlinearity
set(sbioselect(Ecoli.data.Mobj,'Name','n'),'Value',6);
set(sbioselect(Ecoli.data.Mobj,'Name','K_N2|-N2'),'Value',55);

% get the model
Mobj = Ecoli.get_model();


% set simulation time
set(getconfigset(Mobj),'StopTime',1e3);
% we have to follow every species
set(get(getconfigset(Mobj),'RuntimeOptions'),'StatesToLog',{Mobj.Species.Name})
[t,c,names] = sbiosimulate(Mobj);
% set output concentration for Mobj
for i = 1:numel(Mobj.Species)
    Mobj.Species(i).InitialAmount = c(end,i);
end

%% Parameters

% number of divisions in "y" direction -> n+1 grid points
nx1 = 100;
% number of divisions in "x" direction -> n+1 grid points
nx2 = 100;
nt = 1e4; % time steps
L1 = 1e2; % length in "y" direction
T = 1e3;

n_species = numel(Mobj.Species);    

D = zeros(n_species,1);
D(2) = 1;
D(4) = 10;

h = L1/nx1;
dt = T/nt;
L2 = h*nx2;

% number of plots during the simulation
n_plot = 1e2;


noise = 0.01;

%% Boundary conditions

% upper bounary: 'noflux', 'pbc' or 'dirichlet'
b.Boundary_up = 'pbc';
% down bounary: 'noflux', 'pbc' or 'dirichlet'
b.Boundary_down = 'pbc';
% left bounary: 'noflux', 'pbc' or 'dirichlet'
b.Boundary_left = 'pbc';
% right bounary: 'noflux', 'pbc' or 'dirichlet'
b.Boundary_right = 'pbc';

%% initial condition

c  = cell(n_species,1);
for i = 1:n_species
    c{i} = Mobj.Species(i).Value*(ones((nx1+1)*(nx2+1),1)+noise*rand((nx1+1)*(nx2+1),1));
end

%% building the matrices

% Laplace: L = Axx + Ayy;
[Axx,Ayy] = laplace_matrices(nx1,nx2,b,h);

%% Laplace ADI
% ADI method for the diffusion of the amino acids

% store the matrices separately for every diffusion coefficient
% this is not memory efficient, but we need LU decomposition only once
% (fast)

% unique diffusion coefficients and the position for the matrices for the
% individual species:
[D_unique,~,D_pos] = unique(D);
% number of the unique diffusion coefficients
n_diff = length(D_unique);

% M: Lxx, Uxx, Lyy, Uyy, Aexx, Aeyy
M = cell(n_diff,6);

% unit matrix
E = speye((nx1+1)*(nx2+1));

for i = 1:n_diff

    % implicit matrix:
    Aixx = E-D_unique(i)*dt/h^2/2*Axx;
    Aiyy = E-D_unique(i)*dt/h^2/2*Ayy;
    % LU decomposition:
    [Lxx,Uxx] = lu(Aixx);
    [Lyy,Uyy] = lu(Aiyy);

    % explicit matrix:
    Aexx = E+D_unique(i)*dt/h^2/2*Axx;
    Aeyy = E+D_unique(i)*dt/h^2/2*Ayy;

    % store the matrices
    M{i,1} = Lxx;
    M{i,2} = Uxx;
    M{i,3} = Lyy;
    M{i,4} = Uyy;
    M{i,5} = Aexx;
    M{i,6} = Aeyy;

end

clear Axx Ayy E Aixx Aiyy

%% Reaction kinetics settings

% % we have to follow every species
% set(get(getconfigset(Mobj),'RuntimeOptions'),'StatesToLog',{Mobj.Species.Name})
% 
% % we need output time at the end (for COPASI also in zero)
% set(get(getconfigset(Mobj),'SolverOptions'),'OutputTimes',[0,dt]);

% get the constant parameters from the model
p = get_parameters(Mobj,'k_ext');

% generate ODE function from the model
simbio2ode(Mobj,'generated_model','cell');

%% Reaction zone for dirichlet boundary conditions
% if we do not have just zeros on the boundary, we should not calculate the
% reaction kinetics terms  there and we need this part

% % size of the reaction zone
% R = true(nx1+1,nx2+1);
% 
% % Dirichlet boundary conditions
% switch Boundary_up
%     case 'dirichlet'
%         R(1,:) = false;
% end
% switch Boundary_down
%     case 'dirichlet'
%         R(nx1+1,:) = false;
% end
% switch Boundary_left
%     case 'dirichlet'
%         R(:,1) = false;
% end
% switch Boundary_right
%     case 'dirichlet'
%         R(:,nx2+1) = false;
% end
% 
% % get a column vector
% R = R(:);

%% Simulation

% serial number of the followed species
followed = find(contains({Mobj.Species.Name},Ecoli.data.StatesToLog));

% size of the tiledlayout
nplot1 = floor(sqrt(length(followed)));
nplot2 = ceil(sqrt(length(followed)));
while nplot1*nplot2 < length(followed)
    nplot2 = nplot2+1;
end

% store the handles for the plot
h_plot = cell(length(followed),1);

hf = figure;
set(hf,'Position',[104         189        1710         872])
tiledlayout(nplot1,nplot2,'TileSpacing','tight');
for nplot = 1:length(followed)
    nexttile
    h_plot{nplot} = nextsurf(reshape(c{followed(nplot)},[nx1+1,nx2+1]),L1,L2,nx1,nx2,Ecoli.data.StatesToLog{nplot});
end
drawnow

tic
for i = 1:nt

    % === Diffusion ===
    for j = 1:n_species
        % if the diffusion coefficient is not zero
        if D(j) ~= 0
            % impl: L(U*c) = c
            % x:expl, y:impl
            c{j} = M{D_pos(j),3}\(M{D_pos(j),5}*c{j});
            c{j} = M{D_pos(j),4}\c{j};
            % x:impl, y:expl
            c{j} = M{D_pos(j),1}\(M{D_pos(j),6}*c{j});
            c{j} = M{D_pos(j),2}\c{j};
        end
    end

   % === Reactions ===
   dc = generated_model(i*dt, c, p);
   for j = 1:n_species
       c{j} = c{j} + dc{j}*dt;
   end

   % plot the data
   if rem(i,nt/n_plot) == 0

       % update the surf
       for nplot = 1:length(followed)
           set(h_plot{nplot},'ZData',reshape(c{followed(nplot)},[nx1+1,nx2+1]));
       end
       drawnow

       % hatralevo ido:
       disp([int2str(i/nt*100) '%, a varhato befejezes: ' ...
           num2str(toc/i*(nt-i),'%4.2f')])
    end

end
toc
