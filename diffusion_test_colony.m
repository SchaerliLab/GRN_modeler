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
Ecoli = Ecoli.add_node('N1','type1');
Ecoli = Ecoli.add_node('N2','type1');
Ecoli = Ecoli.add_node('N3','type1');
Ecoli = Ecoli.add_regulator('Repression_in','N2','HILL','P_N1');
Ecoli = Ecoli.add_regulator('Repression_in','N3','HILL','P_N2');
Ecoli = Ecoli.add_regulator('Repression_in','N1','HILL','P_N3');
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N1'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N2'}];
Ecoli.data.StatesToLog = [Ecoli.data.StatesToLog, {'P_N3'}];
Ecoli.set('P_N1','InitialAmount',9.800000e+01,'N1'); % Elowitz
Ecoli = Ecoli.add_regulator('Activation_out','N1','HILL','Light');
set(sbioselect(Ecoli.data.Mobj,'Name','Light'),'Constant',false);
addrule(Ecoli.data.Mobj,'Light = A*(1+sin(2*pi/T*time))/2','repeatedAssignment','Name','SinWave');
addparameter(Ecoli.data.Mobj,'A','Value',0.000000e+00,'Units','','Notes','Individual');
addparameter(Ecoli.data.Mobj,'T','Value',1.280000e+02,'Units','','Notes','Individual');
set(get(getconfigset(Ecoli.data.Mobj),'CompileOptions'),'Dimensionalanalysis',false)
set(sbioselect(Ecoli.data.Mobj,'Name','A'),'Value',1.000000e+02);
set(sbioselect(Ecoli.data.Mobj,'Name','T'),'Value',2.560000e+02);
set(getconfigset(Ecoli.data.Mobj),'Stoptime',1e3);
Ecoli.set('P_N1','InitialAmount',1.000000e+02,'N1'); % Elowitz

% get the model
Mobj = Ecoli.get_model();

% delete(Mobj.Reactions)

% add cell growth
addreaction(Mobj, 'null -> C','ReactionRate','k_cell*C*(C_unit-C)');
addparameter(Mobj,'C_unit',1,'Units','molecule'); % just to set up the units
addparameter(Mobj,'k_cell',.1,'Units','1/(molecule*minute)');
set(sbioselect(Mobj.Species,'Name','C'),'Units','molecule')

%% Parameters

% number of divisions in "y" direction -> n+1 grid points
nx1 = 100;
% number of divisions in "x" direction -> n+1 grid points
nx2 = 100;
nt = 1e3; % time steps
L1 = 1e2; % length in "y" direction
T = 1e3;

n_species = numel(Mobj.Species);    

D = zeros(n_species,1);
% cell diffusion coefficient
D(end) = .01;
% reaction zone for cells
cmin = 0*1e-3;
cmax = .99;

h = L1/nx1;
dt = T/nt;
L2 = h*nx2;

% number of plots during the simulation
n_plot = 1e2;


noise = 0.01*0;

%% Boundary conditions

% upper bounary: 'noflux', 'pbc' or 'dirichlet'
Boundary_up = 'pbc';
% down bounary: 'noflux', 'pbc' or 'dirichlet'
Boundary_down = 'pbc';
% left bounary: 'noflux', 'pbc' or 'dirichlet'
Boundary_left = 'pbc';
% right bounary: 'noflux', 'pbc' or 'dirichlet'
Boundary_right = 'pbc';

%% initial condition

sigma = 3;
height = 1;

% center position
center = round([nx1/2+1,nx2/2+1]);
% create gaussian cell distribution
G = gaussian2D(nx1+1, nx2+1, center, sigma, height);
% G = G-cmin;
% G(G<0) = 0;
% % keep the part which is above cmin
% G(G<cmin) = 0;

% center position with indexes
center = sub2ind([nx1+1,nx2+1],center(1),center(2));
c  = cell(n_species,1);
c{end} = reshape(G,[numel(G),1]);
for i = 1:n_species-1
    % c{i} = (c{end}>cmin).*(Mobj.Species(i).Value*(ones((nx1+1)*(nx2+1),1)+noise*rand((nx1+1)*(nx2+1),1)));
    c{i} = (Mobj.Species(i).Value*(ones((nx1+1)*(nx2+1),1)+noise*rand((nx1+1)*(nx2+1),1)));
end
clear G

%% Nabla: N = Ax + Ay;
% y direction:
Ay = spdiags([-ones(nx1+1,1) ones(nx1+1,1)],[-1 1],nx1+1,nx1+1);
switch Boundary_up
    case {'noflux','dirichlet'}
        Ay(1,:) = 0; % no flux
    case 'pbc'
        Ay(1,end) = -1; % pbc
end
switch Boundary_down
    case {'noflux','dirichlet'}
        Ay(end,:) = 0; % no flux
    case 'pbc'
        Ay(end,1) = 1; % pbc
end
Ay = kron(speye(nx2+1),Ay);
switch Boundary_left
        case 'dirichlet'
        Ay(1:nx1+1,:) = 0;
end
switch Boundary_right
        case 'dirichlet'
        Ay(end-(nx1+1)+1:end,:) = 0;
end

% x direction:
Ax = spdiags([-ones((nx1+1)*(nx2+1),1)  ones((nx1+1)*(nx2+1),1)],...
    [-(nx1+1) (nx1+1)],(nx1+1)*(nx2+1),(nx1+1)*(nx2+1));
switch Boundary_left
    case {'noflux','dirichlet'}
        Ax(1:nx1+1,:) = 0; % no flux
    case 'pbc'
        Ax(1:nx1+1,end-(nx1+1)+1:end) = -speye(nx1+1,nx1+1); % pbc
end
switch Boundary_right
    case {'noflux','dirichlet'}
        Ax(end-(nx1+1)+1:end,:) = 0; % no flux
    case 'pbc'
        Ax(end-(nx1+1)+1:end,1:nx1+1) = speye(nx1+1,nx1+1); % pbc
end
switch Boundary_up
    case 'dirichlet'
        Ax(1:nx1+1:end,:) = 0;
end
switch Boundary_down
    case 'dirichlet'
        Ax(nx1+1:nx1+1:end,:) = 0;
end


% with our parameters
Ax = Ax/(2*h);
Ay = Ay/(2*h);

%% Laplace: L = Axx + Ayy;
% y direction:
Ayy = spdiags([ones(nx1+1,1) -2*ones(nx1+1,1) ones(nx1+1,1)],[-1 0 1],nx1+1,nx1+1);
switch Boundary_up
    case 'noflux'
        Ayy(1,2) = 2; % no flux
    case 'pbc'
        Ayy(1,end) = 1; % pbc
    case 'dirichlet'
        Ayy(1,:) = 0;
end
switch Boundary_down
    case 'noflux'
        Ayy(end,end-1) = 2; % no flux
    case 'pbc'
        Ayy(end,1) = 1; % pbc
    case 'dirichlet'
        Ayy(end,:) = 0;
end
Ayy = kron(speye(nx2+1),Ayy);
switch Boundary_left
        case 'dirichlet'
        Ayy(1:nx1+1,:) = 0;
end
switch Boundary_right
        case 'dirichlet'
        Ayy(end-(nx1+1)+1:end,:) = 0;
end

% x direction:
Axx = spdiags([ones((nx1+1)*(nx2+1),1) -2*ones((nx1+1)*(nx2+1),1) ones((nx1+1)*(nx2+1),1)],...
    [-(nx1+1) 0 (nx1+1)],(nx1+1)*(nx2+1),(nx1+1)*(nx2+1));
switch Boundary_left
    case 'noflux'
        Axx(1:nx1+1,nx1+1+1:2*(nx1+1)) = 2*Axx(1:nx1+1,nx1+1+1:2*(nx1+1)); % no flux
    case 'pbc'
        Axx(1:nx1+1,end-(nx1+1)+1:end) = speye(nx1+1,nx1+1); % pbc
    case 'dirichlet'
        Axx(1:nx1+1,:) = 0;
end
switch Boundary_right
    case 'noflux'
        Axx(end-(nx1+1)+1:end,end-2*(nx1+1)+1:end-(nx1+1)) = ...
            2*Axx(end-(nx1+1)+1:end,end-2*(nx1+1)+1:end-(nx1+1)); % no flux
    case 'pbc'
        Axx(end-(nx1+1)+1:end,1:nx1+1) = speye(nx1+1,nx1+1); % pbc
    case 'dirichlet'
        Axx(end-(nx1+1)+1:end,:) = 0;
end
switch Boundary_up
    case 'dirichlet'
        Axx(1:nx1+1:end,:) = 0;
end
switch Boundary_down
    case 'dirichlet'
        Axx(nx1+1:nx1+1:end,:) = 0;
end

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
simbio2ode(Mobj,'generated_model','range');

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

% followed species names
followed_names = Ecoli.data.StatesToLog;
followed_names{end+1} = 'C';
% serial number of the followed species
followed = [find(contains({Mobj.Species.Name},Ecoli.data.StatesToLog)),n_species];


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
    if nplot==length(followed)
        h_plot{nplot} = nextsurf(reshape(c{followed(nplot)},[nx1+1,nx2+1]),L1,L2,nx1,nx2,['[' followed_names{nplot} ']']);
    else
        h_plot{nplot} = nextsurf(reshape(c{followed(length(followed))}.*c{followed(nplot)},[nx1+1,nx2+1]),L1,L2,nx1,nx2,['[C] [' followed_names{nplot} ']']);
    end
end
drawnow

tic
for i = 1:nt

    % active zone
    R = c{end}>cmin & c{end}<=cmax;

    % === Diffusion ===
    % D*dt*Nabla(cell)/cell
    Ncellx = D(end)*dt*(Ax*c{end})./(c{end}+1e-6);
    Ncelly = D(end)*dt*(Ay*c{end})./(c{end}+1e-6);
    % % just when it is > cmin
    % Ncellx = Ncellx.*(c{end}>cmin);
    % Ncelly = Ncelly.*(c{end}>cmin);

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
        else
            % the non diffusible species are in the cell and move with it
            c{j} = c{j}+(Ax*c{j}).*Ncellx+(Ay*c{j}).*Ncelly;
        end
    end
   

   % === Reactions ===
   dc = generated_model(i*dt, c, p, R);
   for j = 1:n_species
       c{j}(R) = c{j}(R) + dc{j}*dt;
   end

   % plot the data
   if rem(i,nt/n_plot) == 0

       % update the surf
       for nplot = 1:length(followed)
           if nplot==length(followed)
               set(h_plot{nplot},'ZData',reshape(c{followed(nplot)},[nx1+1,nx2+1]));
           else
               set(h_plot{nplot},'ZData',reshape(c{followed(length(followed))}.*c{followed(nplot)},[nx1+1,nx2+1]));
           end
       end
       drawnow

       % hatralevo ido:
       disp([int2str(i/nt*100) '%, a varhato befejezes: ' ...
           num2str(toc/i*(nt-i),'%4.2f')])
    end

end
toc
