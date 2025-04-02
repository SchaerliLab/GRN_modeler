function [Ax,Ay] = grad_matrices(nx1,nx2,b,h)
%GRAD_MATRICES create matrices to calculate gradient
% size of the system: (nx1+1)*(nx2+1)
% b is a structure containing the boundary condition
% b.Boundary_up, b.Boundary_down, b.Boundary_left, b.Boundary_right
% with values: 'noflux','dirichlet' or 'pbc'
% h: spatial step size
% Ax, Ay: derivative in x and y directions

% Nabla: N = Ax + Ay;
% y direction:
Ay = spdiags([-ones(nx1+1,1) ones(nx1+1,1)],[-1 1],nx1+1,nx1+1);
switch b.Boundary_up
    case {'noflux','dirichlet'}
        Ay(1,:) = 0; % no flux
    case 'pbc'
        Ay(1,end) = -1; % pbc
end
switch b.Boundary_down
    case {'noflux','dirichlet'}
        Ay(end,:) = 0; % no flux
    case 'pbc'
        Ay(end,1) = 1; % pbc
end
Ay = kron(speye(nx2+1),Ay);
switch b.Boundary_left
        case 'dirichlet'
        Ay(1:nx1+1,:) = 0;
end
switch b.Boundary_right
        case 'dirichlet'
        Ay(end-(nx1+1)+1:end,:) = 0;
end

% x direction:
Ax = spdiags([-ones((nx1+1)*(nx2+1),1)  ones((nx1+1)*(nx2+1),1)],...
    [-(nx1+1) (nx1+1)],(nx1+1)*(nx2+1),(nx1+1)*(nx2+1));
switch b.Boundary_left
    case {'noflux','dirichlet'}
        Ax(1:nx1+1,:) = 0; % no flux
    case 'pbc'
        Ax(1:nx1+1,end-(nx1+1)+1:end) = -speye(nx1+1,nx1+1); % pbc
end
switch b.Boundary_right
    case {'noflux','dirichlet'}
        Ax(end-(nx1+1)+1:end,:) = 0; % no flux
    case 'pbc'
        Ax(end-(nx1+1)+1:end,1:nx1+1) = speye(nx1+1,nx1+1); % pbc
end
switch b.Boundary_up
    case 'dirichlet'
        Ax(1:nx1+1:end,:) = 0;
end
switch b.Boundary_down
    case 'dirichlet'
        Ax(nx1+1:nx1+1:end,:) = 0;
end


% with our parameters
Ax = Ax/(2*h);
Ay = Ay/(2*h);

end

