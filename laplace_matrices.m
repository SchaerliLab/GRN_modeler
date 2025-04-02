function [Axx,Ayy] = laplace_matrices(nx1,nx2,b,h)
% Laplace: L = Axx + Ayy;
% size of the system: (nx1+1)*(nx2+1)
% b is a structure containing the boundary condition
% b.Boundary_up, b.Boundary_down, b.Boundary_left, b.Boundary_right
% with values: 'noflux','dirichlet' or 'pbc'
% h: spatial step size
% Ax, Ay: derivative in x and y directions

% y direction:
Ayy = spdiags([ones(nx1+1,1) -2*ones(nx1+1,1) ones(nx1+1,1)],[-1 0 1],nx1+1,nx1+1);
switch b.Boundary_up
    case 'noflux'
        Ayy(1,2) = 2; % no flux
    case 'pbc'
        Ayy(1,end) = 1; % pbc
    case 'dirichlet'
        Ayy(1,:) = 0;
end
switch b.Boundary_down
    case 'noflux'
        Ayy(end,end-1) = 2; % no flux
    case 'pbc'
        Ayy(end,1) = 1; % pbc
    case 'dirichlet'
        Ayy(end,:) = 0;
end
Ayy = kron(speye(nx2+1),Ayy);
switch b.Boundary_left
        case 'dirichlet'
        Ayy(1:nx1+1,:) = 0;
end
switch b.Boundary_right
        case 'dirichlet'
        Ayy(end-(nx1+1)+1:end,:) = 0;
end

% x direction:
Axx = spdiags([ones((nx1+1)*(nx2+1),1) -2*ones((nx1+1)*(nx2+1),1) ones((nx1+1)*(nx2+1),1)],...
    [-(nx1+1) 0 (nx1+1)],(nx1+1)*(nx2+1),(nx1+1)*(nx2+1));
switch b.Boundary_left
    case 'noflux'
        Axx(1:nx1+1,nx1+1+1:2*(nx1+1)) = 2*Axx(1:nx1+1,nx1+1+1:2*(nx1+1)); % no flux
    case 'pbc'
        Axx(1:nx1+1,end-(nx1+1)+1:end) = speye(nx1+1,nx1+1); % pbc
    case 'dirichlet'
        Axx(1:nx1+1,:) = 0;
end
switch b.Boundary_right
    case 'noflux'
        Axx(end-(nx1+1)+1:end,end-2*(nx1+1)+1:end-(nx1+1)) = ...
            2*Axx(end-(nx1+1)+1:end,end-2*(nx1+1)+1:end-(nx1+1)); % no flux
    case 'pbc'
        Axx(end-(nx1+1)+1:end,1:nx1+1) = speye(nx1+1,nx1+1); % pbc
    case 'dirichlet'
        Axx(end-(nx1+1)+1:end,:) = 0;
end
switch b.Boundary_up
    case 'dirichlet'
        Axx(1:nx1+1:end,:) = 0;
end
switch b.Boundary_down
    case 'dirichlet'
        Axx(nx1+1:nx1+1:end,:) = 0;
end

Axx = Axx/h^2;
Ayy = Ayy/h^2;

end

