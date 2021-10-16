load TopConstMod_Assignment2021.mat;
nelm = size(edof, 1);
ndof = size(coord, 1) * 2;  
npoints = 3;            % nbr of points in each element

% Finding the element coordinates in reference state
coordx = coord(:, 1);
coordy = coord(:, 2);
ex = coordx(enod);
ey = coordy(enod);

ec = zeros(nelm, 2*npoints);
ec(:, 1:2:end-1) = ex;
ec(:, 2:2:end) = ey;

clear coordx coordy

% Plane stress
ptype = 1;
t = 1e-3;
ep = [ptype t];

% For assembling the sparse stiffness matrix
I = reshape(kron(edof(:, 2:end), ones(2*npoints, 1))', [], 1);
J = reshape(kron(edof(:, 2:end), ones(2*npoints, 1)')', [], 1);

% Boundary condition
np = bc(:, 1);