ndof = 8;
nelm = 2;
edof = [1 1 2 3 4 5 6
        2 1 2 5 6 7 8];
enod = [1 1 2 3;
        2 1 3 4];
l = 10e-3;
coord = [0 0
         l 0
         l l
         0 l];
cx = coord(:, 1);
cy = coord(:, 2);
ex = cx(enod(:, 2:end));
ey = cy(enod(:, 2:end));
ec = zeros(nelm, 2*3);
ec(:, 1:2:end) = ex;
ec(:, 2:2:end) = ey;

bc = [1 0
      2 0
      3 1
      4 0
      5 1
      7 0];

clear coordx coordy

% Plane stress
ptype = 1;
t = 1e-3;
ep = [ptype t];

% Boundary condition
np = bc(:, 1);