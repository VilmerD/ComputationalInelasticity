%% Square
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
%% Material
ptype = 1;
t = 1e-3;
ep = [ptype, t];
E = 190e9;
nu = 0.3;
De = hooke(ptype, E, nu);

Kinf = 300e6;
h = 21e9;
sy0 = 230e6;
mp = [Kinf, h, sy0];
update_vars = @(es_old, ep_eff_old, delta_eps) update_variables1(es_old, ...
    ep_eff_old, delta_eps, De, mp);
Dats = @(es_old, dl, ep_eff_old) alg_tan_stiff1(es_old, dl, ep_eff_old, ...
    De, mp);

npoints = 3;
I = reshape(kron(edof(:, 2:end), ones(2*npoints, 1))', [], 1);
J = reshape(kron(edof(:, 2:end), ones(2*npoints, 1)')', [], 1);
K = @(Dats) aKSS(edof, ec, Dats, ep, I, J);

f = @(u, update_variables) aFSS(edof, ec, u, update_variables, ep);

dmax = -1e-4;
nsteps = 10;
dincr = dmax/nsteps;
load = [0 ones(1, nsteps) -1*ones(1, nsteps)]*dincr;
uload = kron(bc(:, 2), load);

%% Load
[P, U] = NRDCp(K, f, uload, bc(:, 1), ndof, nelm, update_vars, Dats);
plot(-(U(5, :)), -(P(5, :)));