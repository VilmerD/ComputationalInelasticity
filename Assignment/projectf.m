%% Loading data
load('project1data_fine.mat');
ep = [2 3e-2];
npoints = 3;        % Number of points in elements
ndof = numel(p);    % number of degrees of freedom
mp = [1e5, 10e6, 0.45];

%% Setup
control = 0;
tau_max = 30e3;
nsteps = 30;
dtau = tau_max/nsteps;

[bc, df, edof, ~, coord, enod, plot_dof] = ...
    TopConstMod(p, t, dtau, 0, ep(2), control);

ec = coord(edof(:, 2:end));
% Makes assembling the (sparse) stiffness matrix faster
I = reshape(kron(edof(:, 2:end), ones(2*npoints, 1))', [], 1);
J = reshape(kron(edof(:, 2:end), ones(2*npoints, 1)')', [], 1);

f = @(u) aFSS(edof, ec, u, ep, mp);
Kt = @(u) aKSS(edof, ec, u, ep, mp, I, J);
%% Computing
fload = kron(df, 0:nsteps);
[P, U] = NRFCSS(Kt, f, fload, bc);
udof = U(plot_dof, :);
Pdof = P(plot_dof, :);

%% Results
f = figure();
ax = nexttile;
hold(ax, 'ON');
plot(udof, Pdof);
xstr = sprintf('$u_{%i}$', plot_dof);
ystr = sprintf('$f_{%i}$', plot_dof);
xlabel(ax, xstr, 'Interpreter', 'Latex');
ylabel(ax, ystr, 'Interpreter', 'Latex');

%% VonMises
es = stresses(edof, ec, U, ep, mp);

eseff = zeros(ndof/2, size(U, 2));
for k = 1:size(U, 2)
    esk = es(:, k);
    for node = 1:ndof/2
        [c0,c1] = find(enod==node);
        eseff(node, k) = sum(esk(c0)/size(c0,1));
    end
end

% Plot
f = figure();
ax = nexttile;
hold(ax, 'ON');

k = nsteps;
eseffk = eseff(:, k);
eseffknod = eseffk(edof(:, 2:2:end));
uk = U(:, k);
eck = ec + uk(edof(:, 2:end));
exk = eck(:, 1:2:end);
eyk = eck(:, 2:2:end);

fill(ax, exk', eyk', eseffknod', 'Linestyle', 'none');

xticks(ax, {});
yticks(ax, {});
colorbar();