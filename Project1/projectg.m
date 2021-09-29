%% Loading data
load('project1data_fine.mat');
ep = [2 3e-2];
npoints = 3;        % Number of points in elements
ndof = numel(p);    % number of degrees of freedom 
mp = [1e5, 10e6, 0.45];

%% Setup Displacement Controlled
control = 1;
u_max = 3e-3;
nsteps = 40;
du = u_max/nsteps;

[bc, ~, edof, ~, coord, enod, plot_dof] = ...
    TopConstMod(p, t, 0, du, ep(2), control);

ec = coord(edof(:, 2:end));
I = reshape(kron(edof(:, 2:end), ones(2*npoints, 1))', [], 1);
J = reshape(kron(edof(:, 2:end), ones(1, 2*npoints))', [], 1);

f = @(u) aFSS(edof, ec, u, ep, mp);
Kt = @(u) aKSS(edof, ec, u, ep, mp, I, J);
%% Computing
uload = kron(bc(:, 2), [ones(1, nsteps) -1*ones(1, nsteps)]);
[P, U] = NRDCSS(Kt, f, uload, bc(:, 1), ndof);
udof = U(plot_dof, :);
Pdof = P(plot_dof, :);

%% Results
f = figure();
ax = nexttile;
hold(ax, 'ON');
plot(ax, udof(:, 1:nsteps), -Pdof(:, 1:nsteps), 'r', ...
    'Displayname', 'Load');
plot(ax, udof(:, (nsteps+1):end), -Pdof(:, (nsteps+1):end), 'g--', ...
    'Displayname', 'Unload');
legend(ax, 'Location', 'northwest');
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