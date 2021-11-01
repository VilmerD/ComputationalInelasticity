%% Load project
square_geometry;

% Linear elastic model
E = 190e9;
nu = 0.3;
Dstar = hooke(ptype, E, nu);

%% Loading scheme
load_steps = 50;                            % Number of loading steps
unload_steps = load_steps;              
total_load = -1e-3;                         % End load
load_increment = total_load/load_steps;     % Incremental load

% Final loading scheme
loading = kron([1 -1 1], ones(1, load_steps))*load_increment;
uload = kron(bc(:, 2), loading);

%% Material 1
mp1 = [300e6 21e9 230e6];
[P1, U1, kappa1] = NRDCp(edof, ec, uload, np, Dstar, @update_variables1, ...
    @alg_tan_stiff1, mp1, ep);

mp2 = [17 0.61 230e6];
[P2, U2, kappa2] = NRDCp(edof, ec, uload, np, Dstar, @update_variables2, ...
    @alg_tan_stiff2, mp2, ep);

%% Plot force-displacement curve
disp_nodes = bc((bc(:, 2) ~= 0), 1);
plot_dof = disp_nodes(1);

figure;
ax = nexttile;
hold(ax, 'ON');
plot(ax, -U1(plot_dof, :)*1e3, -P1(plot_dof, :), 'Displayname', 'Material 1');
plot(ax, -U2(plot_dof, :)*1e3, -P2(plot_dof, :), 'Displayname', 'Material 2');
xlabel('Displacement [mm]');
ylabel('Force [N]');
title('Force vs Displacement for the square');
axis(ax, 'tight');
legend('Location', 'southeast');