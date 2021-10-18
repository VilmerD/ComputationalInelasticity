%% Load project
squareFactory;

% Linear elastic model
E = 190e9;
nu = 0.3;
Dstar = hooke(ptype, E, nu);

%% Loading scheme
load_steps = 50;                            % Number of loading steps
unload_steps = load_steps;              
total_load = -0.5e-3;                         % End load
load_increment = total_load/load_steps;     % Incremental load

% Final loading scheme
loading = [ones(1, load_steps), -1*ones(1, unload_steps)]*load_increment;
uload = kron(bc(:, 2), loading);

%% Material 1
mp1 = [300e6 21e9 230e6];
[P1, U1, kappa1] = NRDCp(edof, ec, uload, np, Dstar, @update_variables1, ...
    @alg_tan_stiff1, mp1, ep);

mp2 = [17 0.61 230e6];
[P2, U2, kappa2] = NRDCp(edof, ec, uload, np, Dstar, @update_variables2, ...
    @alg_tan_stiff2, mp2, ep);

%% Stresses
es01 = yieldStress(kappa1(:, load_steps), enod, edof, @yieldstress1, mp1);
esm1 = yieldStress(kappa1(:, load_steps + unload_steps), enod, edof, ...
@yieldstress1, mp1);

es02 = yieldStress(kappa2(:, load_steps), enod, edof, @yieldstress2, mp2);
esm2 = yieldStress(kappa2(:, load_steps + unload_steps), enod, edof, ...
@yieldstress2, mp2);

%% Plot stresses
sy0 = mp1(3);
qmax = max(abs([es01(:); esm1(:); es02(:); esm2(:)]/sy0));

figure();
tiledlayout(2, 2);

% Material 1
ax1 = nexttile;
fill(ax1, ex', ey', es01'/sy0, 'Linestyle', 'None');
xticks({});
yticks({});

ylabel(ax1, '$\textbf{Material 1}$', 'Interpreter', 'Latex', 'FontSize', 12);
title(ax1, 'Max Load');

ax2 = nexttile;
fill(ax2, ex', ey', esm1'/sy0, 'Linestyle', 'None');
xticks({});
yticks({});
title(ax2, 'Zero Load');

% Material 2
ax3 = nexttile;
fill(ax3, ex', ey', es02'/sy0, 'Linestyle', 'None');
xticks({});
yticks({});

ylabel(ax3, '$\textbf{Material 2}$', 'Interpreter', 'Latex', 'FontSize', 12);

ax4 = nexttile;
fill(ax4, ex', ey', esm2'/sy0, 'Linestyle', 'None');
xticks({});
yticks({});

axis([ax1 ax2 ax3 ax4], 'tight');

cb = colorbar('Ticks', 0:ceil(qmax));
cb.Layout.Tile = 'east';
cb.Label.String = '$\sigma_{eff}/\sigma_{y0}$';
cb.Label.Interpreter = 'Latex';
cb.Label.FontSize = 14;

caxis(ax1, [0 qmax]);
caxis(ax2, [0 qmax]);
caxis(ax3, [0 qmax]);
caxis(ax4, [0 qmax]);

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
title('Force vs Displacement at load node');
axis(ax, 'tight');
legend('Location', 'southeast');