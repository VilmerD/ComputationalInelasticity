% Load data
materialFactory;
load m1data.mat
load m2data.mat

%% Compute von Mises stresses
load_steps = (size(P1, 2)-1)/2;
unload_steps = load_steps;

mp1 = [300e6 21e9 230e6];
es01 = yieldStress(kappa1(:, load_steps), enod, edof, @yieldstress1, mp1);
esm1 = yieldStress(kappa1(:, load_steps + unload_steps), enod, edof, ...
    @yieldstress1, mp1);

mp2 = [17 0.61 230e6];
es02 = yieldStress(kappa2(:, load_steps), enod, edof, @yieldstress2, mp2);
esm2 = yieldStress(kappa2(:, load_steps + unload_steps), enod, edof, ...
    @yieldstress2, mp2);

%% Plot stresses
sy0 = 230e6;
qmax = max(abs([es01(:); esm1(:); es02(:); esm2(:)]/sy0));

figure;
tiledlayout(2, 2);

% Material 1
ax1 = nexttile;
fill(ax1, ex', ey', es01'/sy0, 'Linestyle', 'None');
xticks({});
yticks({});

ylabel(ax1, '$\textbf{Material 1}$', 'Interpreter', 'Latex', 'FontSize', 12);
title(ax1, 'Max Displacement');

ax2 = nexttile;
fill(ax2, ex', ey', esm1'/sy0, 'Linestyle', 'None');
xticks({});
yticks({});
title(ax2, 'Original state');

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

cb = colorbar('Ticks', [1 2 3]);
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