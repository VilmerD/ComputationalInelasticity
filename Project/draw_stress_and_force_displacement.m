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

f = figure;
pw = max(ex(:)) - min(ex(:));
ph = max(ey(:)) - min(ey(:));
h = 300;
w = h*pw/ph;
set(f, 'Position', [50 50 w h]);

axhg = 0.04;
axhm = 0.10;
cbw = 0.01;
cblw = 0.02;
cbh = 1 - 2*axhm;
axh = (1 - axhm*2 - axhg)/2;
axwg = 0.01;
axwm = 0.035;
axw = (1 - axwm*2 - axwg - cbw - cblw)/2;
px0 = axwm;
py0 = 1 - axhm - axh;

px = px0;
py = py0;

% Material 1
ax(1) = axes('Units', 'Normalized', 'Position', [px py axw, axh]);
fill(ax(1), ex', ey', es01'/sy0, 'Linestyle', 'None');
xticks({});
yticks({});

ylabel(ax(1), '$\textbf{Material 1}$', 'Interpreter', 'Latex', 'FontSize', 12);
title(ax(1), 'Max Displacement', 'Fontsize', 14);
px = px + axw + axwg;

ax(2) = axes('Units', 'Normalized', 'Position', [px py axw, axh]);
fill(ax(2), ex', ey', esm1'/sy0, 'Linestyle', 'None');
xticks({});
yticks({});
title(ax(2), 'Original state', 'Fontsize', 14);
py = py - axh - axhg;
px = px0;

% Material 2
ax(3) = axes('Units', 'Normalized', 'Position', [px py axw, axh]);
fill(ax(3), ex', ey', es02'/sy0, 'Linestyle', 'None');
xticks({});
yticks({});

ylabel(ax(3), '$\textbf{Material 2}$', 'Interpreter', 'Latex', 'FontSize', 12);
px = px + axw + axwg;

ax(4) = axes('Units', 'Normalized', 'Position', [px py axw, axh]);
fill(ax(4), ex', ey', esm2'/sy0, 'Linestyle', 'None');
xticks({});
yticks({});


cb = colorbar('Ticks', [1 2 3]);
xcb = axwm + 2*(axw + axwg);
ycb = axhm;
set(cb, 'Position', [xcb ycb cbw cbh]);
cb.Label.String = '$\sigma_{eff}/\sigma_{y0}$';
cb.Label.Interpreter = 'Latex';
cb.FontSize = 13;
cb.Label.FontSize = 16;

for axi = ax; caxis(axi, [0 qmax]), axis(axi, 'tight'); end
colormap('hot');
figure_name = '/Figures/m1vsm2full.png';
figure_path = [pwd figure_name];
saveas(f, figure_path, 'png');
%% Plot force-displacement curve
disp_nodes = bc((bc(:, 2) ~= 0), 1);
plot_dof = disp_nodes(1);

f = figure;
ax = nexttile;
hold(ax, 'ON');
plot(ax, -U1(plot_dof, :)*1e3, -P1(plot_dof, :), 'Displayname', 'Material 1');
plot(ax, -U2(plot_dof, :)*1e3, -P2(plot_dof, :), 'Displayname', 'Material 2');
xlabel('Displacement [mm]');
ylabel('Force [N]');
title('Force vs Displacement at load node');
axis(ax, 'tight');
legend('Location', 'southeast');
figure_name = '/Figures/fvsd.png';
figure_path = [pwd figure_name];
saveas(f, figure_path, 'png');