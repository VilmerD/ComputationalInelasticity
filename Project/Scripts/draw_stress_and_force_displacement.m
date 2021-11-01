% Load data
geometry;
load m1data.mat
load m2data.mat

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

% Creating axes
ax(1) = axes('Units', 'Normalized', 'Position', [px py axw, axh]);
px = px + axw + axwg;
ax(2) = axes('Units', 'Normalized', 'Position', [px py axw, axh]);
py = py - axh - axhg;
px = px0;
ax(3) = axes('Units', 'Normalized', 'Position', [px py axw, axh]);
px = px + axw + axwg;
ax(4) = axes('Units', 'Normalized', 'Position', [px py axw, axh]);
xticks(ax, {}); yticks(ax, {});
hold(ax, 'on');

%%% Material 1 %%%
% Finding displacements
uk = U1(:, load_steps);
edk = uk(edof(:, 2:end));
exk = ex + edk(:, 1:2:end);
eyk = ey + edk(:, 2:2:end);

% Plotting
fill(ax(1), ex', ey', [0.45 0.45 0.95], 'Linestyle', 'None');
fill(ax(1), exk', eyk', es01'/sy0, 'Linestyle', 'None');
ylabel(ax(1), '$\textbf{Material 1}$', 'Interpreter', 'Latex', 'FontSize', 12);
title(ax(1), 'Fully loaded', 'Fontsize', 14);

% Finding displacements
uk = U1(:, end);
edk = uk(edof(:, 2:end));
exk = ex + edk(:, 1:2:end);
eyk = ey + edk(:, 2:2:end);

% Plotting
fill(ax(2), ex', ey', [0.45 0.45 0.95], 'Linestyle', 'None');
fill(ax(2), exk', eyk', esm1'/sy0, 'Linestyle', 'None');
title(ax(2), 'Fully unloaded', 'Fontsize', 14);

%%% Material 2 %%%
% Finding displacements
uk = U2(:, load_steps);
edk = uk(edof(:, 2:end));
exk = ex + edk(:, 1:2:end);
eyk = ey + edk(:, 2:2:end);

% Plotting
fill(ax(3), ex', ey', [0.45 0.45 0.95], 'Linestyle', 'None');
fill(ax(3), exk', eyk', es02'/sy0, 'Linestyle', 'None');
ylabel(ax(3), '$\textbf{Material 2}$', 'Interpreter', 'Latex', 'FontSize', 12);

% Finding displacements
uk = U2(:, end);
edk = uk(edof(:, 2:end));
exk = ex + edk(:, 1:2:end);
eyk = ey + edk(:, 2:2:end);
fill(ax(4), ex', ey', [0.45 0.45 0.95], 'Linestyle', 'None');
fill(ax(4), exk', eyk', esm2'/sy0, 'Linestyle', 'None');

% Colorbar
cb = colorbar('Ticks', [1 2 3]);
xcb = axwm + 2*(axw + axwg);
ycb = axhm;
set(cb, 'Position', [xcb ycb cbw cbh]);
cb.Label.String = '$\sigma_{eff}/\sigma_{y0}$';
cb.Label.Interpreter = 'Latex';
cb.FontSize = 13;
cb.Label.FontSize = 16;

% Fixing axes and saving figure
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
title('Force vs Displacement for the steel profile');
axis(ax, 'tight');
legend('Location', 'southeast');
figure_name = '/Figures/fvsd.png';
figure_path = [pwd figure_name];
saveas(f, figure_path, 'png');