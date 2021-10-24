% Load data
materialFactory;
load m1data.mat
load m2data.mat

%% Distinct load points
steps = [25 50 75 100]+1;
nsteps = numel(steps);
[ha, pos] = tight_subplot(2, 4, [0.03, 0.02], [0.05 0.12], [0.05 0.05]);
for k = 1:nsteps
    step = steps(k);
    ax1 = ha(k);
    hold(ax1, 'ON');
    pload = int32((50 - abs(step-51))/50*100);
    title(ax1, sprintf('%i%%', pload));
    ax2 = ha(nsteps + k);
    hold(ax2, 'ON');
    
    % m1
    kstep = kappa1(:, step);
    dofp = find(kstep > 0);
    dofe = (1:nelm)'; dofe(dofp) = [];
    exel = ex(dofe, [1 2 3]); eyel = ey(dofe, [1 2 3]);
    expl = ex(dofp, [1 2 3]); eypl = ey(dofp, [1 2 3]);
    fill(ax1, exel', eyel', 1);
    fill(ax1, expl', eypl', -1);
    colormap(ax1, 'hot');
    xticks(ax1, {});
    yticks(ax1, {});

    % m2
    kstep = kappa2(:, step);
    dofp = find(kstep > 0);
    dofe = (1:nelm)'; dofe(dofp) = [];
    exel = ex(dofe, [1 2 3]); eyel = ey(dofe, [1 2 3]);
    expl = ex(dofp, [1 2 3]); eypl = ey(dofp, [1 2 3]);
    fill(ax2, exel', eyel', 1);
    fill(ax2, expl', eypl', -1);
    colormap(ax2, 'hot');
    xticks(ax2, {});
    yticks(ax2, {});
    axis([ax1, ax2], 'tight');
end
ax1 = ha(1);
ylabel(ax1, '$\textbf{Material 1}$', 'Interpreter', 'Latex', 'FontSize', 12);
ax2 = ha(nsteps+1);
ylabel(ax2, '$\textbf{Material 2}$', 'Interpreter', 'Latex', 'FontSize', 12);
sgtitle('Plastic response at different load levels');

%% Heatmap
kappa1load = kappa1(:, 1:50);
level1 = (1 - sum(kappa1load > 0, 2)/size(kappa1load, 2))*100;
level1avg = enod_average(level1, enod);

kappa2load = kappa2(:, 1:50);
level2 = (1 - sum(kappa2load > 0, 2)/size(kappa2load, 2))*100;
level2avg = enod_average(level2, enod);

figure;
tiledlayout(2, 1);
ax1 = nexttile;
fill(ax1, ex', ey', level1avg', 'Linestyle', 'None');
xticks({});
yticks({});
title("% of max load yielding starts");

ax2 = nexttile;
fill(ax2, ex', ey', level2avg', 'Linestyle', 'None');
xticks({});
yticks({});

cb = colorbar('Ticks', [0 0.25 0.5 0.75 1]*100);
cb.Layout.Tile = 'east';

colormap('Hot');
axis([ax1, ax2], 'tight');