% Load data
materialFactory;
load m1data.mat
load m2data.mat

%%
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