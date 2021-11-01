% Load data
geometry;
load m1data.mat
load m2data.mat

%% Distinct load points
f = figure;
steps = [25 50 75 100]+1;
nsteps = numel(steps);
[ha, pos] = tight_subplot(2, 4, [0.03, 0.02], [0.05 0.12], [0.05 0.05]);
gc = [1 1 1];
tc = [0 0 0]; 
for k = 1:nsteps
    step = steps(k);
    ax1 = ha(k);
    ax1.Color = gc;
    ax1.YAxis.Color = gc;
    ax1.XAxis.Color = gc;
    hold(ax1, 'ON');
    pload = int32((50 - abs(step-51))/50*100);
    title(ax1, sprintf('%i%%', pload));
    ax2 = ha(nsteps + k);
    ax2.Color = gc;
    ax2.YAxis.Color = gc;
    ax2.XAxis.Color = gc;
    hold(ax2, 'ON');
    
    % m1
    kstep = kappa1(:, step);
    dofp = find(kstep > 0);
    dofe = (1:nelm)'; dofe(dofp) = [];
    exel = ex(dofe, [1 2 3]); eyel = ey(dofe, [1 2 3]);
    expl = ex(dofp, [1 2 3]); eypl = ey(dofp, [1 2 3]);
    fill(ax1, exel', eyel', 1, 'Linestyle', 'None');
    fill(ax1, expl', eypl', 0, 'Linestyle', 'None');
    
    kprev = kappa1(:, step - 1);
    dofpw = find(kstep - kprev > 0);
    expw = ex(dofpw, [1 2 3]); eypw = ey(dofpw, [1 2 3]); 
    fill(ax1, expw', eypw', -1, 'Linestyle', 'None');
    
    colormap(ax1, 'winter');
    xticks(ax1, {});
    yticks(ax1, {});

    % m2
    kstep = kappa2(:, step);
    dofp = find(kstep > 0);
    dofe = (1:nelm)'; dofe(dofp) = [];
    exel = ex(dofe, [1 2 3]); eyel = ey(dofe, [1 2 3]);
    expl = ex(dofp, [1 2 3]); eypl = ey(dofp, [1 2 3]);
    fill(ax2, exel', eyel', 1, 'Linestyle', 'None');
    fill(ax2, expl', eypl', 0, 'Linestyle', 'None');
    
    kprev = kappa2(:, step - 1);
    dofpw = find(kstep - kprev > 0);
    expw = ex(dofpw, [1 2 3]); eypw = ey(dofpw, [1 2 3]); 
    fill(ax2, expw', eypw', -1, 'Linestyle', 'None');
    
    colormap(ax2, 'winter');
    xticks(ax2, {});
    yticks(ax2, {});
    
    axis([ax1, ax2], 'tight');
end
ax1 = ha(1);
ylabel(ax1, '$\textbf{Material 1}$', 'Interpreter', 'Latex', 'FontSize', 12);
ax1.YAxis.Label.Color = tc;
ax2 = ha(nsteps+1);
ylabel(ax2, '$\textbf{Material 2}$', 'Interpreter', 'Latex', 'FontSize', 12);
ax2.YAxis.Label.Color = tc;
sgtitle('Plastic response at different load levels');