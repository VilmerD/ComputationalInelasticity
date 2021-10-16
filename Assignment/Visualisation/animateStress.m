function F = animateStress(ex, ey, U, eseff, edof)
nf = size(eseff, 2);
f = figure();
ax = nexttile;
F = cell2struct(cell(nf, 2), {'cdata', 'colormap'}, 2);
sform = '';
cax = [min(eseff(:)) max(eseff(:))];
sax = [min(ex(:)) max(ex(:)) min(ey(:)) max(ey(:))];
for i = 1:nf
    
    % Extract effective stresses in nodes
    eseffi = eseff(:, i);
    esfin = eseffi(edof(:, 2:2:end));
    
    s = sprintf(sform);
    
    % Plot using fill function
    if exist('P', 'var')
        % Displace
        ui = U(:, i);
        edxi = ui(edof(:, 2:2:end));
        edyi = ui(edof(:, 3:2:end));
        exi = ex + edxi;
        eyi = ey + edyi;
        
        % Update Data
        setPatchData(P, 'CData', esfin);
        setPatchData(P, 'XData', exi);
        setPatchData(P, 'YData', eyi);
    else
        P = fill(ax, ex', ey', esfin', 'EdgeColor', 'none');
        % Adjusting labels and such
        axis image;
        ax.XAxis.TickLength = [0 0];
        ax.YAxis.TickLength = [0 0];
        xticklabels({});
        yticklabels({});
        
        r = (sax(4) - sax(3))/(sax(2) - sax(1));
        xres = 600;
        yres = xres*r;
        set(f, 'Position', [700 900 xres yres]);
        cbar = colorbar(ax, 'East');
        caxis(ax, cax);
        axis(ax, sax);
        hold(ax, 'ON');
    end
    drawnow;
    
    %     if exist('A', 'var')
    %         A.String = s;
    %     else
    %         A = annotation('textbox', [0.13 0.87 0.1 0.1], 'string', s, ...
    %             'FitBoxToText', 'on');
    %     end
    
    % Save frame and clear figure for next plot
    figure(f);
    frame = getframe(f);
    F(i) = frame;
end
name = 'P1Seff.avi';
v = VideoWriter(name);
v.FrameRate = 4;
open(v);
writeVideo(v, F);
close(v);
end