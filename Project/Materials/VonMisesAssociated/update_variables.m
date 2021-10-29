function [es, dl, ep_eff] = update_variables(es_old, ep_eff_old, ...
    delta_eps, Dstar, sy)
% Compute trial stress, and trial effective stress
est = es_old + Dstar*delta_eps;
P = [2 -1 0; -1 2 0; 0 0 6]/3;
esteff = sqrt(3/2*(est'*P*est));

% Check if the trial stress violates the yield condition
if esteff - sy(ep_eff_old) > 0
    es2 = @(dl) (eye(3) + 3/2*Dstar*P*dl/sy(ep_eff_old + dl))\est;
    cond = @(dl) 3/2*es2(dl)'*P*es2(dl) - sy(ep_eff_old + dl)^2;
    
    % Solve for dl
    I = [0, 0.1];
    dl = fzero(cond, I);
    
    ep_eff = ep_eff_old + dl;
    es = es2(dl);
else
    % If not, the response is elastic
    dl = 0;
    
    ep_eff = ep_eff_old;
    es = est;
end
end