function [es, dl, peteff] = update_variables(sigma_old, ep_eff_old, ...
    delta_eps, Dstar, sy)

% Compute trial stress, and trial effective stress
est = sigma_old + Dstar*delta_eps;
P = [2 -1 0; -1 2 0; 0 0 6]/3;
esteff = sqrt(3/2*(est'*P*est));

% Check if the trial stress violates the yield condition
if esteff - sy(ep_eff_old) > 0
    M = @(dl) inv(eye(3) + 3/2*Dstar*P*dl/sy(ep_eff_old + dl));
    cond = @(dl) 3/2*(M(dl)*est)'*P*(M(dl)*est) - sy(ep_eff_old + dl)^2;
    dl0 = ep_eff_old;
    dl = fzero(cond, dl0);
    
    peteff = ep_eff_old + dl;
    es = M(dl)*est;
else
    dl = 0;
    
    peteff = ep_eff_old;
    es = est;
end
end