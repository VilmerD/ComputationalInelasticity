function [es, dl, peteff] = update_variables(es_old, ep_eff_old, ...
    delta_eps, Dstar, sy)

% Compute trial stress, and trial effective stress
est = es_old + Dstar*delta_eps;
P = [2 -1 0; -1 2 0; 0 0 6]/3;
esteff = sqrt(3/2*(est'*P*est));

% Check if the trial stress violates the yield condition
if esteff - sy(ep_eff_old) > 0
    Ms = @(dl) (eye(3) + 3/2*Dstar*P*dl/sy(ep_eff_old + dl))\est;
    cond = @(dl) 3/2*Ms(dl)'*P*Ms(dl) - sy(ep_eff_old + dl)^2;
    dl0 = 0.01;
    dl = fzero(cond, dl0);
    
    peteff = ep_eff_old + dl;
    es = Ms(dl);
else
    dl = 0;
    
    peteff = ep_eff_old;
    es = est;
end
end