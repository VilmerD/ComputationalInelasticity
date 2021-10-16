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
    
    % Solve for dl
    I = findInterval(cond, 0, 1);
    [dl, ~, exitflag] = fzero(cond, I);
    if exitflag ~= 1
        errstruct.message = 'Could not solve plastic problem';
        error(errstruct);
    end
    
    peteff = ep_eff_old + dl;
    es = Ms(dl);
else
    dl = 0;
    
    peteff = ep_eff_old;
    es = est;
end
end

function I = findInterval(f, x1, x20)
s1 = sign(f(x1));

x2 = x20;
s2 = sign(f(x2));
q = 2;
while s1*s2 > 0
    x2 = x2*q;
    s2 = sign(f(x2));
end
I = [x1 x2];
end