function [es, dl, peteff] = update_variables1(es_old, ep_eff_old, ...
    delta_eps, Dstar, mp)
sy = @(x) yieldstress1(x, mp);
[es, dl, peteff] = update_variables(es_old, ep_eff_old, delta_eps, Dstar, sy);
end