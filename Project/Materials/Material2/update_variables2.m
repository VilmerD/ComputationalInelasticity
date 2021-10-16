function [es, dl, ep_eff_k] = update_variables2(es_old, ep_eff_old, ...
    deps, De, mp)
sy = @(x) yieldstress2(x, mp);
[es, dl, ep_eff_k] = update_variables(es_old, ep_eff_old, deps, De, sy);
end