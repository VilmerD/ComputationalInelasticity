function [es, dl, peteff] = update_variables1(es_old, peteff_old, ...
    deps, De, mp)
sy = @(x) yieldstress1(x, mp);
[es, dl, peteff] = update_variables(es_old, peteff_old, deps, De, sy);
end