function [es, dl, peteff] = update_variables2(es_old, peteff_old, ...
    deps, De, mp)
sy = @(x) yieldstress2(x, mp);
[es, dl, peteff] = update_variables(es_old, peteff_old, deps, De, sy);
end