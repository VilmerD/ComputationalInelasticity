function Dat = alg_tan_stiff(es, dl, Dstar, da)
% If the response is elastic, the elastic material tangent is returned
if (dl == 0)
    Dat = Dstar;
    return
end

% Otherwise the algorithmic material tangent is returned
[dfds, ddfdsds] = fgrad(es);
Da = inv(inv(Dstar) + dl*ddfdsds);
A = dfds'*Da*dfds + da;
Dat = Da - 1/A*Da*(dfds*dfds')*Da;
end

% Helper function that computes the gradient of f wrt the stresses
function [dfds, ddfdsds] = fgrad(es)
P = [2 -1 0; -1 2 0; 0 0 6]/3;
s = P*es;

eseff = sqrt(3/2*es'*s);

dfds = 3/2*P*es/eseff;

ddfdsds = 3/2*(P/eseff - (3/2)*(s*s')/eseff^3);
end