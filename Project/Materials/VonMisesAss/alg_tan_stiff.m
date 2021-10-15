function Dat = alg_tan_stiff(es, dl, De, da)
[dfds, ddfdsds] = fgrad(es);
Da = inv(inv(De) + dl*ddfdsds);
H = -1*da;
A = dfds'*Da*dfds - H;
Dat = Da - 1/A*Da*(dfds*dfds')*Da;
end

function [dfds, ddfdsds] = fgrad(es)
P = [2 -1 0; -1 2 0; 0 0 6]/3;
s = P*es;

eseff = sqrt(3/2*es'*s);

dfds = 3/2*P*es/eseff;

ddfdsds = 3/2*(P/eseff - (3/2)*(s*s')/eseff^3);
end