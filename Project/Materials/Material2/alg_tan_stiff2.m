function Dat = alg_tan_stiff2(es_old, dl, ep_eff, Dstar, mp)
Dat = alg_tan_stiff(es_old, dl, Dstar, hardrate2(ep_eff, mp));
end