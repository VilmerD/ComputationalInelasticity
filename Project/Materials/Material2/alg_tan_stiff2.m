function Dat = alg_tan_stiff2(es, dl, epeff, De, mp)
Dat = alg_tan_stiff(es, dl, De, hardrate2(epeff, mp));
end