function Dat = alg_tan_stiff1(es, dl, epeff, De, mp)
Dat = alg_tan_stiff(es, dl, De, hardrate1(epeff, mp));
end