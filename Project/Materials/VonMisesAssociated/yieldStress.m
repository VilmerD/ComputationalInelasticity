function es_eff_elm_avg = yieldStress(kappa, enod, edof, sy, mp)
nelm = size(edof, 1);
es_eff_elm = zeros(nelm, 1);
for elm = 1:nelm
    es_eff_elm(elm) = sy(kappa(elm), mp);
end

es_eff_elm_avg = enod_average(es_eff_elm, enod);
end