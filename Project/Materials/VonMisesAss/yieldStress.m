function es_eff_elm_avg = yieldStress(kappa, enod, edof, sy, mp)
nelm = size(edof, 1);
es_eff_nod = zeros(nelm, 1);
for elm = 1:nelm
    es_eff_nod(elm) = sy(kappa(elm), mp);
end

ndof = max(edof(:));
es_eff_nod_avg = zeros(ndof/2, 1);
for nod = 1:ndof/2
    [c0, ~] = find(enod == nod);
    es_eff_nod_avg(nod, 1) = sum(es_eff_nod(c0))/size(c0, 1);
end

es_eff_elm_avg = es_eff_nod_avg(enod(:, 2:end));
end