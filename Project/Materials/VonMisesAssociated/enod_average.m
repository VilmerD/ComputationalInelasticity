function qavg = enod_average(q_elm, enod)
% ENOD_AVERAGE(q_elm, enod) averages an elemental quantity over the nodes
% and is used for plotting with the fill function.
nnod = max(enod(:));
nod_avg = zeros(nnod, 1);
for nod = 1:nnod
    [c0, ~] = find(enod == nod);
    nod_avg(nod, 1) = sum(q_elm(c0))/size(c0, 1);
end
qavg = nod_avg(enod);