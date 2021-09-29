function es = stresses(edof, ec, u, ep, mp)
[ndof, nsteps] = size(u);

es = zeros(ndof, nsteps);
for k = 1:nsteps
    uk = u(:, k);
    for elm = 1:size(edof, 1)
        % Extracting dofs, coords, disps
        elmdof = edof(elm, 2:end);
        ed = uk(elmdof);
        exk = ec(elm, 1:2:end);
        eyk = ec(elm, 2:2:end);
        
        % Computing strains and then stresses
        [~, etk] = plants(exk, eyk, ep, 0, ed');
        esk = my_stress(etk, mp);
        
        esv = vonMises(esk);
        es(elm, k) = esv;
    end
end
end