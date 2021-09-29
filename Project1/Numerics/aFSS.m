function fint = aFSS(edof, ec, u, ep, mp)
ndof = numel(u);

fint = zeros(ndof, 1);
for elm = 1:size(edof, 1)
    % Extracting dofs, coords, disps
    elmdof = edof(elm, 2:end);
    ed = u(elmdof);
    exk = ec(elm, 1:2:end);
    eyk = ec(elm, 2:2:end);
    
    % Computing strains and then stresses
    [~, etk] = plants(exk, eyk, ep, 0, ed');
    esk = my_stress(etk, mp);
    
    % Computing elemt forces
    efk = plantf(exk, eyk, ep, esk');
    fint(elmdof) = fint(elmdof) + efk;
end