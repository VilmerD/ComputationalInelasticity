function K = aKSS(edof, ec, u, ep, mp, I, J)
nelm = size(edof, 1);
npoints = 3;
nne = (2 * npoints)^2;  % Number of entries in stiffness matrix per element
eq = [0; 0];

X = zeros(nelm*nne, 1);
for elm = 1:nelm
    % Extracting dofs, coords, disps
    elmdof = edof(elm, 2:end);
    ed = u(elmdof);
    exk = ec(elm, 1:2:end);
    eyk = ec(elm, 2:2:end);
    
    % Computing strains and then tangent material
    [~, etk] = plants(exk, eyk, ep, 0, ed');
    Dk = my_tangent(etk, mp);
    
    [Ke, ~] = plante(exk, eyk, ep, Dk, eq);
    k0 = ((elm - 1)*nne + 1); ke = (elm*nne);   % Indices of entries
    X(k0:ke) = Ke;
end

K = sparse(I, J, X);