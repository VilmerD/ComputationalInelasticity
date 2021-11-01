function [P, U, kappa] = NRDCp(edof, ec, uload, np, Dstar, update_vars, ...
    Dats, mp, ep)
% NRFCSS loads the structure and computes the displacement with a
% Newton-Raphson scheme

% Relative tolerance
rtol = 1e-6;

% Shape properties
nelm = size(edof, 1);
ndof = max(edof(:));
ex = ec(:, 1:2:end-1);
ey = ec(:, 2:2:end);

% Boundary condition
nf = 1:ndof;
nf(np) = [];

% Stiffness matrix stuff
npoints = 3;
nne = (2 * npoints)^2;  % Number of entries in stiffness matrix per element
eq = [0; 0];
I = reshape(kron(edof(:, 2:end), ones(2*npoints, 1))', [], 1);
J = reshape(kron(edof(:, 2:end), ones(2*npoints, 1)')', [], 1);

% Initializing quantities
nsteps = size(uload, 2);
P = zeros(ndof, nsteps);
U = zeros(ndof, nsteps);
kappa = zeros(nelm, nsteps);

es_old = zeros(3, nelm);
et_old = zeros(3, nelm);
ep_eff_old = zeros(nelm, 1);
for k = 1:nsteps
    %%% Initializing load quantities %%%
    res = zeros(ndof, 1);
    uk = U(:, k);
    fk = zeros(ndof, 1);
    du_tot = zeros(ndof, 1);
    dl = zeros(nelm, 1);
    es = zeros(3, nelm);
    ep_eff = zeros(nelm, 1);
    
    %%% Loading and correcting %%%
    nsteps = 0;
    while nsteps == 0 || norm(res(nf)) > rtol*norm(res(np))
        % Find load
        bci = uload(:, k) - du_tot(np);
        
        % Stiffness matrix
        X = zeros(nne*nelm, 1);
        for elm = 1:nelm
            % Extracting dofs, coords, disps
            exk = ec(elm, 1:2:end);
            eyk = ec(elm, 2:2:end);
            
            % Computing strains and then tangent material
            Dk = Dats(es(:, elm), dl(elm), ep_eff(elm), Dstar, mp);
            
            [Ke, ~] = plante(exk, eyk, ep, Dk, eq);
            k0 = 1 + (elm - 1)*nne;
            ke = elm * nne;
            X(k0:ke) = Ke;
        end
        K = sparse(I, J, X);
        
        % Solving lienar system
        du = solveq(K, -res, [np bci]);
        uk = uk + du;
        du_tot = du_tot + du;
        
        % Updating quantities after each step, and checking for convergance
        fint = zeros(ndof, 1);
        [~, et] = plants(ex, ey, ep, 0, uk(edof(:, 2:end)));
        et = et';
        for elm = 1:nelm
            % Extracting dofs, coords, disps
            elmdof = edof(elm, 2:end);
            exk = ec(elm, 1:2:end);
            eyk = ec(elm, 2:2:end);
            
            % Computing strains and then stresses
            delta_et = et(:, elm) - et_old(:, elm);
            [esk, dlk, ep_eff_k] = update_vars(es_old(:, elm), ...
                ep_eff_old(elm), delta_et, Dstar, mp);
            es(:, elm) = esk;
            dl(elm) = dlk;
            ep_eff(elm) = ep_eff_k;
            
            % Computing element forces
            efk = plantf(exk, eyk, ep, esk');
            fint(elmdof) = fint(elmdof) + efk;
        end
        res = fint - fk;
        nsteps = nsteps + 1;
    end
    
    %%% Finnally updating converged quantites %%%
    P(:, k+1) = fint;
    U(:, k+1) = uk;
    kappa(:, k+1) = ep_eff;
    es_old = es;
    et_old = et;
    ep_eff_old = ep_eff;
    
    % Display some information
    fprintf('Step %i\n#lineqs: %i\n', k, nsteps);
end
end