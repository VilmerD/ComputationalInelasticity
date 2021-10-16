function [P, U] = NRDCp(edof, ec, I, J, uload, np, Dstar, update_vars, ...
    Dats, mp, ep)
% NRFCSS loads the structure and computes the displacement with a
% Newton-Raphson scheme
%
% update_vars   is assumed to be a function of
%                   - es_old
%                   - ep_eff_old
%                   - delta_eps
%
%               meaning Dstar and mp are already filled in
%
% Dats          is assumed to be a function of
%                   - es_old
%                   - dl
%                   - ep_eff_old
%               meaning Dstar and np are already filled in

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

% Initializing quantities
nsteps = size(uload, 2);
P = zeros(ndof, nsteps);
U = zeros(ndof, nsteps);

es_old = zeros(3, nelm);
et_old = zeros(3, nelm);
ep_eff_old = zeros(nelm, 1);

es = zeros(3, nelm);
et = zeros(3, nelm);
ep_eff = zeros(nelm, 1);
dl = zeros(nelm, 1);
for k = 1:nsteps
    %%% Taking initial step %%%
    res = zeros(ndof, 1);
    uk = U(:, k);
    fk = P(:, k);
    du_tot = zeros(ndof, 1);
    
    %%% Correcting %%%
    nsteps = 0;
    while nsteps == 0 || norm(res(nf)) > rtol*norm(res(np))
        % Find load
        bci = uload(:, k) - du_tot(np);
        % New algorithmic tangent each correction step
        
        % Stiffness matrix
        X = zeros(nelm*nne, 1);
        for elm = 1:nelm
            % Extracting dofs, coords, disps
            exk = ec(elm, 1:2:end);
            eyk = ec(elm, 2:2:end);
            
            % Computing strains and then tangent material
            Dk = Dats(es_old(:, elm), dl(elm), ep_eff_old(elm), Dstar, mp);
            
            [Ke, ~] = plante(exk, eyk, ep, Dk, eq);
            k0 = ((elm - 1)*nne + 1); ke = (elm*nne);   % Indices of entries
            X(k0:ke) = Ke(:);
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
        for elm = 1:size(edof, 1)
            % Extracting dofs, coords, disps
            elmdof = edof(elm, 2:end);
            
            % Computing strains and then stresses
            [esk, dlk, ep_eff_k] = update_vars(es_old(:, elm), ...
                ep_eff_old(elm), et(:, elm) - et_old(:, elm), Dstar, mp);
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
    P(:, k) = fint;
    U(:, k) = uk;
    es_old = es;
    et_old = et;
    ep_eff_old = ep_eff;
    
    % Display some information
    fprintf('Step %i\n#lineqs: %i\n', k, nsteps);
end
end