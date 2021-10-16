function [P, U] = NRDCp_old(K, f, uload, np, ndof, nelm, update_vars, Dats)
% NRFCSS loads the structure and computes the displacement with a
% Newton-Raphson scheme
% K             is assumed to be a function of
%                   - Dats (function of elm)
% 
%               meaning edof, ec, ep, I and J are already filled in
%
% f             is assumed to be a function of
%                   - u
%                   - update variables (function of et and elm)
% 
%               meaning edof, ec, and ep are already filled in
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

% Boundary condition
nf = 1:ndof;
nf(np) = [];

% Initializing quantities
nsteps = size(uload, 2);
P = zeros(ndof, nsteps);
U = zeros(ndof, nsteps);

es_old = zeros(3, nelm);
et_old = zeros(3, nelm);
ep_eff_old = zeros(nelm, 1);

dl = zeros(nelm, 1);
for k = 1:nsteps
    %%% Taking initial step %%%
    % Update function stays constant in a given loading step as
    % it depends on the previously converges quantities    
    update_vars_k = @(etk, elm) update_vars(...
        es_old(:, elm), ep_eff_old(elm), etk - et_old(:, elm));
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
        Dk = @(k) Dats(es_old(:, k), dl(k), ep_eff_old(k));
        Kt = K(Dk);
        du = solveq(Kt, -res, [np bci]);
        uk = uk + du;
        du_tot = du_tot + du;
        
        % Updating quantities after each step, and checking for convergance
        [fint, et, es, dl, ep_eff] = f(uk, update_vars_k);
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