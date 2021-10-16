function [P, U] = NRFCSS(K, f, fload, bc)
% NRFCSS loads the structure and computes the displacement with a
% Newton-Raphson scheme
rtol = 1e-6;

% Initializing quantities
[ndof, nsteps] = size(fload);
P = zeros(ndof, nsteps);
U = zeros(ndof, nsteps);

% Free/Prescribed nodes
np = bc(:, 1);
nf = 1:ndof;
nf(np) = [];

uk = U(:, 1);
for k = 1:nsteps
    % Taking initial step
    fk = fload(:, k);
    fint = f(uk);
    res = fint - fk;
    Kt = K(uk);
    du = solveq(Kt, -res, bc);
    
    % Updating and computing residual
    uk = uk + du;
    fint = f(uk);
    res = (fint - fk);
    
    % Correcting
    while norm(res(nf)) > rtol*norm(res(np))
        Kt = K(uk);
        du = solveq(Kt, -res, bc);
        uk = uk + du;
        fint = f(uk);
        res = (fint - fk);
    end
    
    P(:, k) = fint;
    U(:, k) = uk;
end
end