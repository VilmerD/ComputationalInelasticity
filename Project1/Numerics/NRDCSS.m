function [P, U] = NRDCSS(K, f, uload, np, ndof)
% NRFCSS loads the structure and computes the displacement with a
% Newton-Raphson scheme
rtol = 1e-6;

% Initializing quantities
nsteps = size(uload, 2);
P = zeros(ndof, nsteps);
U = zeros(ndof, nsteps);

% Free/Prescribed nodes
nf = 1:ndof;
nf(np) = [];
bcz = [np zeros(size(np))];

uk = U(:, 1);   fk = P(:, 1);
for k = 1:nsteps
    % Taking initial step
    bck = [np uload(:, k)];
    fint = f(uk);
    res = fint - fk;
    Kt = K(uk);
    du = solveq(Kt, -res, bck);
    
    % Updating and computing residual
    uk = uk + du;
    fint = f(uk);
    res = (fint - fk);
    
    % Correcting
    while norm(res(nf)) > rtol*norm(res(np))
        Kt = K(uk);
        du = solveq(Kt, -res, bcz);
        uk = uk + du;
        fint = f(uk);
        res = (fint - fk);
    end
    
    P(:, k) = fint;
    U(:, k) = uk;
end
end