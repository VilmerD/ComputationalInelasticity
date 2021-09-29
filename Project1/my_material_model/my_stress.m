function sig = my_stress(eps, mp)
% Defining parameters
B = mp(1);
E = mp(2);
nu = mp(3);

% Bulk- and shear modulus
G0 = E/(2*(1+nu));
K = E/(3*(1-2*nu));

% eps = [eps_11 eps_22 2*eps_12]
epskk = eps(1) + eps(2) + 0;
Jt = ((eps(1) - epskk/3)^2 + (eps(2) - epskk/3)^2 + ...
    (0-epskk/3)^2 + 2*(eps(3)/2)^2)/2;
G = G0*exp(B*Jt);

% sig = [sig_11 sig_22 sig_33 sig_12]
ds = epskk*(3*K - 2*G)/3;
sig = [2*G*eps(1) + ds; 2*G*eps(2) + ds; ds; 2*G*eps(3)/2];
end