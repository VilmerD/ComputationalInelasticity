function D = my_tangent(eps, mp)
% Material parameters
B = mp(1);
E = mp(2);
nu = mp(3);

% Bulk- and shear modulus
G0 = E/(2*(1+nu));
K = E/(3*(1-2*nu));

% eps = [eps_11 eps_22 2*eps_12]
epskk = eps(1) + eps(2) + 0;
e11 = eps(1) - epskk/3;
e22 = eps(2) - epskk/3;
e33 = 0      - epskk/3;
e12 = eps(3)/2;
Jt = (e11^2 + e22^2 + e33^2 + 2*e12^2)/2;
G = G0*exp(B*Jt);

d = (3*K - 2*G)/3;

D1111 = 2*G*(B*e11*e11 + 1.0) + d;
D1122 = 2*G*(B*e11*e22 + 0.0) + d;
D1112 = 2*G*(B*e11*e12 + 0.0) + 0;
D2222 = 2*G*(B*e22*e22 + 1.0) + d;
D2212 = 2*G*(B*e22*e12 + 0.0) + 0;
D1212 = 2*G*(B*e12*e12 + 0.5) + 0;
D = [D1111 D1122 D1112; D1122 D2222 D2212; D1112 D2212 D1212];
end