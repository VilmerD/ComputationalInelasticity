% Testing Es
load('testDats');

%% Material 1
mp = [Kinf, hpar, sigma_y0];
Dat = alg_tan_stiff1(sigma_1, dlambda_1, ep_eff_1, Dstar, mp);

eD = norm(Dat - Dats_1, 'fro')/norm(Dat, 'fro');

grade = '';
if (eD) < 1e3*eps
    grade = 'Passed';
else
    grade = 'Failed';
end

fprintf(['Material 1: ', grade, '\n']);

%% Material 2
mp = [alpha, npar, sigma_y0];
Dat = alg_tan_stiff2(sigma_2, dlambda_2, ep_eff_2, Dstar, mp);

eD = norm(Dat - Dats_2, 'fro')/norm(Dat, 'fro');

grade = '';
if (eD) < 1e3*eps
    grade = 'Passed';
else
    grade = 'Failed';
end

fprintf(['Material 1: ', grade, '\n']);