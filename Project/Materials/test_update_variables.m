% Testing Es
load('testEs');

%% Material 1
mp = [Kinf, hpar, sigma_y0];
[es, dl, peteff] = update_variables1(sigma_old, ep_eff_old, delta_eps, ...
    Dstar, mp);

ees = norm(es - sigma_1)/norm(es);
edl = abs(dl - dlambda_1)/abs(dl);
epeteff = norm(peteff - ep_eff_1)/norm(peteff);

grade = '';
if (ees + edl + epeteff) < 100*eps
    grade = 'Passed';
else
    grade = 'Failed';
end

fprintf(['Material 1: ', grade, '\n']);

%% Material 2
[es, dl, peteff] = update_variables2(sigma_old, ep_eff_old, delta_eps, ...
    Dstar, [alpha npar sigma_y0]);

ees = norm(es - sigma_2)/norm(es);
edl = abs(dl - dlambda_2)/abs(dl);
epeteff = norm(peteff - ep_eff_2)/norm(peteff);

grade = '';
if (ees + edl + epeteff) < 1e-6
    grade = 'Passed';
else
    grade = 'Failed';
end

fprintf(['Material 2: ', grade, '\n']);