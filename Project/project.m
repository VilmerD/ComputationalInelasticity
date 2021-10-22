%% Load project
materialFactory;

% Linear elastic model
E = 190e9;
nu = 0.3;
Dstar = hooke(ptype, E, nu);

%% Loading scheme
load_steps = 50;                            % Number of loading steps
unload_steps = load_steps;              
total_load = -1e-3;                         % End load
load_increment = total_load/load_steps;     % Incremental load

% Final loading scheme
loading = [ones(1, load_steps), -1*ones(1, unload_steps)]*load_increment;
uload = kron(bc(:, 2), loading);

%% Simulations
% Material 1
mp1 = [300e6 21e9 230e6];
[P1, U1, kappa1] = NRDCp(edof, ec, uload, np, Dstar, @update_variables1, ...
    @alg_tan_stiff1, mp1, ep);
save('m1data', 'P1', 'U1', 'kappa1');

% Material 2
mp2 = [17 0.61 230e6];
[P2, U2, kappa2] = NRDCp(edof, ec, uload, np, Dstar, @update_variables2, ...
    @alg_tan_stiff2, mp2, ep);
save('m2data', 'P2', 'U2', 'kappa2');