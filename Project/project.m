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

%% Material 1
mp = [300e6 21e9 230e6];
[P, U] = NRDCp(edof, ec, I, J, uload, np, Dstar, @update_variables1, ...
    @alg_tan_stiff1, mp, ep);