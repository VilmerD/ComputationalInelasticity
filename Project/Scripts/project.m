%% Load project
geometry;

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

% Material 2
mp2 = [17 0.61 230e6];
[P2, U2, kappa2] = NRDCp(edof, ec, uload, np, Dstar, @update_variables2, ...
    @alg_tan_stiff2, mp2, ep);

%% Compute von Mises stresses
load_steps = (size(P1, 2)-1)/2;
unload_steps = load_steps;

mp1 = [300e6 21e9 230e6];
es01 = yieldStress(kappa1(:, load_steps), enod, edof, @yieldstress1, mp1);
esm1 = yieldStress(kappa1(:, load_steps + unload_steps), enod, edof, ...
    @yieldstress1, mp1);

mp2 = [17 0.61 230e6];
es02 = yieldStress(kappa2(:, load_steps), enod, edof, @yieldstress2, mp2);
esm2 = yieldStress(kappa2(:, load_steps + unload_steps), enod, edof, ...
    @yieldstress2, mp2);

%% Save data
folder = pwd;
if ispc
    folder = [folder, '\Data\'];
elseif isunix || ismac
    folder = [folder, '/Data/']; 
end
data1_path = [folder, 'm1data'];
data2_path = [folder, 'm2data'];
save(data1_path, 'P1', 'U1', 'kappa1', 'mp1', 'es01', 'esm1');
save(data2_path', 'P2', 'U2', 'kappa2', 'mp2', 'es02', 'esm2');